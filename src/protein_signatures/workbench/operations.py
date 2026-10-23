"""Named workbench operations; all analytical outputs belong to new jobs."""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from protein_signatures.errors import InputValidationError
from protein_signatures.io_utils import write_json_atomic, write_tsv_atomic

from .analysis import align_family, index_fasta, run_command, search_sequences, write_fasta
from .catalogue import open_catalogue, resolve_ids, retrieve_clusters, source_view
from .config import WorkbenchConfig, bounded_number
from .jobs import resolve_input

LOGGER = logging.getLogger(__name__)


def execute_job(*, config: WorkbenchConfig, request: dict, output: Path) -> None:
    """Dispatch one validated named operation with explicit allowed parameters.

    Args:
        config: Trusted worker profile.
        request: Saved job request and input references.
        output: New private job output directory.

    Raises:
        InputValidationError: If a request is unsupported or malformed.
    """
    operation = request["operation"]
    parameters = request["parameters"]
    allowed = {
        "resolve": {"identifiers"},
        "retrieve": {"cluster_ids"},
        "search": {"identity", "coverage", "evalue", "maximum_hits"},
        "screen": {"identity", "coverage", "evalue", "deduplicate", "taxon_ids"},
        "align": {"selected_ids", "tree"},
        "annotate": set(),
        "score": {"comparison_ids", "include_within_comparison"},
        "prepare": set(),
        "build_search": set(),
    }
    if (
        operation not in allowed
        or not isinstance(parameters, dict)
        or set(parameters) - allowed[operation]
    ):
        raise InputValidationError("Unsupported operation or parameter in the job request.")
    Path(config.workspace, "temporary").mkdir(parents=True, exist_ok=True, mode=0o700)
    if operation == "resolve":
        frame = resolve_ids(config=config, values=parameters.get("identifiers", []))
        frame.to_csv(path_or_buf=output / "mapping.tsv", sep="\t", index=False)
    elif operation == "retrieve":
        retrieve_clusters(config=config, values=parameters.get("cluster_ids", []), output=output)
    elif operation == "build_search":
        build_search_database(config=config, output=output)
    elif operation == "search":
        search_representatives(
            config=config,
            source=resolve_input(config=config, request=request, name="seeds"),
            parameters=parameters,
            output=output,
        )
    elif operation == "screen":
        screen_family(
            config=config,
            source=resolve_input(config=config, request=request, name="family"),
            seeds=resolve_input(config=config, request=request, name="seeds"),
            parameters=parameters,
            output=output,
        )
    elif operation == "align":
        align_family(
            config=config,
            source=resolve_input(config=config, request=request, name="family"),
            parameters=parameters,
            output=output,
        )
    elif operation == "annotate":
        annotate_family(
            config=config,
            source=resolve_input(config=config, request=request, name="family"),
            output=output,
        )
    elif operation == "score":
        feature_path = (
            resolve_input(config=config, request=request, name="features")
            if "features" in request["sources"]
            else None
        )
        score_family(
            config=config,
            source=resolve_input(config=config, request=request, name="family"),
            features=feature_path,
            parameters=parameters,
            output=output,
        )
    elif operation == "prepare":
        prepare_campaign(
            config=config,
            source=resolve_input(config=config, request=request, name="family"),
            seeds=resolve_input(config=config, request=request, name="seeds"),
            output=output,
        )
    LOGGER.info("Finished workbench operation %s", operation)


def build_search_database(*, config: WorkbenchConfig, output: Path) -> None:
    """Build a DIAMOND database from an operator-supplied global representative FASTA.

    Args:
        config: Readable centroid_fasta and installed DIAMOND.
        output: Owned job directory on the data disk.

    Raises:
        InputValidationError: If the source is unavailable or output is missing.
    """
    if not config.centroid_fasta or not Path(config.centroid_fasta).is_file():
        raise InputValidationError(
            (
                "Configure centroid_fasta with the release's global representative "
                "FASTA. The sequence Parquet is not a DIAMOND database."
            )
        )
    run_command(
        config=config,
        output=output,
        arguments=[
            config.diamond,
            "makedb",
            "--in",
            config.centroid_fasta,
            "--db",
            str(output / "representatives"),
            "--threads",
            str(config.threads),
        ],
    )
    destination = output / "representatives.dmnd"
    if not destination.is_file() or destination.stat().st_size == 0:
        raise InputValidationError("DIAMOND did not create a representative search database.")
    (output / "README.txt").write_text(
        data=(
            f'Set centroid_database = "{destination}" in your workbench TOML '
            "before submitting new searches.\n"
            "The source FASTA must use the release's exact global centroid IDs.\n"
        ),
        encoding="utf-8",
    )


def search_representatives(
    *, config: WorkbenchConfig, source: Path, parameters: dict, output: Path
) -> None:
    """Find bounded seed-to-centroid similarities, without claiming exact membership.

    Args:
        config: Local global-representative DIAMOND database.
        source: Seed protein FASTA.
        parameters: Search controls.
        output: Job output directory.
    """
    if not config.centroid_database or not Path(config.centroid_database).is_file():
        raise InputValidationError(
            (
                "Global homology search needs centroid_database (.dmnd). Configure the "
                "existing representative database, or build it from centroid_fasta in "
                "Start here."
            )
        )
    maximum = bounded_number(
        value=parameters.get("maximum_hits", 50),
        name="Targets per seed",
        minimum=1,
        maximum=1000,
        integer=True,
    )
    with open_catalogue(config=config, with_index=False) as connection:
        index_fasta(
            connection=connection, path=source, name="seeds", maximum=config.max_seed_sequences
        )
        write_fasta(
            connection=connection,
            query="SELECT sequence_id,sequence FROM seeds ORDER BY sequence_id",
            path=output / "seeds.fasta",
        )
    hits = search_sequences(
        config=config,
        query=output / "seeds.fasta",
        database=config.centroid_database,
        parameters=parameters,
        output=output,
        maximum_hits=maximum,
    )
    with open_catalogue(config=config) as connection:
        source_view(connection=connection, path=str(hits), name="hits")
        connection.execute(
            "CREATE VIEW requested AS SELECT DISTINCT sseqid AS cluster_id FROM hits"
        )
        result = connection.sql(
            query=(
                "SELECT h.*,h.sseqid AS "
                "cluster_id,c.member_count,'SIMILARITY_TO_CENTROID' AS "
                "mapping_evidence,CASE WHEN c.cluster_id IS NULL THEN 'HIT_UNINDEXED' "
                "ELSE 'HIT_INDEXED' END AS status FROM hits h LEFT JOIN clusters c ON "
                "h.sseqid=c.cluster_id ORDER BY h.qseqid,CAST(h.bitscore AS DOUBLE) "
                "DESC,h.sseqid"
            )
        )
        result.write_csv(file_name=str(output / "representative_hits.tsv"), sep="\t", header=True)
    write_json_atomic(
        path=output / "search_scope.json",
        value={
            "maximum_targets_per_seed": maximum,
            "exhaustive": False,
            "interpretation": (
                "Retained sequence-similarity hits; no-hit is not evidence of absence "
                "from the tree of life. Homologues are not automatically orthologues or"
                " functional positives."
            ),
        },
    )


def screen_family(
    *, config: WorkbenchConfig, source: Path, seeds: Path, parameters: dict, output: Path
) -> None:
    """Recheck every retrieved member against seeds with a complete selection audit.

    Args:
        config: DIAMOND installation and input bounds.
        source: Retrieved unaligned family FASTA.
        seeds: Original query FASTA.
        parameters: Similarity, taxonomy and exact-sequence deduplication controls.
        output: Owned job directory.
    """
    deduplicate = parameters.get("deduplicate", True)
    if not isinstance(deduplicate, bool):
        raise InputValidationError("deduplicate must be a Boolean.")
    taxa = parameters.get("taxon_ids", [])
    if not isinstance(taxa, list):
        raise InputValidationError("taxon_ids must be a list of exact NCBI taxonomy identifiers.")
    taxa = [
        bounded_number(
            value=value, name="Taxon ID", minimum=1, maximum=10_000_000_000, integer=True
        )
        for value in taxa
    ]
    if taxa and not config.taxonomy:
        raise InputValidationError("Taxonomy filtering needs a reviewed sequence taxonomy source.")
    with open_catalogue(config=config, with_index=False) as connection:
        total = index_fasta(
            connection=connection, path=source, name="family", maximum=config.max_sequences
        )
        index_fasta(
            connection=connection, path=seeds, name="seeds", maximum=config.max_seed_sequences
        )
        write_fasta(
            connection=connection,
            query="SELECT sequence_id,sequence FROM family ORDER BY sequence_id",
            path=output / "family_input.fasta",
        )
        write_fasta(
            connection=connection,
            query="SELECT sequence_id,sequence FROM seeds ORDER BY sequence_id",
            path=output / "seeds.fasta",
        )
        run_command(
            config=config,
            output=output,
            arguments=[
                config.diamond,
                "makedb",
                "--in",
                str(output / "seeds.fasta"),
                "--db",
                str(output / "seed_database"),
            ],
        )
        hits = search_sequences(
            config=config,
            query=output / "family_input.fasta",
            database=str(output / "seed_database.dmnd"),
            parameters=parameters,
            output=output,
            maximum_hits=1,
        )
        source_view(connection=connection, path=str(hits), name="hits")
        if (
            connection.execute(
                "SELECT qseqid FROM hits EXCEPT SELECT sequence_id FROM family"
            ).fetchone()
            or connection.execute(
                "SELECT sseqid FROM hits EXCEPT SELECT sequence_id FROM seeds"
            ).fetchone()
        ):
            raise InputValidationError(
                "DIAMOND output refers to a protein outside the submitted inputs."
            )
        predicate = "TRUE"
        if taxa:
            columns = source_view(connection=connection, path=config.taxonomy, name="taxonomy")
            reviewed = "AND mapping_status='REVIEWED'" if "mapping_status" in columns else ""
            connection.execute(
                (
                    "CREATE TABLE taxon_filter AS SELECT sequence_id FROM taxonomy GROUP BY"
                    " sequence_id HAVING count(DISTINCT taxon_id)=1 AND "
                    "bool_and(TRY_CAST(taxon_id AS BIGINT) IN (SELECT unnest(?))) "
                )
                + reviewed.replace(
                    "mapping_status='REVIEWED'", "bool_and(mapping_status='REVIEWED')"
                ),
                [taxa],
            )
            predicate = "f.sequence_id IN (SELECT sequence_id FROM taxon_filter)"
        connection.execute(
            (
                "CREATE TABLE passing AS SELECT f.*,h.sseqid AS "
                "best_seed_id,CAST(h.pident AS DOUBLE) AS "
                "identity_percent,CAST(h.query_coverage_percent AS DOUBLE) AS "
                "member_coverage_percent,CAST(h.subject_coverage_percent AS DOUBLE) AS "
                "seed_coverage_percent,CAST(h.evalue AS DOUBLE) AS evalue,row_number() "
                "OVER(PARTITION BY f.sequence_sha256 ORDER BY f.sequence_id) AS "
                "duplicate_rank FROM family f JOIN hits h ON f.sequence_id=h.qseqid "
                "WHERE "
            )
            + predicate
            + (
                " QUALIFY row_number() OVER(PARTITION BY f.sequence_id ORDER BY "
                "CAST(h.bitscore AS DOUBLE) DESC,h.sseqid)=1"
            )
        )
        keep = "duplicate_rank=1" if deduplicate else "TRUE"
        connection.execute("CREATE VIEW selected AS SELECT * FROM passing WHERE " + keep)
        count = write_fasta(
            connection=connection,
            query="SELECT sequence_id,sequence FROM selected ORDER BY sequence_id",
            path=output / "candidates.fasta",
        )
        audit_query = (
            "SELECT "
            "f.sequence_id,f.sequence_sha256,p.best_seed_id,p.identity_percent,"
            "p.member_coverage_percent,p.seed_coverage_percent,p.evalue,CASE"
            " WHEN s.sequence_id IS NOT NULL THEN 'SELECTED' WHEN p.sequence_id IS "
            "NOT NULL THEN 'EXACT_SEQUENCE_DUPLICATE' WHEN h.qseqid IS NOT NULL "
            "THEN 'TAXONOMY_FILTERED_OR_UNKNOWN' ELSE 'NO_RETAINED_HIT' END AS "
            "selection_status FROM family f LEFT JOIN passing p USING(sequence_id) "
            "LEFT JOIN selected s USING(sequence_id) LEFT JOIN (SELECT DISTINCT "
            "qseqid FROM hits) h ON f.sequence_id=h.qseqid ORDER BY f.sequence_id"
        )
        connection.sql(query=audit_query).write_csv(
            file_name=str(output / "candidate_audit.tsv"), sep="\t", header=True
        )
        write_json_atomic(
            path=output / "selection_summary.json",
            value={
                "input_members": total,
                "selected_members": count,
                "exact_sequence_deduplication": deduplicate,
                "exact_taxon_ids": taxa,
                "scope": "Sequence-similarity candidates; functional labels not inferred",
            },
        )


def annotate_family(*, config: WorkbenchConfig, source: Path, output: Path) -> None:
    """Scan a bounded family against a local pressed Pfam HMM database.

    Args:
        config: hmmscan executable and pressed Pfam-A.hmm.
        source: Selected protein FASTA.
        output: Owned job output.

    Raises:
        InputValidationError: If database, IDs or HMMER evidence are invalid.
    """
    if not config.pfam_hmm or not Path(config.pfam_hmm).is_file():
        raise InputValidationError(
            "Configure a local Pfam-A.hmm database, prepared using hmmpress."
        )
    with open_catalogue(config=config, with_index=False) as connection:
        index_fasta(
            connection=connection, path=source, name="family", maximum=config.max_seed_sequences
        )
        write_fasta(
            connection=connection,
            query="SELECT sequence_id,sequence FROM family ORDER BY sequence_id",
            path=output / "sequences.fasta",
        )
        sizes = dict(
            connection.execute("SELECT sequence_id,length(sequence) FROM family").fetchall()
        )
    run_command(
        config=config,
        output=output,
        arguments=[
            config.hmmscan,
            "--cpu",
            str(config.threads),
            "--cut_ga",
            "--domtblout",
            str(output / "pfam.domtblout"),
            config.pfam_hmm,
            str(output / "sequences.fasta"),
        ],
    )
    domain_rows, feature_rows, seen = [], [], set()
    with (output / "pfam.domtblout").open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split(maxsplit=22)
            if len(fields) < 22:
                raise InputValidationError("HMMER returned a malformed domain row.")
            identifier, domain = fields[3], fields[1].split(".")[0]
            start, end = int(fields[17]), int(fields[18])
            if (
                identifier not in sizes
                or not 1 <= start <= end <= sizes[identifier]
                or not domain.startswith("PF")
            ):
                raise InputValidationError(
                    "HMMER returned an invalid protein, Pfam accession or interval."
                )
            seen.add(identifier)
            domain_rows.append(
                {
                    "protein_id": identifier,
                    "domain_authority": "Pfam",
                    "assessment_status": "ASSESSED_WITH_HIT",
                    "domain_id": domain,
                    "domain_name": fields[0],
                    "start": start,
                    "end": end,
                    "score": fields[13],
                    "e_value": fields[12],
                    "evidence_source": "HMMER hmmscan --cut_ga",
                    "evidence_reference": config.pfam_hmm,
                }
            )
            feature_rows.append(
                {
                    "protein_id": identifier,
                    "feature_type": "PFAM_DOMAIN",
                    "feature_id": domain,
                    "evidence_status": "ASSESSED_WITH_FEATURE",
                }
            )
    for identifier in sizes.keys() - seen:
        domain_rows.append(
            {
                "protein_id": identifier,
                "domain_authority": "Pfam",
                "assessment_status": "ASSESSED_NO_HIT",
                "domain_id": "",
                "domain_name": "",
                "start": "",
                "end": "",
                "score": "",
                "e_value": "",
                "evidence_source": "HMMER hmmscan --cut_ga",
                "evidence_reference": config.pfam_hmm,
            }
        )
    write_tsv_atomic(
        path=output / "domains.tsv",
        fieldnames=(
            "protein_id",
            "domain_authority",
            "assessment_status",
            "domain_id",
            "domain_name",
            "start",
            "end",
            "score",
            "e_value",
            "evidence_source",
            "evidence_reference",
        ),
        records=domain_rows,
    )
    unique = pd.DataFrame(
        data=feature_rows, columns=["protein_id", "feature_type", "feature_id", "evidence_status"]
    ).drop_duplicates()
    unique.to_csv(path_or_buf=output / "query_features.tsv", sep="\t", index=False)


def score_family(
    *, config: WorkbenchConfig, source: Path, features: Path | None, parameters: dict, output: Path
) -> None:
    """Apply frozen reference signatures; unassessed structures and domains stay unknown.

    Args:
        config: Completed reference result and bounded scoring inputs.
        source: New protein FASTA.
        features: Optional positive Pfam feature assessments from a completed annotation job.
        parameters: Selected comparisons and optional within-comparison tier.
        output: Owned job directory.
    """
    from protein_signature_app.backend import (
        query_dataframe,
        resolve_database_for_app,
        result_inventory_identity,
    )
    from protein_signatures.query import query_signature_evidence

    if not config.resource:
        raise InputValidationError(
            "Set resource to the completed reference result in the operator profile."
        )
    resource = Path(config.resource)
    database = resolve_database_for_app(
        resource=resource, inventory_identity=result_inventory_identity(resource=resource)
    )
    with open_catalogue(config=config, with_index=False) as connection:
        index_fasta(
            connection=connection, path=source, name="family", maximum=config.max_seed_sequences
        )
        write_fasta(
            connection=connection,
            query="SELECT sequence_id,sequence FROM family ORDER BY sequence_id",
            path=output / "sequences.fasta",
        )
    for table in ("signatures", "comparisons"):
        frame = query_dataframe(database=database, sql=f"SELECT * FROM {table}")
        frame.to_csv(path_or_buf=output / f"reference_{table}.tsv", sep="\t", index=False)
    selected = parameters.get("comparison_ids", [])
    include = parameters.get("include_within_comparison", False)
    if (
        not isinstance(selected, list)
        or any(not isinstance(v, str) for v in selected)
        or not isinstance(include, bool)
    ):
        raise InputValidationError(
            "Supply comparison_ids as a list and include_within_comparison as a Boolean."
        )
    query_signature_evidence(
        signatures_path=output / "reference_signatures.tsv",
        comparisons_path=output / "reference_comparisons.tsv",
        query_fasta=output / "sequences.fasta",
        query_features_path=features,
        output_dir=output / "signature_evidence",
        comparison_ids=tuple(selected),
        include_within_comparison=include,
    )


def prepare_campaign(*, config: WorkbenchConfig, source: Path, seeds: Path, output: Path) -> None:
    """Prepare a deduplicated expansion and curation starter without assigning functions.

    Args:
        config: Workspace and complete input bounds.
        source: Selected DeepClust candidates.
        seeds: Original seed proteins.
        output: New job output directory.
    """
    with open_catalogue(config=config, with_index=False) as connection:
        index_fasta(connection=connection, path=source, name="family", maximum=config.max_sequences)
        index_fasta(
            connection=connection, path=seeds, name="seeds", maximum=config.max_seed_sequences
        )
        connection.execute(
            (
                "CREATE TABLE combined AS SELECT 'seed__'||sequence_id AS "
                "protein_id,sequence_id AS source_id,sequence,sequence_sha256,'SEED' AS"
                " source_role,0 AS priority FROM seeds UNION ALL SELECT "
                "'deepclust__'||sequence_id,sequence_id,sequence,sequence_sha256,'CANDIDATE',1"
                " FROM family"
            )
        )
        connection.execute(
            (
                "CREATE TABLE selected AS SELECT * FROM combined QUALIFY row_number() "
                "OVER(PARTITION BY sequence_sha256 ORDER BY priority,source_id)=1"
            )
        )
        write_fasta(
            connection=connection,
            query="SELECT protein_id,sequence FROM selected ORDER BY protein_id",
            path=output / "expanded.fasta",
        )
        connection.sql(
            query=(
                "SELECT "
                "c.protein_id,c.source_id,c.source_role,c.sequence_sha256,s.protein_id "
                "AS retained_protein_id FROM combined c JOIN selected s "
                "USING(sequence_sha256) ORDER BY c.protein_id"
            )
        ).write_csv(file_name=str(output / "sequence_provenance.tsv"), sep="\t", header=True)
        connection.sql(
            query=(
                "SELECT protein_id,'' AS label_id,'PROPOSED' AS "
                "curation_status,'UNREVIEWED' AS evidence_status,'DeepClust expansion' "
                "AS evidence_source,'' AS evidence_reference,'' AS "
                "component_role,'Review function, domains, structure and independent "
                "homology units before analysis' AS curation_reason FROM selected ORDER"
                " BY protein_id"
            )
        ).write_csv(
            file_name=str(output / "label_assignments_to_review.tsv"), sep="\t", header=True
        )
    (output / "README.txt").write_text(
        data=(
            "EXPANDED PROTEIN ANALYSIS — REVIEW BEFORE A NEW CAMPAIGN\n\n"
            "expanded.fasta contains exact-sequence-deduplicated seeds and "
            "candidates. Seed sequences take precedence.\nsequence_provenance.tsv "
            "preserves every source ID and duplicate-to-retained mapping.\nFill "
            "label_assignments_to_review.tsv using your chosen profile, including "
            "reviewed background controls.\nSimilarity and DeepClust membership do "
            "not establish function. No new protein is automatically a target or "
            "control.\nRecompute homology/independence blocks across seeds AND "
            "expanded proteins (for example using OrthoFinder).\nNever split members"
            " of the same related family across discovery and validation, or call "
            "this external validation.\nAdditional relatives may increase protein "
            "counts without increasing independent sample size.\nUse a new campaign "
            "ID and output directory. The reference campaign remains unchanged.\n\n"
            "After review, use protein-signatures initialise --help to create a "
            "campaign with explicit labels, domains and homology resources.\nUse "
            "protein-signatures query-signatures --help to inspect frozen reference"
            " evidence without retraining.\n"
        ),
        encoding="utf-8",
    )
