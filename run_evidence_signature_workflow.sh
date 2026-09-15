#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
    cat <<'EOF'
Usage:
  ./run_evidence_signature_workflow.sh \
    --work-dir /persistent/path/campaign \
    --campaign-id protein_type_campaign \
    --sequences-fasta /absolute/path/proteins.faa \
    --profile e3 \
    --evidence-rules e3 \
    [--protein-metadata /absolute/path/protein_metadata.tsv] \
    [--domains /absolute/path/domains.tsv] \
    [--structures /absolute/path/structures.tsv] \
    [--external-annotations /absolute/path/external_annotations.tsv] \
    [--seed-assignments /absolute/path/reviewed_seed_labels.tsv] \
    [--orthofinder-resource /absolute/path/resource | \
     --orthofinder-results /absolute/path/Results] \
    [--accept-provisional-evidence-labels] \
    [--submit-slurm --slurm-account barton --slurm-partition barton \
     --slurm-memory 128G --slurm-time 2-00:00:00] \
    [--threads 24]

Without --accept-provisional-evidence-labels the Snakemake DAG stops after the
checksummed evidence bundle. Rerun the same command with the explicit acceptance
flag to initialise, analyse and verify the complete provisional campaign.

The profile and rules may describe any protein type. The built-in E3 pair is the
default. OrthoFinder 2.5.5/3 output is consumed but never generated. Structural
evidence remains enabled by default; provide verified structures or AlphaFold
accessions as required by the selected profile.
EOF
}

require_value() {
    if (( $# < 2 )) || [[ -z "${2:-}" || "${2:-}" == --* ]]; then
        echo "Option $1 requires a non-empty value." >&2
        exit 2
    fi
}

validate_slurm_time() {
    local value="$1"
    if [[ "${value}" =~ ^([0-9]+)-([0-9]{2}):([0-9]{2}):([0-9]{2})$ ]]; then
        (( 10#${BASH_REMATCH[2]} <= 23 )) || return 1
        (( 10#${BASH_REMATCH[3]} <= 59 )) || return 1
        (( 10#${BASH_REMATCH[4]} <= 59 )) || return 1
        return 0
    fi
    if [[ "${value}" =~ ^([0-9]+):([0-9]{2}):([0-9]{2})$ ]]; then
        (( 10#${BASH_REMATCH[2]} <= 59 )) || return 1
        (( 10#${BASH_REMATCH[3]} <= 59 )) || return 1
        return 0
    fi
    return 1
}

WORK_DIR=""
CAMPAIGN_ID=""
SEQUENCES_FASTA=""
PROFILE="e3"
EVIDENCE_RULES="e3"
PROTEIN_METADATA=""
REVIEW_CONTEXT=""
DOMAINS=""
STRUCTURES=""
TEMPLATE_LABELS=""
SEED_ASSIGNMENTS=""
SEED_CATALOGUE=""
EXTERNAL_ANNOTATIONS=""
REDUNDANCY_CLUSTERS=""
FEATURES=""
STRUCTURE_COMPARISONS=""
STRUCTURAL_ALIGNMENT_RESOURCE=""
ALPHAFOLD_ACCESSIONS=""
ORTHOFINDER_RESOURCE=""
ORTHOFINDER_RESULTS=""
ORTHOFINDER_GROUP_TYPE="HOG"
ORTHOFINDER_HIERARCHY_NODE="N0"
ORTHOFINDER_HIERARCHY_NODE_SEEN="false"
ORTHOFINDER_RUN_ID="evidence_labelling"
ENABLE_ALPHAFOLD="false"
ENABLE_FOLDSEEK="true"
FOLDSEEK_MAXIMUM_HITS="1000"
ACCEPT_PROVISIONAL="false"
RESUME="false"
DRY_RUN="false"
CONDA_ENVIRONMENT="protein_signature_analysis"
THREADS="1"
MEMORY_MB="128000"
RUNTIME_MINUTES="2880"
LOG_LEVEL="INFO"
SUBMIT_SLURM="false"
SLURM_ACCOUNT="barton"
SLURM_PARTITION="barton"
SLURM_MEMORY="128G"
SLURM_TIME="2-00:00:00"
SLURM_LOG_DIR=""
SKIP_ENVIRONMENT_SYNC="false"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --work-dir) require_value "$@"; WORK_DIR="${2:-}"; shift 2 ;;
        --campaign-id) require_value "$@"; CAMPAIGN_ID="${2:-}"; shift 2 ;;
        --sequences-fasta) require_value "$@"; SEQUENCES_FASTA="${2:-}"; shift 2 ;;
        --profile) require_value "$@"; PROFILE="${2:-}"; shift 2 ;;
        --evidence-rules) require_value "$@"; EVIDENCE_RULES="${2:-}"; shift 2 ;;
        --protein-metadata) require_value "$@"; PROTEIN_METADATA="${2:-}"; shift 2 ;;
        --review-context) require_value "$@"; REVIEW_CONTEXT="${2:-}"; shift 2 ;;
        --domains) require_value "$@"; DOMAINS="${2:-}"; shift 2 ;;
        --structures) require_value "$@"; STRUCTURES="${2:-}"; shift 2 ;;
        --template-labels) require_value "$@"; TEMPLATE_LABELS="${2:-}"; shift 2 ;;
        --seed-assignments) require_value "$@"; SEED_ASSIGNMENTS="${2:-}"; shift 2 ;;
        --seed-catalogue) require_value "$@"; SEED_CATALOGUE="${2:-}"; shift 2 ;;
        --external-annotations) require_value "$@"; EXTERNAL_ANNOTATIONS="${2:-}"; shift 2 ;;
        --redundancy-clusters) require_value "$@"; REDUNDANCY_CLUSTERS="${2:-}"; shift 2 ;;
        --features) require_value "$@"; FEATURES="${2:-}"; shift 2 ;;
        --structure-comparisons) require_value "$@"; STRUCTURE_COMPARISONS="${2:-}"; shift 2 ;;
        --structural-alignment-resource)
            require_value "$@"; STRUCTURAL_ALIGNMENT_RESOURCE="${2:-}"; shift 2 ;;
        --alphafold-accessions) require_value "$@"; ALPHAFOLD_ACCESSIONS="${2:-}"; shift 2 ;;
        --orthofinder-resource) require_value "$@"; ORTHOFINDER_RESOURCE="${2:-}"; shift 2 ;;
        --orthofinder-results) require_value "$@"; ORTHOFINDER_RESULTS="${2:-}"; shift 2 ;;
        --orthofinder-group-type)
            require_value "$@"; ORTHOFINDER_GROUP_TYPE="${2:-}"; shift 2 ;;
        --orthofinder-hierarchy-node)
            require_value "$@"
            ORTHOFINDER_HIERARCHY_NODE="${2:-}"
            ORTHOFINDER_HIERARCHY_NODE_SEEN="true"
            shift 2
            ;;
        --orthofinder-run-id) require_value "$@"; ORTHOFINDER_RUN_ID="${2:-}"; shift 2 ;;
        --enable-alphafold) ENABLE_ALPHAFOLD="true"; shift ;;
        --disable-foldseek) ENABLE_FOLDSEEK="false"; shift ;;
        --foldseek-maximum-hits)
            require_value "$@"; FOLDSEEK_MAXIMUM_HITS="${2:-}"; shift 2 ;;
        --accept-provisional-evidence-labels) ACCEPT_PROVISIONAL="true"; shift ;;
        --resume) RESUME="true"; shift ;;
        --dry-run) DRY_RUN="true"; shift ;;
        --conda-environment) require_value "$@"; CONDA_ENVIRONMENT="${2:-}"; shift 2 ;;
        --threads) require_value "$@"; THREADS="${2:-}"; shift 2 ;;
        --memory-mb) require_value "$@"; MEMORY_MB="${2:-}"; shift 2 ;;
        --runtime-minutes) require_value "$@"; RUNTIME_MINUTES="${2:-}"; shift 2 ;;
        --log-level) require_value "$@"; LOG_LEVEL="${2:-}"; shift 2 ;;
        --submit-slurm) SUBMIT_SLURM="true"; shift ;;
        --slurm-account) require_value "$@"; SLURM_ACCOUNT="${2:-}"; shift 2 ;;
        --slurm-partition) require_value "$@"; SLURM_PARTITION="${2:-}"; shift 2 ;;
        --slurm-memory) require_value "$@"; SLURM_MEMORY="${2:-}"; shift 2 ;;
        --slurm-time) require_value "$@"; SLURM_TIME="${2:-}"; shift 2 ;;
        --slurm-log-dir) require_value "$@"; SLURM_LOG_DIR="${2:-}"; shift 2 ;;
        --skip-environment-sync) SKIP_ENVIRONMENT_SYNC="true"; shift ;;
        --help|-h) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
done

if [[ -z "${WORK_DIR}" || -z "${CAMPAIGN_ID}" || -z "${SEQUENCES_FASTA}" ]]; then
    echo "--work-dir, --campaign-id and --sequences-fasta are required." >&2
    usage >&2
    exit 2
fi
if [[ "${WORK_DIR}" != /* || "${WORK_DIR}" == "/" || \
        "${WORK_DIR}" == "${HOME:-}" ]]; then
    echo "--work-dir must be a safe absolute persistent directory." >&2
    exit 2
fi
if [[ ! -s "${SEQUENCES_FASTA}" || "${SEQUENCES_FASTA}" != /* ]]; then
    echo "--sequences-fasta must be an existing non-empty absolute file." >&2
    exit 2
fi
for number in "${THREADS}" "${MEMORY_MB}" "${RUNTIME_MINUTES}" \
        "${FOLDSEEK_MAXIMUM_HITS}"; do
    if [[ ! "${number}" =~ ^[1-9][0-9]*$ ]]; then
        echo "Thread, memory, runtime and maximum-hit values must be positive integers." >&2
        exit 2
    fi
done
if [[ ! "${CAMPAIGN_ID}" =~ ^[A-Za-z0-9][A-Za-z0-9_.:+/-]{0,254}$ ]]; then
    echo "--campaign-id is not a valid identifier." >&2
    exit 2
fi
if [[ ! "${ORTHOFINDER_GROUP_TYPE}" =~ ^(HOG|LEGACY_ORTHOGROUP)$ ]]; then
    echo "--orthofinder-group-type must be HOG or LEGACY_ORTHOGROUP." >&2
    exit 2
fi
if [[ "${ORTHOFINDER_GROUP_TYPE}" == "LEGACY_ORTHOGROUP" ]]; then
    if [[ "${ORTHOFINDER_HIERARCHY_NODE_SEEN}" == "true" ]]; then
        echo "Do not supply --orthofinder-hierarchy-node for legacy orthogroups." >&2
        exit 2
    fi
    ORTHOFINDER_HIERARCHY_NODE=""
elif [[ -z "${ORTHOFINDER_HIERARCHY_NODE}" ]]; then
    echo "HOG membership requires --orthofinder-hierarchy-node." >&2
    exit 2
fi
if [[ -n "${ORTHOFINDER_RESOURCE}" && -n "${ORTHOFINDER_RESULTS}" ]]; then
    echo "Choose either --orthofinder-resource or --orthofinder-results." >&2
    exit 2
fi
if [[ "${ENABLE_ALPHAFOLD}" == "true" && -z "${ALPHAFOLD_ACCESSIONS}" ]]; then
    echo "--enable-alphafold requires --alphafold-accessions." >&2
    exit 2
fi
if [[ ! "${CONDA_ENVIRONMENT}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ || \
        ! "${LOG_LEVEL}" =~ ^(DEBUG|INFO|WARNING|ERROR|CRITICAL)$ ]]; then
    echo "Invalid Conda environment name or log level." >&2
    exit 2
fi

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly WORKFLOW_STATE="${WORK_DIR}/workflow_state/evidence"
readonly SLURM_WORKER="${SCRIPT_DIR}/slurm/run_completed_e3_workflow.sbatch"
SLURM_LOG_DIR="${SLURM_LOG_DIR:-${WORK_DIR}/slurm_logs}"

INPUT_FILES=(
    "${PROTEIN_METADATA}" "${REVIEW_CONTEXT}" "${DOMAINS}" "${STRUCTURES}"
    "${TEMPLATE_LABELS}" "${SEED_ASSIGNMENTS}" "${SEED_CATALOGUE}"
    "${EXTERNAL_ANNOTATIONS}" "${REDUNDANCY_CLUSTERS}" "${FEATURES}"
    "${STRUCTURE_COMPARISONS}" "${ALPHAFOLD_ACCESSIONS}"
)
for input_file in "${INPUT_FILES[@]}"; do
    if [[ -n "${input_file}" && ("${input_file}" != /* || ! -s "${input_file}") ]]; then
        echo "Optional input must be an existing non-empty absolute file: ${input_file}" >&2
        exit 2
    fi
done
for input_dir in "${STRUCTURAL_ALIGNMENT_RESOURCE}" "${ORTHOFINDER_RESOURCE}" \
        "${ORTHOFINDER_RESULTS}"; do
    if [[ -n "${input_dir}" && ("${input_dir}" != /* || ! -d "${input_dir}") ]]; then
        echo "Optional resource must be an existing absolute directory: ${input_dir}" >&2
        exit 2
    fi
done
for specification in "${PROFILE}" "${EVIDENCE_RULES}"; do
    if [[ "${specification}" == */* && \
            ("${specification}" != /* || ! -s "${specification}") ]]; then
        echo "Custom profile/rules must be existing non-empty absolute files." >&2
        exit 2
    fi
done
if [[ "${PROFILE}" != */* && \
        ! "${PROFILE}" =~ ^[A-Za-z0-9][A-Za-z0-9_.:+-]{0,254}$ ]]; then
    echo "Built-in --profile contains unsafe characters." >&2
    exit 2
fi
if [[ "${EVIDENCE_RULES}" != */* && \
        ! "${EVIDENCE_RULES}" =~ ^[A-Za-z0-9][A-Za-z0-9_.:+-]{0,254}$ ]]; then
    echo "Built-in --evidence-rules contains unsafe characters." >&2
    exit 2
fi

WORKER_ARGUMENTS=(
    --work-dir "${WORK_DIR}"
    --campaign-id "${CAMPAIGN_ID}"
    --sequences-fasta "${SEQUENCES_FASTA}"
    --profile "${PROFILE}"
    --evidence-rules "${EVIDENCE_RULES}"
    --orthofinder-group-type "${ORTHOFINDER_GROUP_TYPE}"
    --orthofinder-hierarchy-node "${ORTHOFINDER_HIERARCHY_NODE}"
    --orthofinder-run-id "${ORTHOFINDER_RUN_ID}"
    --foldseek-maximum-hits "${FOLDSEEK_MAXIMUM_HITS}"
    --conda-environment "${CONDA_ENVIRONMENT}"
    --threads "${THREADS}"
    --memory-mb "${MEMORY_MB}"
    --runtime-minutes "${RUNTIME_MINUTES}"
    --log-level "${LOG_LEVEL}"
)
OPTION_NAMES=(
    --protein-metadata --review-context --domains --structures --template-labels
    --seed-assignments --seed-catalogue --external-annotations --redundancy-clusters
    --features --structure-comparisons --alphafold-accessions
    --structural-alignment-resource --orthofinder-resource --orthofinder-results
)
OPTION_VALUES=(
    "${PROTEIN_METADATA}" "${REVIEW_CONTEXT}" "${DOMAINS}" "${STRUCTURES}"
    "${TEMPLATE_LABELS}" "${SEED_ASSIGNMENTS}" "${SEED_CATALOGUE}"
    "${EXTERNAL_ANNOTATIONS}" "${REDUNDANCY_CLUSTERS}" "${FEATURES}"
    "${STRUCTURE_COMPARISONS}" "${ALPHAFOLD_ACCESSIONS}"
    "${STRUCTURAL_ALIGNMENT_RESOURCE}" "${ORTHOFINDER_RESOURCE}" "${ORTHOFINDER_RESULTS}"
)
for index in "${!OPTION_NAMES[@]}"; do
    if [[ -n "${OPTION_VALUES[$index]}" ]]; then
        WORKER_ARGUMENTS+=("${OPTION_NAMES[$index]}" "${OPTION_VALUES[$index]}")
    fi
done
[[ "${ENABLE_ALPHAFOLD}" == "true" ]] && WORKER_ARGUMENTS+=(--enable-alphafold)
[[ "${ENABLE_FOLDSEEK}" == "false" ]] && WORKER_ARGUMENTS+=(--disable-foldseek)
[[ "${ACCEPT_PROVISIONAL}" == "true" ]] && \
    WORKER_ARGUMENTS+=(--accept-provisional-evidence-labels)
[[ "${RESUME}" == "true" ]] && WORKER_ARGUMENTS+=(--resume)
[[ "${DRY_RUN}" == "true" ]] && WORKER_ARGUMENTS+=(--dry-run)
[[ "${SKIP_ENVIRONMENT_SYNC}" == "true" ]] && \
    WORKER_ARGUMENTS+=(--skip-environment-sync)

if [[ "${SUBMIT_SLURM}" == "true" ]]; then
    command -v sbatch >/dev/null 2>&1 || {
        echo "sbatch is unavailable; submit from a Slurm login node." >&2
        exit 2
    }
    if [[ ! "${SLURM_MEMORY}" =~ ^[1-9][0-9]*[KMGT]?$ ]] || \
            ! validate_slurm_time "${SLURM_TIME}"; then
        echo "Invalid Slurm memory or time request." >&2
        exit 2
    fi
    if [[ ! "${SLURM_ACCOUNT}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ || \
            ! "${SLURM_PARTITION}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]]; then
        echo "Invalid Slurm account or partition." >&2
        exit 2
    fi
    mkdir -p -- "${SLURM_LOG_DIR}"
    submission="$(
        env -u SLURM_CPUS_PER_TASK sbatch \
            --parsable \
            --job-name=protein_signature_evidence \
            "--account=${SLURM_ACCOUNT}" \
            "--partition=${SLURM_PARTITION}" \
            "--mem=${SLURM_MEMORY}" \
            "--time=${SLURM_TIME}" \
            "--cpus-per-task=${THREADS}" \
            "--export=ALL,PROTEIN_SIGNATURE_REQUESTED_CPUS=${THREADS}" \
            "--chdir=${SCRIPT_DIR}" \
            "--output=${SLURM_LOG_DIR}/protein_signature_evidence_%j.out" \
            "--error=${SLURM_LOG_DIR}/protein_signature_evidence_%j.err" \
            "${SLURM_WORKER}" "${SCRIPT_DIR}/run_evidence_signature_workflow.sh" \
            "${WORKER_ARGUMENTS[@]}"
    )"
    job_id="${submission%%;*}"
    [[ "${job_id}" =~ ^[0-9]+$ ]] || {
        echo "sbatch returned an unexpected job identifier: ${submission}" >&2
        exit 2
    }
    echo "Submitted evidence-led workflow job ${job_id}."
    echo "Stdout: ${SLURM_LOG_DIR}/protein_signature_evidence_${job_id}.out"
    echo "Stderr: ${SLURM_LOG_DIR}/protein_signature_evidence_${job_id}.err"
    echo "Monitor: squeue --job ${job_id}"
    exit 0
fi

if ! command -v conda >/dev/null 2>&1; then
    echo "conda is required but was not found on PATH." >&2
    exit 2
fi
mkdir -p -- "${WORK_DIR}" "${WORKFLOW_STATE}"
if [[ "${SKIP_ENVIRONMENT_SYNC}" != "true" ]]; then
    if conda run --name "${CONDA_ENVIRONMENT}" python --version >/dev/null 2>&1; then
        conda env update --file "${SCRIPT_DIR}/environment.yml" \
            --name "${CONDA_ENVIRONMENT}" --prune
    else
        conda env create --file "${SCRIPT_DIR}/environment.yml" \
            --name "${CONDA_ENVIRONMENT}"
    fi
    conda run --name "${CONDA_ENVIRONMENT}" python -m pip install \
        --no-deps --force-reinstall --editable "${SCRIPT_DIR}"
fi
conda run --name "${CONDA_ENVIRONMENT}" snakemake --version >/dev/null

SNAKEMAKE_COMMAND=(
    conda run --no-capture-output --name "${CONDA_ENVIRONMENT}"
    snakemake
    --snakefile "${SCRIPT_DIR}/workflow/EvidenceSnakefile"
    --directory "${WORKFLOW_STATE}"
    --profile "${SCRIPT_DIR}/profiles/local"
    --cores "${THREADS}"
    --jobs 1
    --rerun-incomplete
    --printshellcmds
    all
    --config
    "evidence_work_dir=${WORK_DIR}"
    "evidence_campaign_id=${CAMPAIGN_ID}"
    "evidence_sequences_fasta=${SEQUENCES_FASTA}"
    "evidence_profile=${PROFILE}"
    "evidence_rules=${EVIDENCE_RULES}"
    "evidence_protein_metadata=${PROTEIN_METADATA}"
    "evidence_review_context=${REVIEW_CONTEXT}"
    "evidence_domains=${DOMAINS}"
    "evidence_structures=${STRUCTURES}"
    "evidence_template_labels=${TEMPLATE_LABELS}"
    "evidence_seed_assignments=${SEED_ASSIGNMENTS}"
    "evidence_seed_catalogue=${SEED_CATALOGUE}"
    "evidence_external_annotations=${EXTERNAL_ANNOTATIONS}"
    "evidence_redundancy_clusters=${REDUNDANCY_CLUSTERS}"
    "evidence_features=${FEATURES}"
    "evidence_structure_comparisons=${STRUCTURE_COMPARISONS}"
    "evidence_structural_alignment_resource=${STRUCTURAL_ALIGNMENT_RESOURCE}"
    "evidence_alphafold_accessions=${ALPHAFOLD_ACCESSIONS}"
    "evidence_orthofinder_resource=${ORTHOFINDER_RESOURCE}"
    "evidence_orthofinder_results=${ORTHOFINDER_RESULTS}"
    "evidence_orthofinder_group_type=${ORTHOFINDER_GROUP_TYPE}"
    "evidence_orthofinder_hierarchy_node=${ORTHOFINDER_HIERARCHY_NODE}"
    "evidence_orthofinder_run_id=${ORTHOFINDER_RUN_ID}"
    "evidence_enable_alphafold=${ENABLE_ALPHAFOLD}"
    "evidence_enable_foldseek=${ENABLE_FOLDSEEK}"
    "evidence_foldseek_maximum_hits=${FOLDSEEK_MAXIMUM_HITS}"
    "evidence_accept_provisional=${ACCEPT_PROVISIONAL}"
    "evidence_resume=${RESUME}"
    "evidence_threads=${THREADS}"
    "evidence_memory_mb=${MEMORY_MB}"
    "evidence_runtime_minutes=${RUNTIME_MINUTES}"
    "evidence_log_level=${LOG_LEVEL}"
)
[[ "${DRY_RUN}" == "true" ]] && SNAKEMAKE_COMMAND+=(--dry-run)
"${SNAKEMAKE_COMMAND[@]}"

if [[ "${ACCEPT_PROVISIONAL}" == "true" ]]; then
    echo "Completed and verified provisional evidence-led result: ${WORK_DIR}/result"
else
    echo "Evidence bundle ready; analysis was not started: ${WORK_DIR}/evidence_label_bundle"
    echo "Inspect the TSV/XLSX audits, then rerun with --accept-provisional-evidence-labels."
fi
