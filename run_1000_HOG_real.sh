
# RUN ALL AT ONCE
cd /gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_E3_protac/protein-signature-analysis

git pull --ff-only origin main

RUN_ROOT="/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_E3_protac/analysis/e3_end_to_end_runs/grant_aligned_corrected_expression_structural_all1972_v0_16_0_20260909"

SIGNATURE_PRODUCTION_WORK="/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_E3_protac/analysis/protein_signature_runs/e3_all1972_evidence_provisional_v0_1_0_20260915"

#  UNLOCK FIST

cd /gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_E3_protac/protein-signature-analysis

RUN_ROOT="/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_E3_protac/analysis/e3_end_to_end_runs/grant_aligned_corrected_expression_structural_all1972_v0_16_0_20260909"

SIGNATURE_PRODUCTION_WORK="/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_E3_protac/analysis/protein_signature_runs/e3_all1972_evidence_provisional_v0_1_0_20260915"

E3_STATE_DIR="${SIGNATURE_PRODUCTION_WORK}/workflow_state/e3"

EVIDENCE_LABELS="${SIGNATURE_PRODUCTION_WORK}/evidence_label_bundle/label_assignments.tsv"

EVIDENCE_APPROVAL="${E3_STATE_DIR}/02_label_review/EVIDENCE_PROVISIONAL.REVIEW_APPROVED.json"

ls -lh -- "${EVIDENCE_LABELS}" "${EVIDENCE_APPROVAL}"

conda run --no-capture-output --name protein_signature_analysis \
  snakemake \
  --snakefile "${PWD}/workflow/E3Snakefile" \
  --directory "${E3_STATE_DIR}" \
  --unlock \
  --config \
  "e3_run_root=${RUN_ROOT}" \
  "e3_work_dir=${SIGNATURE_PRODUCTION_WORK}" \
  "e3_campaign_id=e3_all1972_evidence_provisional_20260915" \
  "e3_profile=e3" \
  "e3_minimum_mean_plddt=50" \
  "e3_reviewed_labels=${EVIDENCE_LABELS}" \
  "e3_review_approval=${EVIDENCE_APPROVAL}" \
  "e3_automated_test_labels=false" \
  "e3_evidence_led_labels=true" \
  "e3_accept_provisional_evidence_labels=true" \
  "e3_evidence_rules=e3" \
  "e3_threads=24"


./run_completed_e3_workflow.sh \
  --phase all \
  --run-root "${RUN_ROOT}" \
  --work-dir "${SIGNATURE_PRODUCTION_WORK}" \
  --campaign-id "e3_all1972_evidence_provisional_20260915" \
  --evidence-led-labels \
  --accept-provisional-evidence-labels \
  --evidence-rules e3 \
  --minimum-mean-plddt 50 \
  --submit-slurm \
  --slurm-account barton \
  --slurm-partition barton \
  --slurm-memory 128G \
  --slurm-time 2-00:00:00 \
  --threads 24



# RUN IN STAGES


RUN_ROOT="/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_E3_protac/analysis/e3_end_to_end_runs/grant_aligned_corrected_expression_structural_all1972_v0_16_0_20260909"
SIGNATURE_WORK="/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_E3_protac/analysis/protein_signature_runs/e3_all1972_v0_1_0_20260914"

./run_completed_e3_workflow.sh \
  --phase prepare \
  --run-root "${RUN_ROOT}" \
  --work-dir "${SIGNATURE_WORK}" \
  --campaign-id "e3_all1972_signatures_20260914" \
  --minimum-mean-plddt 50 \
  --submit-slurm \
  --slurm-account barton \
  --slurm-partition barton \
  --slurm-memory 64G \
  --slurm-time 04:00:00 \
  --threads 18



#This verifies and reuses prepared_inputs, and adopts the already copied
#reviewed_label_assignments.tsv without overwriting it. Curate that TSV before
#approval. Every populated target class requires its profile-resolved matched
#control label.

./run_completed_e3_workflow.sh \
  --phase approve \
  --work-dir "${SIGNATURE_WORK}" \
  --curator "Peter Thorpe" \
  --review-note "Reviewed E3 subclasses, component roles and matched controls"

#After approval, submit the remaining DAG:

./run_completed_e3_workflow.sh \
  --phase all \
  --run-root "${RUN_ROOT}" \
  --work-dir "${SIGNATURE_WORK}" \
  --campaign-id "e3_all1972_signatures_20260914" \
  --submit-slurm \
  --slurm-account barton \
  --slurm-partition barton \
  --slurm-memory 128G \
  --slurm-time 2-00:00:00 \
  --threads 24




