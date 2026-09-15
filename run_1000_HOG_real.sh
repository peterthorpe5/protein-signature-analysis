RUN_ROOT="/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_E3_protac/analysis/e3_end_to_end_runs/grant_aligned_corrected_expression_structural_all1972_v0_16_0_20260909"

SIGNATURE_WORK="/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_E3_protac/analysis/protein_signature_runs/e3_all1972_v0_1_0_20260914"

./run_completed_e3_workflow.sh \
  --phase prepare \
  --run-root "${RUN_ROOT}" \
  --work-dir "${SIGNATURE_WORK}" \
  --minimum-mean-plddt 50 \
  --submit-slurm \
  --slurm-account barton \
  --slurm-partition general \
  --slurm-memory 28G \
  --slurm-time 04:00:00 \
  --threads 8 \
  --slurm-dry-run




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




  RUN_ROOT="/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_E3_protac/analysis/e3_end_to_end_runs/grant_aligned_corrected_expression_structural_all1972_v0_16_0_20260909"

SIGNATURE_TEST_WORK="/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_E3_protac/analysis/protein_signature_runs/e3_all1972_automated_smoke_test_v0_1_0_20260915"

./run_completed_e3_workflow.sh \
  --phase all \
  --run-root "${RUN_ROOT}" \
  --work-dir "${SIGNATURE_TEST_WORK}" \
  --campaign-id "e3_all1972_automated_smoke_test_20260915" \
  --automated-test-labels \
  --test-target-label ALL \
  --test-samples-per-class 20 \
  --minimum-mean-plddt 50 \
  --submit-slurm \
  --slurm-account barton \
  --slurm-partition barton \
  --slurm-memory 128G \
  --slurm-time 2-00:00:00 \
  --threads 24

#######################################
RUN_ROOT="/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_E3_protac/analysis/e3_end_to_end_runs/grant_aligned_corrected_expression_structural_all1972_v0_16_0_20260909"

SIGNATURE_WORK="/gpfs/uod-scale-01/cluster/gjb_lab/pthorpe001/2026_E3_protac/analysis/protein_signature_runs/e3_all1972_v0_1_0_20260914"

REVIEWED_LABELS="${SIGNATURE_WORK}/reviewed_label_assignments.tsv"

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

