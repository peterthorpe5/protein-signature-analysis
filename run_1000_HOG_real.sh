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
  --slurm-memory 128G \
  --slurm-time 04:00:00 \
  --threads 8 \
  --slurm-dry-run
