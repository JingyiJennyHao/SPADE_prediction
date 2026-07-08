# GMM HPC Run

The original `gmm_sim.Rmd` is unchanged. The HPC run uses these new files:

- `gmm_core.R`: shared simulation and GMM functions copied from the notebook.
- `gmm_start_worker.R`: runs one start point and appends one row to a loop CSV.
- `gmm_start_worker_lsf.sh`: LSF array wrapper that maps `LSB_JOBINDEX` to `start_id`.
- `gmm_hpc_driver.R`: submits each loop, waits until 100 converged starts are recorded, picks the best beta, then submits the next loop.

## Run on an LSF cluster

From this directory:

```bash
Rscript gmm_hpc_driver.R \
  --seed 1 \
  --loops 5 \
  --starts 120 \
  --target 100 \
  --out-dir hpc_runs/seed1
```

Optional arguments:

```bash
Rscript gmm_hpc_driver.R \
  --seed 1 \
  --loops 5 \
  --starts 120 \
  --target 100 \
  --queue normal \
  --bsub-extra "-M 8000 -R 'rusage[mem=8000]'" \
  --poll-seconds 60 \
  --max-wait-hours 72 \
  --out-dir hpc_runs/seed1
```

## Outputs

For loop `k`, all 120 start-point jobs append to one file:

```text
hpc_runs/seed1/loop_0k_start_results_seed1.csv
```

The driver advances once at least `--target` converged rows exist in that loop file. It writes the best beta from each loop to:

```text
hpc_runs/seed1/base_beta_loop{k+1}.txt
```

The per-loop best beta and objective are also summarized in:

```text
hpc_runs/seed1/loop_summary_seed1.csv
```

No per-start `.rds` files are written.
