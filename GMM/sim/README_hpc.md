# GMM HPC Run

The original `gmm_sim.Rmd` is unchanged. The HPC run uses these new files:

- `gmm_core.R`: shared simulation and GMM functions copied from the notebook.
- `gmm_start_worker.R`: runs one start point and appends one row to a loop CSV.
- `gmm_start_worker_lsf.sh`: LSF array wrapper that maps `LSB_JOBINDEX` to `start_id`.
- `gmm_hpc_driver.R`: submits each loop, waits until 100 converged starts are recorded, picks the best beta, then submits the next loop.

## Run on an LSF cluster

From this directory:

```bash
./run_loop.sh
```

Default `run_loop.sh` settings are:

```text
SEED=1
LOOPS=20
STARTS=120
TARGET=100
TOL=1e-7
QUEUE=serial
WALL_TIME=14:00
MEMORY_GB=8
R_MODULE=R/4.4.0
R_LIBS_USER_PATH=$HOME/R/x86_64-pc-linux-gnu-library/4.4
```

Optional arguments:

```bash
Rscript gmm_hpc_driver.R \
  --seed 1 \
  --loops 20 \
  --starts 120 \
  --target 100 \
  --tol 1e-7 \
  --queue serial \
  --wall-time 14:00 \
  --memory-gb 8 \
  --r-module R/4.4.0 \
  --r-libs-user "$HOME/R/x86_64-pc-linux-gnu-library/4.4" \
  --poll-seconds 60 \
  --max-wait-hours 72 \
  --out-dir hpc_runs/seed1
```

The driver submits each loop with an LSF script using `#BSUB` directives, matching the working `hpc_long_loop` style.

`--loops` is the maximum number of loops. If `--tol` is positive, the driver stops early when
`max(abs(beta_new - beta_old)) <= tol`.

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

Logs are grouped by loop:

```text
hpc_runs/seed1/logs/driver_seed1.log
hpc_runs/seed1/logs/loop_01_seed1.out
hpc_runs/seed1/logs/loop_01_seed1.err
```

No per-start `.rds` files are written.
