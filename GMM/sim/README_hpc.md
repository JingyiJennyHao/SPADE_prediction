# GMM HPC Run

The original `gmm_sim.Rmd` is unchanged. The HPC run uses these new files:

- `gmm_core.R`: shared simulation and GMM functions copied from the notebook.
- `gmm_start_worker.R`: runs one start point and appends one row to a loop CSV.
- `gmm_start_worker_lsf.sh`: LSF array wrapper that maps `LSB_JOBINDEX` to `start_id`.
- `gmm_hpc_driver.R`: submits each loop, waits until enough completed start rows are recorded or no jobs are left, picks the best converged beta, then submits the next loop.
- `submit_controller_lsf.sh`: submits the driver/controller itself as an LSF job so it keeps polling and submitting later loops after logout.

## Run on an LSF cluster

From this directory:

```bash
./run_loop.sh
```

If you do not want to keep a tmux session open, submit the controller itself:

```bash
./submit_controller_lsf.sh
```

Default `run_loop.sh` settings are:

```text
SEED=1
LOOPS=20
START_LOOP=1
STARTS=120
TARGET=70
TARGET_TYPE=completed
TOL=1e-7
QUEUE=serial
WALL_TIME=14:00
MEMORY_GB=8
R_MODULE=R/4.4.0
R_LIBS_USER_PATH=$HOME/R/x86_64-pc-linux-gnu-library/4.4
REQUIRE_CONVERGED_BEST=1
```

Optional arguments:

```bash
Rscript gmm_hpc_driver.R \
  --seed 1 \
  --loops 20 \
  --starts 120 \
  --target 70 \
  --target-type completed \
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

`--loops` is the maximum number of loops. `--target` is the number of target rows needed to move to the next loop. By default `--target-type completed`, so target rows are any completed rows in the loop CSV. The next beta is still chosen only from `convergence == 0` rows unless `REQUIRE_CONVERGED_BEST=0` is set for debugging. If `--tol` is positive, the driver stops early when
`max(abs(beta_new - beta_old)) <= tol`.

The driver also proceeds when the submitted array job is no longer active, even if the target was not reached. This handles loops where fewer than `--target` starts write usable rows.

To resume after existing loop results, set `START_LOOP`:

```bash
START_LOOP=3 ./run_loop.sh
```

The driver checks prior loop CSVs and beta files, logs how many results already exist, and starts from the requested loop.

For a small HPC plumbing smoke test, submit 10 starts and move to the next loop after one completed result:

```bash
STARTS=10 TARGET=1 TARGET_TYPE=completed LOOPS=3 REQUIRE_CONVERGED_BEST=0 OUT_DIR=hpc_runs/smoke_seed1_10 ./run_loop.sh
```

To run the same smoke test without tmux:

```bash
STARTS=10 TARGET=1 TARGET_TYPE=completed LOOPS=3 REQUIRE_CONVERGED_BEST=0 OUT_DIR=hpc_runs/smoke_seed1_10 ./submit_controller_lsf.sh
```

`REQUIRE_CONVERGED_BEST=0` is only for debugging the submit/wait/resume procedure with any finite objective row. For the usual simulation, keep the default `REQUIRE_CONVERGED_BEST=1`.

## Outputs

For loop `k`, all 120 start-point jobs append to one file:

```text
hpc_runs/seed1/loop_0k_start_results_seed1.csv
```

The driver advances once at least `--target` rows of `--target-type` exist in that loop file, or once the submitted job is no longer active. It writes the best converged beta from each loop to:

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
