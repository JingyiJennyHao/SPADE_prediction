# GMM HPC Run

The original `gmm_sim.Rmd` is unchanged. The HPC run uses these new files:

- `gmm_core.R`: shared simulation and GMM functions copied from the notebook.
- `gmm_start_worker.R`: runs one start point and appends one row to a loop CSV.
- `gmm_start_worker_lsf.sh`: LSF array wrapper that maps `LSB_JOBINDEX` to `start_id`.
- `gmm_hpc_driver.R`: submits loop 1 until enough converged start rows are recorded, then submits one exact-beta optimization for each later loop.
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
J=200
LOOPS=6
START_LOOP=1
STARTS=10
TARGET=2
TARGET_TYPE=converged
TOL=1e-7
QUEUE=serial
WALL_TIME=14:00
MEMORY_GB=8
R_MODULE=R/4.4.0
R_LIBS_USER_PATH=$HOME/R/x86_64-pc-linux-gnu-library/4.4
START_SD=0.1
MAXIT=8000
REQUIRE_CONVERGED_BEST=1
RUN_INFERENCE=1
```

The optimization procedure is now:

- Loop 1 submits `STARTS` independent perturbed start points and waits for
  `TARGET` converged results. The defaults are 10 starts and 2 converged
  results.
- Each later loop submits exactly one optimization, initialized at the beta
  selected from the previous loop. No random perturbation is added in loops 2+
  (the `START_SD` setting therefore only affects loop 1).
- From loop 2 onward, the optimizer does not need `convergence == 0`: any
  result with a finite objective and finite beta estimate can provide the
  starting point for the next loop.
- `MAXIT` is passed to `optim(control = list(maxit = MAXIT))` for every loop.

For example, to use 80 loop-1 starts, require 30 converged results, and allow
12,000 optimizer iterations:

```bash
STARTS=80 TARGET=30 MAXIT=12000 ./run_loop.sh
```

The equivalent R options are `--starts 80 --target 30 --maxit 12000`.

Optional arguments:

```bash
Rscript gmm_hpc_driver.R \
  --seed 1 \
  --J 200 \
  --loops 6 \
  --starts 10 \
  --target 2 \
  --target-type completed \
  --tol 1e-7 \
  --start-sd 0.1 \
  --maxit 8000 \
  --queue serial \
  --wall-time 14:00 \
  --memory-gb 8 \
  --r-module R/4.4.0 \
  --r-libs-user "$HOME/R/x86_64-pc-linux-gnu-library/4.4" \
  --poll-seconds 60 \
  --max-wait-hours 72 \
  --out-dir hpc_runs/seed1
```

After the final loop (or tolerance-based early stop), the driver uses that
loop's selected beta to run the existing beta and prediction inference. It
writes `inference_seed<seed>.csv` in the seed output directory.

To submit seeds 1 through 10 with the default `J=200`, 10 loop-1 starts, and
advance after 2 converged results:

```bash
./submit_10_seeds_lsf.sh
```

The driver submits each loop with an LSF script using `#BSUB` directives, matching the working `hpc_long_loop` style.

`--loops` is the maximum number of loops. `--starts` and `--target` control
only loop 1: the number of random starts and the number of converged rows
needed to move to loop 2. Later loops always use one start and require only a
finite beta estimate and objective value; their `convergence` code is not used
as a gate. The next beta from loop 1 is selected from `convergence == 0` rows
unless `REQUIRE_CONVERGED_BEST=0` is set for debugging. The next beta from loops
2+ is selected from all usable finite beta rows. If `--tol` is positive, the driver stops early when
`max(abs(beta_new - beta_old)) <= tol`.

The driver stops with an error if loop 1 finishes before its required target is
reached, or if a later loop finishes without producing a usable beta estimate.
This prevents the chain from advancing without a beta value to pass forward.

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

`REQUIRE_CONVERGED_BEST=0` affects only loop 1 and is useful for debugging the
submit/wait/resume procedure with any finite objective row. Later loops always
accept a usable finite beta estimate. For the usual simulation, keep the
default `REQUIRE_CONVERGED_BEST=1`.

## Outputs

For loop 1, the `STARTS` start-point jobs append to one file. Each later loop
has one start-point job and appends to its own file:

```text
hpc_runs/seed1/loop_0k_start_results_seed1.csv
```

The driver advances once the required target exists in the loop file and
writes the selected beta from each loop to:

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
