# Submit GMM simulations on the cluster

The HPC workflow now uses one independent LSF task for one `(seed,
start point)` pair. Each task runs loop 1 through `LOOPS` sequentially, using
its own optimized beta as the starting beta for the next loop. There is no
controller job coordinating the loops; only a dependent seed-level finalizer
runs after the start-point tasks finish.

Relevant files:

- `gmm_core.R`: shared simulation, GMM, optimization, and inference functions.
- `gmm_start_worker.R`: runs all loops for one seed and one start point.
- `gmm_finalize_seed.R`: selects the lowest-objective final beta and runs
  inference after the start-point jobs finish.
- `submit_controller_lsf.sh`: submits the start-point array for one seed. The
  historical filename is retained for compatibility.
- `submit_10_seeds_lsf.sh`: submits one start-point array per seed.
- `gmm_start_worker_lsf.sh`: runs one start-to-finish task directly.

## Submit jobs

```bash
cd /home/jhao3/SPADE_prediction/GMM/sim
chmod u+x submit_10_seeds_lsf.sh submit_controller_lsf.sh
```

For 20 start points per seed and 20 loops:

```bash
SEED_FIRST=1 \
SEED_LAST=10 \
J=200 \
LOOPS=20 \
STARTS=20 \
START_SD=0.1 \
MAXIT=8000 \
RUN_ROOT=hpc_runs/J200_seeds1_10_starts20_loops20 \
./submit_10_seeds_lsf.sh
```

The defaults are `STARTS=20` and `LOOPS=6`. `STARTS` can be changed, but the
requested workflow is 20 start points per seed.

For each seed, the submission is one LSF array with one task per start point:

```text
gmm_s1[1-20]
gmm_s2[1-20]
...
```

Each task runs:

```text
start point -> loop 1 -> loop 2 -> ... -> loop LOOPS
```

Loop 1 starts from a seed- and start-specific random perturbation of
`true_beta`. Later loops start from that same task's previous optimized beta;
different start points do not select or share betas during execution.

After all start-point tasks for a seed finish, one dependent finalizer job
reads the final loop CSV, selects the row with the smallest finite objective
value across start points, and uses that row's beta for inference. The selected
row is saved as `final_selection_seedX.csv`, and inference is saved as
`inference_seedX.csv`.

## Settings

- `SEED_FIRST`, `SEED_LAST`: simulation seed range.
- `J`: number of clusters/hospitals in each simulation.
- `LOOPS`: number of GMM weight-update loops run by every task.
- `STARTS`: number of start-point tasks per seed; default `20`.
- `START_SD`: random perturbation size for loop 1.
- `MAXIT`: maximum `optim()` iterations per loop.
- `RUN_INFERENCE`: run final inference after selecting the lowest-objective
  beta; default `1`.
- `INFERENCE_PROB`: inference interval probability; default `0.95`.
- `RUN_ROOT`: output directory for the experiment. Use a new root for a new
  experiment.

## Results

Each seed has its own folder, for example:

```text
hpc_runs/J200_seeds1_10_starts20_loops20/seed1/
  loop_01_start_results_seed1.csv
  loop_02_start_results_seed1.csv
  ...
  loop_20_start_results_seed1.csv
  final_selection_seed1.csv
  inference_seed1.csv
  ...
  logs/start_1.out
  logs/start_1.err
```

Every loop CSV contains one row per completed start point. The `start_id`
column identifies the task that produced the row. `final_selection_seed1.csv`
records the minimum-objective start point and the beta used for inference.

Monitor tasks with:

```bash
bjobs
bjobs -l
tail -f hpc_runs/J200_seeds1_10_starts20_loops20/seed1/logs/start_1.out
```
