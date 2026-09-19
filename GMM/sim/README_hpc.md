# Two-stage simulation: one replication on the cluster

The current workflow is:

1. Simulate one dataset and pre-generate a candidate pool for each start.
   Workers accept the first numerically valid initial vector from their pool.
2. Stage 1 array: each start independently iterates
   `beta_k -> score(beta_k) -> W(beta_k) -> optim with that W -> beta_(k+1)`.
   Check `max(abs(beta_(k+1) - beta_k)) < outer_tol` after each step.
3. Rank the converged paths by their final `optim$value`; retain at most 10.
4. Stage 2 array: independently refine each candidate by alternating c and beta.
   Every objective evaluation recomputes W at the effective beta, exactly as
   in the real-data `Q_beta_c` implementation.
5. Select the smallest final Stage 2 objective; calculate final W and inference.

There is no shared-W outer loop around the two stages, and no identity-W
initialization. The moment functions, DGP, true-beta defaults, numerical scores,
W estimator (including its existing ridge), and inference calculations are
supplied by the existing `gmm_sim.Rmd` and remain unchanged.

## Submit one replication

Edit `config_single_replication.R`, then run from the cluster's sim directory:

```bash
bash submit_two_stage_lsf.sh config_single_replication.R
bjobs
```

This submits one controller, which submits one Stage 1 array and later one
Stage 2 array. Both arrays use one CPU per element. Actual concurrency depends
on resources and account limits. A dependent barrier waits for each entire
array to end. The controller holds one CPU while waiting; its wall time covers
both stages and queue waits. Compute nodes must permit `bsub`, the barrier needs
a scheduler slot, and all jobs must access the same shared directory.

For sequential execution inside an allocated session:

```bash
Rscript run_two_stage_sim.R config_single_replication.R
```

Only static syntax checks have been performed locally. No simulations or LSF
jobs have been run as part of these edits. Older already-running jobs do not
switch algorithms. Finish or stop the old controller and its child arrays/
barriers before replacing their files, then submit with a fresh output path.
Stopping a controller does not automatically cancel its children; submission
IDs are retained in the per-stage directories.

## Configuration and convergence

- `true_beta`: 14-element DGP truth, also used for beta errors and coverage.
- `start_sd`: scalar or 14 positive finite standard deviations. The example sets
  `rep(2, 14)`; change each coordinate independently if desired. Without bounds,
  proposals are `true_beta + rnorm(14, 0, start_sd)`.
- `start_lower`, `start_upper`: optional 14-element absolute uniform bounds for
  legacy configurations. Supplying both selects uniform proposals instead of sd.
- `start_max_attempts`: positive integer, default 50 including the first proposal.
  Before dispatch, all candidate pools are generated and saved. Identical seed
  and configuration reproduce proposals independently of scheduling order.
  Changing the attempt limit can change later starts' random draws.
  Initial beta, all observation scores, W and initial objective must be finite
  and computable. Reject an invalid proposal and try the next one; reuse the
  accepted W for the first optimization. Exhaustion records
  `initial_checks_exhausted` and excludes the start. Later path errors or
  nonconvergence do not trigger resampling. No model/objective changes are made.
  Accepted starts follow the proposal distribution conditional on passing checks;
  rejection histories are retained for interpretation of experiments.
- `n_starts`: number of independent Stage 1 paths, default 20.
- `loops`: maximum beta-W iterations **per Stage 1 path**, default 20.
- `outer_tol`: path convergence threshold, default 1e-5 (absolute max beta change).
- `maxit`: maximum iterations for each Stage 1 `optim`, default 8000; its
  default relative objective tolerance is about 1.49e-8.
- `cb_max_iter`: maximum c-beta alternating rounds, default 20.
- `cb_tol`: all three absolute changes (raw beta max norm, c, objective) must be
  below this threshold, default 1e-5, starting from round 2.
- `cb_beta_maxit`: maximum iterations of each Stage 2 beta optimization,
  default 5000; `reltol` remains 1e-8.
- `c_interval`: c search interval, default [-3, 3]. The real-data default
  `optimize` position tolerance (about 1.22e-4) is preserved.

Stage 1 eligibility uses path convergence and finite final beta/value, not
merely a single `optim$convergence == 0`. Inner optimizer codes are recorded
separately. Paths that reach `loops` without meeting `outer_tol` are excluded;
1-9 converged paths are all refined, and zero converged paths stops with
saved diagnostics. Stage 1 ranking uses the literal last `optim$value`, whose
W was computed at the beta *before* that optimization; it is not recomputed
for ranking.

Stage 2 initializes raw beta at the candidate and c at zero. Its effective
beta is `raw_beta + c*v_c`, with direction indices 1, 4, 7, 8, 11, 14.
Each evaluation uses `F(theta) = Q(theta; W_hat(theta))`. A finite result at
the alternating cap still participates in final selection, with convergence
false. A refinement error/nonfinite result stops with diagnostics; missing
Stage 2 worker output files also stop aggregation. Stage 1 is tolerant of
missing, unreadable or malformed result files after the entire array ends:
these starts are logged and excluded, with original IDs preserved. Remaining
converged paths proceed to Stage 2; zero eligible paths still stops.

## Reading actual accepted initial points

`input.rds$initials` and the main diagnostic `initials` contain FIRST proposals,
not necessarily accepted starts. `candidate_pools[[start_id]][[attempt]]` stores
all proposals. The Stage 1 input also saves `start_settings`.

After workers finish, run this R code from their result directory:

```r
files <- sort(list.files(pattern = "^start_[0-9]{4}\\.rds$"))
rows <- lapply(files, function(file) {
  x <- readRDS(file)
  row <- data.frame(start_id = x$start_id,
    accepted_attempt = x$accepted_attempt,
    attempts = length(x$initial_attempts), stop_reason = x$stop_reason)
  row[paste0("b", 1:14)] <- as.list(x$initial)
  row
})
if (length(rows)) print(do.call(rbind, rows), row.names = FALSE)
# Detailed rejection reasons and proposed vectors for start 1:
x <- readRDS("start_0001.rds")
x$initial_attempts
```

Actual accepted beta is `starts[[id]]$initial` in the main diagnostics. No accepted
point is represented by 14 NAs; missing job results also have unknown/NA initial
points. Acceptance is not a claim that the subsequent optimization will converge.

## Logs and diagnostics

Controller log: `logs/two_stage_<jobID>.out` (and `.err`).
It reports array IDs, converged path count and the selected candidate.

Stage 1: `<diagnostics_file>.starts/loop_01_<unique-id>/` contains `input.rds`,
`submission.txt`, per-start `.out`/`.err` and `start_0001.rds` etc. Each start
log prints its beta-W iteration, objective, beta change, inner optimizer code,
and path convergence. Missing/unreadable/invalid outputs are listed in `unavailable_results.csv`.
The controller logs submitted, received, unavailable, failed and converged counts.
`failed` includes unavailable outputs plus returned optimizer/path errors;
nonconverged paths without errors are not counted as failed.
The directory's `loop_01` is a historical label for the
single array dispatch, not the start's internal iteration number.

Stage 2: `<diagnostics_file>.refinements/loop_01_<unique-id>/` contains input,
submission IDs, per-candidate logs and `candidate_0001.rds` etc. Candidate
indices are Stage 1 ranks and results retain the original start ID. Each
candidate log prints a completion summary; full alternating history is in RDS.

For `out_csv = "single_replication/inference_seed1.csv"`:

- Inference CSV: beta estimates, variances, beta coverage, all 17 original PI
  coverage fields. `outer_converged` describes the winning candidate's Stage 1
  path, `final_loop` its Stage 1 iteration count, `cb_converged` its Stage 2
  alternating convergence, and `beta_conv` its final Stage 2 inner optimizer code.
- `.stages.csv`: two rows total, Stage 1 winner and Stage 2 winner, with start
  IDs, original objectives, objective definitions, dynamic objectives F(beta),
  dynamic-objective improvement, elapsed seconds (including scheduler waits),
  iteration counts, full betas, signed errors, L2 errors and candidate counts.
  The literal Stage 1 and Stage 2 objective columns have different W semantics;
  compare `dynamic_objective` for a common function. This extra evaluation does
  not alter Stage 1 ranking. It is not a separate Stage-1-only coverage run.
- `.diagnostics.rds`: settings, initial vectors, `starts` (full independent
  paths with each W_used, beta, objective, beta change, inner status, and final
  path convergence), `unavailable_start_ids`, `failed_start_ids`, `n_received`,
  `top_start_ids`, `stage1`, `refined` (original objective,
  initial dynamic objective, final objective/c/raw and effective beta, iteration
  count, convergence/status, full alternating history), selected candidate rank,
  `stage2`, comparison table, final W. Checkpoints survive later inference errors.

Inference CSV appends; diagnostic and stage files overwrite. Use fresh output
paths for each run. Rendering the Rmd does not launch its example.

The legacy scripts below are not called by this entry and do not implement
this two-stage workflow.

---

# Legacy per-start cluster workflow (not the two-stage entry)

# Submit GMM simulations on the cluster

The HPC workflow now uses one independent LSF task for one `(seed,
start point)` pair. Each task runs up to `LOOPS` sequentially, using its own
optimized beta as the starting beta for the next loop. There is no
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

Each task runs up to:

```text
start point -> loop 1 -> loop 2 -> ... -> loop LOOPS
```

The path can stop early from loop 2 onward when
`max(abs(beta_new - beta_old)) < OUTER_TOL`.

Loop 1 starts from a seed- and start-specific random perturbation of
`true_beta`. Later loops start from that same task's previous optimized beta;
different start points do not select or share betas during execution.

After all start-point tasks for a seed finish, one dependent finalizer job
reads the final loop CSV, selects the row with the smallest finite objective
value across start points, and uses that row's beta for inference. The selected
row is saved as `final_selection_seedX.csv`, and inference is saved as
`inference_seedX.csv`.

The finalizer runs after the array ends even if one or more start-point tasks
exit because of an LSF or other job-level failure. It uses all available rows
with finite beta and objective values; it stops only if no valid final beta is
available.

## Settings

- `SEED_FIRST`, `SEED_LAST`: simulation seed range.
- `J`: number of clusters/hospitals in each simulation.
- `LOOPS`: number of GMM weight-update loops run by every task.
- `OUTER_TOL`: outer-loop convergence threshold for
  `max(abs(beta_new - beta_old))`; default `1e-5`. The check starts at loop 2.
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
  final_start_results_seed1.csv
  final_selection_seed1.csv
  inference_seed1.csv
  ...
  logs/start_1.out
  logs/start_1.err
```

Every loop CSV contains one row per completed start point. The `start_id`
column identifies the task that produced the row. Each start point also writes
one final row to `final_start_results_seed1.csv`; this row records its actual
final loop, whether outer convergence was reached, and its final objective
value. `final_selection_seed1.csv` records the minimum-objective start point
and the beta used for inference.

Monitor tasks with:

```bash
bjobs
bjobs -l
tail -f hpc_runs/J200_seeds1_10_starts20_loops20/seed1/logs/start_1.out
```
