# Two-stage simulation with separate inference

## Confirmed workflow

One replication generates one dataset shared by both stages. Stage 1 submits
20 independent starts in an LSF array. Only after all Stage 1 jobs end does the
controller rank results and submit the Stage 2 candidate array.

### Stage 1

- Propose beta0[j] = true_beta[j] + Normal(0, sd[j]^2). Example sd is 2 for all
  14 parameters; each is configurable. Legacy uniform bounds remain supported.
- Pre-generate and save up to `start_max_attempts=50` candidates per start.
  Workers check finite beta, all scores, and Q(beta0; I42). There is NO
  covariance inversion in the initial check. Reject and try the next proposal
  only if this check fails; exhausted starts are excluded with diagnostics.
- First optim uses identity W. From optimization 2, compute W from the previous
  beta, hold W fixed during optim, and update beta. Every start has its own path.
- From optimization 2, stop when max(abs(beta_new-beta_old)) < `outer_tol=1e-5`.
  `loops=20` includes the first identity-W optimization. Each optim uses
  `maxit=8000` and the existing default relative tolerance. Later errors or
  max-loop nonconvergence do not trigger re-drawing starts.
- Filter on path convergence plus finite beta/value, then rank the literal last
  optim$value. Keep at most 10; fewer converged paths means fewer candidates.
  Inner optim status is recorded separately, not used as path convergence.
- Select the Stage 1 winner, recompute W at its final beta, and run/save the
  existing inference BEFORE dispatching Stage 2. Zero eligible starts stops
  with diagnostics and neither inference runs.

### Stage 2

Each selected candidate runs in its own array job. Initialize b at Stage 1 beta
and c=0. Effective theta=b+c*v_c; v_c has ones at 1,4,7,8,11,14. Alternate:

1. Fix b, optimize c over [-3,3].
2. Fix c, optimize b (Nelder-Mead, maxit=5000, reltol=1e-8).

EVERY evaluation uses F(theta)=Q(theta; W_hat(theta)), recomputing W from theta.
From alternating round 2, convergence is ONLY abs(c_new-c_old) < `cb_tol=1e-5`.
Beta/objective changes remain diagnostics, not stopping criteria. At most
`cb_max_iter=20` rounds; unconverged candidates are excluded from selection.
The existing c-search tolerance is unchanged.

Missing, corrupt, errored or nonfinite candidate results are logged and skipped.
Among c-converged finite candidates, select the smallest final F(theta), then
recompute W from this selected theta and run the same inference calculations.
If none qualify, mark Stage 2 failed; retain Stage 1 estimates and inference.
Scheduler submission/barrier errors also retain already saved Stage 1 results.

All active scores use score_ll_row_closed(); ridge=0 with no automatic fallback.
Inference retains its numerical Jacobian of the score. The DGP, moments,
objective formula and interval calculations are unchanged. Each stage uses the
same simulated data and the same starting RNG state for interval sampling;
inference does not advance the optimization's RNG state. Inference errors are
recorded separately from estimation failures, without inventing intervals.

## Submit one replication

Upload all changed code before submitting, and use a fresh output path. Do not
overwrite scripts being used by an older active controller/array.

```bash
cd ~/SPADE_prediction/GMM/sim
bash submit_two_stage_lsf.sh config_single_replication.R
bjobs
```

Edit R settings inside `config_single_replication.R`, not directly in bash.
The controller holds one CPU while waiting; child arrays and barrier jobs also
require slots. Compute nodes must permit bsub, and files must be on shared
storage. Actual concurrency depends on scheduler resources. A killed controller
does not automatically cancel children; array/barrier IDs are saved in logs.
The direct `Rscript run_two_stage_sim.R config_single_replication.R` entry uses
sequential execution unless GMM_STAGE1_BACKEND=lsf is set (it controls both stages).
Only static syntax/interface checks have been performed; no simulations or jobs
were run during these edits.

## Output files

Set `out_dir = "results/seed1"` in the config. One directory per replication:

```text
results/seed1/
  stage1.csv
  stage2.csv
  comparison.csv
  status.csv
  diagnostics.rds
  starts/
  refinements/
```

Stage 1 inference is saved before Stage 2. Stage 2 inference is written only on
success. Both retain the existing beta/variance/coverage and 17 PI coverage
columns. `comparison.csv` contains estimates, errors, objective definitions,
timings and statuses; `status.csv` records estimation/inference failures.
`diagnostics.rds` contains checkpoints, proposals/rejections, full optimization
histories, selections and final weighting matrices. Stage 1 failure is recorded
here before stopping. Final summary CSVs require the controller to finish.

Each worker directory contains a unique dispatch subdirectory with input.rds,
submission.txt, individual result RDS and logs. Missing or invalid results are
listed in unavailable_results.csv. Controller logs remain in
`logs/two_stage_<jobID>.out` and `.err`.

Legacy `out_csv` configs still work: their containing directory becomes out_dir;
the old basename is ignored. Explicit out_dir takes precedence. Old results are
not renamed. Use a fresh directory for each run: inference appends, summary and
diagnostic files overwrite. A custom diagnostics_file also determines the parent
of worker directories. Old analysis scripts matching inference_seedX.csv require
adjustment for stage1.csv/stage2.csv.

## Inspect accepted initial points

input.rds$initials contains FIRST proposals. Actual accepted beta0 is
`starts[[id]]$initial` in diagnostics or `$initial` in each start result.
`accepted_attempt` and `initial_attempts` retain attempts/rejections; no accepted
point or unavailable job has NA initial beta. Candidate pools are generated
centrally; identical seed/config reproduces them independently of job order.
Accepted points follow the proposal distribution conditional on passing checks.

Run R code below from a Stage 1 result directory (after loading R):

```r
files <- sort(list.files(pattern = "^start_[0-9]{4}\\.rds$"))
rows <- lapply(files, function(file) {
  x <- readRDS(file)
  row <- data.frame(start_id=x$start_id, accepted_attempt=x$accepted_attempt,
                    stop_reason=x$stop_reason)
  row[paste0("b", 1:14)] <- as.list(x$initial)
  row
})
if (length(rows)) print(do.call(rbind, rows), row.names=FALSE)
```

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
