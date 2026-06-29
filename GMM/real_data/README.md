# GMM Estimation

This folder implements the GMM step based on the score vectors from `hpc_long_loop/composite_utils.R`.

For a supplied `beta_0`, the code:

1. Computes each row score `grad_beta ll_row(i, t, beta_0)`.
2. Stacks the time-specific scores into one vector `S_i` for each group `i`.
3. Estimates `V = cov(S_i)` across groups.
4. Sets `W = (V + ridge I)^-1`, with a generalized inverse fallback.
5. Tries 100 beta starts by default and keeps the fit minimizing `sum_i S_i(beta)^T W S_i(beta)`.

Example:

```sh
cd /path/to/SPADE_prediction/GMM
Rscript gmm_runner.R \
  --groups 200 \
  --seed 1 \
  --data /path/to/data.rds \
  --beta0-file /path/to/stage1_beta_summary.rds \
  --out /path/to/gmm_result.rds
```

Useful options:

- `--n-starts 100`: number of starting points.
- `--start-scale 0.1`: normal perturbation scale around `beta_0`.
- `--ridge 1e-6`: ridge added before inverting `V`.
- `--uncentered-v 1`: use `crossprod(S) / J` instead of the centered covariance.

## Iterative HPC workflow

Use `run_gmm_iterations.sh` on an LSF cluster when each GMM beta update should be one HPC iteration.
Each iteration:

1. Reads the previous iteration's best beta as `beta0`.
2. Submits an array job with 120 starts by default.
3. Proceeds once at least 100 task-level `.rds` results are available.
4. Cancels remaining array elements if the job is still active.
5. Saves the best converged beta as `summary/beta_hatXX.rds`.
6. Uses that best beta as `beta0` for the next iteration.

Example:

```sh
cd /path/to/SPADE_prediction/GMM
chmod +x submit_gmm_iteration.sh run_gmm_iterations.sh

INITIAL_BETA_FILE=/path/to/stage1_beta_summary.rds \
DATA_PATH=/path/to/data.rds \
GROUP_COUNT=200 \
MAX_ITERATIONS=20 \
./run_gmm_iterations.sh
```

Default controller settings:

- `TOTAL_TASKS=120`: submit 120 start-point tasks per iteration.
- `MIN_RESULTS=100`: advance once 100 result files exist.
- `START_SCALE=0.1`: perturb starts around the previous best beta.
- `BETA_TOL=0`: run until `MAX_ITERATIONS`; set a positive value to stop when the max absolute beta change is below that tolerance.
- `KEEP_RUN_RDS=0`: delete task-level RDS files after each iteration summary; set `KEEP_RUN_RDS=1` for debugging.

The controller log is:

```sh
results/iterations/gmm_iteration_controller.log
```

To resume at iteration 5 after iterations 1 through 4 already completed:

```sh
START_ITERATION=5 MAX_ITERATIONS=20 ./run_gmm_iterations.sh
```
