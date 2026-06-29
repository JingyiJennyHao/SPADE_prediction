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
cd /path/to/SPADE_prediction/hpc_long_loop/GMM
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
