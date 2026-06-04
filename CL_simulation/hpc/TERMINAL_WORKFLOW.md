# CL Simulation HPC Terminal Workflow

This workflow runs the `CL_simulation` algorithm on an LSF cluster. It does not modify `CL_simulation/composite_likelihood.R`; the runnable files are copied/refactored into this `hpc` folder.

## 1. Get the GitHub Folder on the Cluster

If the repository is not on the cluster yet:

```bash
git clone https://github.com/JingyiJennyHao/SPADE_prediction.git
cd SPADE_prediction/CL_simulation/hpc
```

If the repository already exists:

```bash
cd ~/SPADE_prediction
git pull
cd CL_simulation/hpc
```

Make the loop script executable:

```bash
chmod +x run_cl_loop.sh
```

## 2. Tiny Smoke Test

Use only 1 round, 2 beta starts, 2 weight starts, and a smaller `optim(maxit)` so the test does not run for hours:

```bash
ROUNDS=1 BETA_TASKS=2 WEIGHT_TASKS=2 OPTIM_MAXIT=1000 ./run_cl_loop.sh
```

After it starts, check jobs with:

```bash
bjobs
```

Expected outputs:

```text
results/round01/beta_results_round01.csv
results/round01/weight_results_round01.csv
results/round01/best_beta_round01.rds
results/round01/best_weight_round01.rds
results/round01/logs/
```

The beta CSV should have 2 result rows. The weight CSV should have 2 result rows.

## 3. Full Run

The production defaults are 20 rounds, 100 beta starts per round, 100 weight starts per round, and `optim(maxit=8000)`:

```bash
./run_cl_loop.sh
```

Equivalent explicit command:

```bash
ROUNDS=20 BETA_TASKS=100 WEIGHT_TASKS=100 OPTIM_MAXIT=8000 ./run_cl_loop.sh
```

## 4. Useful Overrides

Change the number of simulated hospitals/groups:

```bash
GROUP_COUNT=200 ./run_cl_loop.sh
```

Use a different result folder:

```bash
RESULTS_ROOT=$PWD/results_test ROUNDS=1 BETA_TASKS=2 WEIGHT_TASKS=2 OPTIM_MAXIT=1000 ./run_cl_loop.sh
```

Use a different R module or R library path:

```bash
R_MODULE=R/4.4 R_LIBS_USER_PATH=$HOME/R/x86_64-pc-linux-gnu-library/4.4 ./run_cl_loop.sh
```

Use a specific LSF queue:

```bash
QUEUE_NAME=normal ./run_cl_loop.sh
```

## 5. Inspect Results

For round 1:

```bash
ls results/round01
head results/round01/beta_results_round01.csv
head results/round01/weight_results_round01.csv
```

Read the selected best beta:

```bash
Rscript -e 'x <- readRDS("results/round01/best_beta_round01.rds"); print(x$objective_value); print(x$beta_hat)'
```

Read the selected best alpha and weights:

```bash
Rscript -e 'x <- readRDS("results/round01/best_weight_round01.rds"); print(x$objective_value); print(x$alpha_hat); print(x$weights_hat)'
```

Inspect logs:

```bash
ls results/round01/logs
tail results/round01/logs/*.err
```

## 6. Resume or Rerun

The loop skips a round stage if its summary file already exists:

```text
best_beta_roundXX.rds
best_weight_roundXX.rds
```

To resume an interrupted run, use the same command again:

```bash
./run_cl_loop.sh
```

To rerun a round from scratch, remove that round's result folder first:

```bash
rm -rf results/round05
./run_cl_loop.sh
```

Only remove result folders you intentionally want to recompute.
