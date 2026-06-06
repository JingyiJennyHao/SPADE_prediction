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

Use only 1 round, 2 beta starts, 2 weight starts, and a smaller `optim(maxit)` so the test does not run for hours. The minimum-result thresholds default to the number of submitted tasks when fewer than 100 tasks are submitted:

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
results/run_cl_loop_controller.log
```

The beta CSV should have 2 result rows. The weight CSV should have 2 result rows.

## 3. Full Run

The production defaults are 20 rounds, 130 beta starts per round, 130 weight starts per round, minimum 100 finished rows per stage, and `optim(maxit=8000)`. The controller moves forward as soon as each stage has enough CSV rows and cancels any remaining array elements.

```bash
./run_cl_loop.sh
```

Equivalent explicit command:

```bash
ROUNDS=20 BETA_TASKS=130 BETA_MIN_RESULTS=100 WEIGHT_TASKS=130 WEIGHT_MIN_RESULTS=100 OPTIM_MAXIT=8000 ./run_cl_loop.sh
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

Submit 130 starts but move forward once 100 results are written:

```bash
BETA_TASKS=130 BETA_MIN_RESULTS=100 WEIGHT_TASKS=130 WEIGHT_MIN_RESULTS=100 ./run_cl_loop.sh
```

Change how often the controller reports progress:

```bash
POLL_SECONDS=30 ./run_cl_loop.sh
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

Each stage now has one shared `.out` and one shared `.err` per array job, not one file per array element:

```text
results/round01/logs/beta_<jobid>.out
results/round01/logs/beta_<jobid>.err
results/round01/logs/weight_<jobid>.out
results/round01/logs/weight_<jobid>.err
```

The controller progress log is:

```bash
tail -f results/run_cl_loop_controller.log
```

Example controller messages:

```text
Round 1 beta-step: 72/130 results written; need 100 to move forward.
Round 1 beta-step: get result threshold reached (100/130); moving forward.
Round 1: summarizing beta results from results/round01/beta_results_round01.csv.
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

## 7. Save weight and beta result

```bash
module load R/4.4
cd ~/SPADE_prediction/CL_simulation/hpc
Rscript -e '
files <- sort(Sys.glob("results/round*/best_beta_round*.rds"))
x <- lapply(files, readRDS)
names(x) <- basename(files)
saveRDS(x, "results/all_best_beta_results.rds")
cat("Saved", length(x), "beta RDS files to results/all_best_beta_results.rds\n")
'
```

```bash
Rscript -e '
files <- sort(Sys.glob("results/round*/best_weight_round*.rds"))
x <- lapply(files, readRDS)
names(x) <- basename(files)
saveRDS(x, "results/all_best_weight_results.rds")
cat("Saved", length(x), "weight RDS files to results/all_best_weight_results.rds\n")
'
```
