# Manual HPC Workflow

This setup is for running each stage by hand while keeping review easy.

## What you get
- One shared `.out` and one shared `.err` file per submitted stage batch
- One `.csv` file that combines all task results for that batch
- You can inspect the CSV, choose the best row, and then decide the next stage yourself

## Automated 20-Round Loop
Use `run_loop.sh` when you want the full beta-weight loop to run automatically.
This command runs up to 20 rounds and stops early if the weight convergence tolerance is reached.

```bash
cd ~/SPADE/hpc_long_loop
chmod +x run_loop.sh submit_stage1_manual.sh submit_stage2_manual.sh
MAX_ROUNDS=20 ./run_loop.sh
```

Default automated-loop settings:
- each beta step submits 130 starts
- the first 30 beta starts are random
- the remaining 100 beta starts use the previous round's beta when available
- round 1 uses fixed starts for those 100 non-random starts
- the controller proceeds after 100 beta result files are available
- each weight step submits 15 starts and proceeds after 10 result files are available

The main controller log is:

```bash
results/loop/run_loop_controller.log
```

Watch it while the loop is running:

```bash
tail -f results/loop/run_loop_controller.log
```

Each round also has stage logs such as:

```bash
results/loop/round01/stage1/logs/s1_round01.out
results/loop/round01/stage1/logs/s1_round01.err
results/loop/round01/stage2/logs/s2_round01.out
results/loop/round01/stage2/logs/s2_round01.err
```

To resume from a later round, keep the final maximum round number and set `START_ROUND`.
For example, resume at round 5 and run through at most round 20:

```bash
cd ~/SPADE/hpc_long_loop
START_ROUND=5 MAX_ROUNDS=20 ./run_loop.sh
```

The automated loop keeps summary CSV files plus compact `beta_hatXX.rds` and `weights_hatXX.rds` files.
It deletes task-level `runs/task_*.rds` files after each stage is summarized.
If you are debugging and want to keep every task-level RDS file:

```bash
cd ~/SPADE/hpc_long_loop
KEEP_RUN_RDS=1 MAX_ROUNDS=20 ./run_loop.sh
```

## Stage 1
Example: 130 starts, with 30 random starts and 100 fixed starts, all hospitals.
When `BETA_REF_FILE` is not set, the non-random starts use the fixed starting value.

```bash
cd ~/SPADE/hpc_long_loop
chmod +x submit_stage1_manual.sh submit_stage2_manual.sh
GROUP_COUNT=210 TOTAL_TASKS=130 RANDOM_TASKS=30 JOB_NAME=s1_round1 ./submit_stage1_manual.sh
```

This writes:
- result files under `results/manual_stage1/runs/`
- one shared log at `results/manual_stage1/logs/s1_round1.out`
- one shared err at `results/manual_stage1/logs/s1_round1.err`

After the jobs finish, collect one review CSV:

```bash
cd ~/SPADE/hpc_long_loop
module load R/4.4
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.4"
Rscript collect_stage1_results.R \
  --input-dir results/manual_stage1/runs \
  --out-csv results/manual_stage1/summary/s1_round1.csv
```

If you want the best beta saved as `.rds` too:

```bash
Rscript stage1_summarize.R \
  --input-dir results/manual_stage1/runs \
  --out results/manual_stage1/summary/beta_hat.rds
```

For a later beta round, use the previous round's beta and weights. The first 30 tasks are still random starts, and the other 100 tasks start near the previous beta estimate:

```bash
cd ~/SPADE/hpc_long_loop
RESULTS_DIR=$HOME/SPADE/hpc_long_loop/results/manual_stage1_round2 \
BETA_REF_FILE=$HOME/SPADE/hpc_long_loop/results/manual_stage1/summary/beta_hat.rds \
WEIGHTS_FILE=$HOME/SPADE/hpc_long_loop/results/manual_stage2/summary/weights_hat.rds \
GROUP_COUNT=210 TOTAL_TASKS=130 RANDOM_TASKS=30 JOB_NAME=s1_round2 ./submit_stage1_manual.sh
```

## Stage 2
Example: 10 starts using a chosen `beta_hat.rds`.

```bash
cd ~/SPADE/hpc_long_loop
BETA_FILE=$HOME/SPADE/hpc_long_loop/results/manual_stage1/summary/beta_hat.rds \
GROUP_COUNT=210 TOTAL_TASKS=10 JOB_NAME=s2_round1 ./submit_stage2_manual.sh
```

Then collect one review CSV:

```bash
cd ~/SPADE/hpc_long_loop
module load R/4.4
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.4"
Rscript collect_stage2_results.R \
  --input-dir results/manual_stage2/runs \
  --out-csv results/manual_stage2/summary/s2_round1.csv
```

If you want the best weights saved as `.rds` too:

```bash
Rscript stage2_summarize.R \
  --input-dir results/manual_stage2/runs \
  --out results/manual_stage2/summary/weights_hat.rds
```

## Notes
- Shared logs are easier to find, but array-task messages will be interleaved.
- The CSV is the main artifact for comparing all starts in one place.
- If you want separate folders per round, override `RESULTS_DIR=...` when you submit.
- The automated `run_loop.sh` deletes task-level `runs/task_*.rds` files after the CSV and compact summary RDS are created.
- Manual runs keep task-level RDS files until you delete them yourself, for example `rm -f results/manual_stage1/runs/task_*.rds`.
- Use `KEEP_RUN_RDS=1 ./run_loop.sh` when debugging and you want the automated loop to keep task-level RDS files.
