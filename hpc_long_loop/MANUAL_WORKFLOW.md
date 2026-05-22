# Manual HPC Workflow

This setup is for running each stage by hand while keeping review easy.

## What you get
- One shared `.out` and one shared `.err` file per submitted stage batch
- One `.csv` file that combines all task results for that batch
- You can inspect the CSV, choose the best row, and then decide the next stage yourself

## Stage 1
Example: 100 starts, 50 random and 50 fixed, all hospitals.

```bash
cd ~/SPADE/hpc_long_loop
chmod +x submit_stage1_manual.sh submit_stage2_manual.sh
GROUP_COUNT=210 TOTAL_TASKS=100 RANDOM_TASKS=50 JOB_NAME=s1_round1 ./submit_stage1_manual.sh
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
