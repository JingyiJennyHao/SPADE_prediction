## Algorithm

Step 0.1: let weight = w^{(0)}, (start with w^{(0) = (1/3, 1/3, 1/3)). Generating 100 starting points of beta, optimize the composite likelihood, to get beta^{(0)} in term of objective value (the objective function is the negative composite log likelihood, so minimize it) 

Step 0.2: Fixed beta^{(0)}, then generating 100 starting points of weights, optimize the composite likelihood to get beta^{(1)} in term of `objective_var_fun` which is A-criteria, minimizing sum of variance to get weight = w^{(1)}

Step i.1: With fixed w^{(i)}, to get beta^{(i)}

Step i.2: With fixed beta^{(i)}, to get w^{(i+1)}

## simulation 1: beta convergence

`revised_composite_likelihood.R`: two loops

Submit the job:

``` bash
J=500
START=1
END=100
RESULT_ROOT="$HOME/SPADE_prediction_results"
OUT_CSV="${RESULT_ROOT}/results/sim_results_J${J}_tasks${START}_${END}.csv"

bsub -J "cl_sim[${START}-${END}]" \
  -W 60:00 \
  -n 1 \
  -R "rusage[mem=8GB]" \
  -o "${RESULT_ROOT}/logs/cl_sim_J${J}_tasks${START}_${END}_%J.out" \
  -e "${RESULT_ROOT}/logs/cl_sim_J${J}_tasks${START}_${END}_%J.err" \
  "module load R/4.4 && export R_LIBS_USER=$HOME/R/x86_64-pc-linux-gnu-library/4.4 && Rscript revised_composite_likelihood.R ${OUT_CSV} \$((1000 + LSB_JOBINDEX)) ${J}"
```
## simulation 2: weights convergence

All files are in the folder hpc


