Algorithm

Step 0.1: let weight = w^{(0)}, (start with w^{(0) = (1/3, 1/3, 1/3)). Generating 100 starting points of beta, optimize the composite likelihood, to get beta^{(0)} in term of objective value (the objective function is the negative composite log likelihood, so minimize it) 

Step 0.2: Fixed beta^{(0)}, then generating 100 starting points of weights, optimize the composite likelihood to get beta^{(1)} in term of `objective_var_fun` which is A-criteria, minimizing sum of variance to get weight = w^{(1)}

Step i.1: With fixed w^{(i)}, to get beta^{(i)}

Step i.2: With fixed beta^{(i)}, to get w^{(i+1)}

