# frecipes (name will likely change to hydrorecipes as it is designed to replace it)

WARNING: This package is in early stages of development and is likely to change dramatically. Names of steps and thier implementations are likely to change.

This package is based on [recipes](https://recipes.tidymodels.org). The goals of the package are to increased the computational speed, decrease memory consumption, increase consistency between steps, and decrease some boiler plate code for step additions. The first three goals are likely to be achieved but the fourth might not. Part of this is due to providing both the *R6* and "standard" R ways to run the code.


It diverges in a few ways:

- based on [R6](https://r6.r-lib.org)
- focus on decreasing memory usage
- focus on speed
- steps tailored to groundwater applications
- decrease the number of dependencies and foreign functions
- more flexible output options (list, matrix, data.frame, data.table, tibble)
- statistically less robust 

- API changes
  - uses *terms* instead of *...* for variable selection and selections are wrapped
  in `c()` when more than one is required.
  - *R6* and standard R interfaces 
  - additional functions for pulling 



steps to do:

- step_temporary_deployment
- step_be_*
