# version 0.0.0.9000
The first development version of the package



# version 0.0.1
* Changed internals of lrt_zi_poisson and lrt_vanilla_poisson:
  - Instead of storing random null tables and LRT statistics through lapply, they are now computed inside a for loop, and p-value is updated on the fly.
  
* Removed pbapply from import; added progress

# version 0.0.1.1
* updated estimation of omega's in lrt_zi_poisson, to reflect grouping structure among drugs, if present, as provided through drug_class_idx.

# version 0.1.1.2
* added option 'grouped_est_omega' (logical argument) in lrt_zi_poisson to specify whehter or not do grouped estimation of omega in extended LRT
* added a check in lrt_zi_poisson and lrt_vanilla_poisson to ensure that drug_class_idx has the same columns as the input contin_table

# version 0.0.2.0
* added options "test_omega", "pval_ineq_strict", "skip_null_omega_estimation", and "use_gamma_smooth_omega" to lrt_zi_poisson()
  - test_omega: performs pseudo LRT for the zi parameter omega
  - pval_ineq_strict: determines whether to use >= or > in the definition of the p-values of the test
  - skip_null_omega_estimation: determines whether or not to skip estimation of the omega on the null tables for lambda test. Defaults to TRUE, in which case the test effectively reduces to a parametric bootstrap pseudo-likelihood ratio test. NOTE: this changes the default behavior of the function. In previous versions the same null omegas were used in all null tables, which can currently be achieved by setting skip_null_omega = FALSE 
  - use_gamma_smooth_omega: determines whether or not to use a gamma prior (smoothing) while estimating the omegas. Defaults to FALSE. NOTE: this changes the default behavior of the function. In previous versions a gamma prior was always used, which in current version can be achieved by setting use_gamma_smooth_omega = FALSE 
