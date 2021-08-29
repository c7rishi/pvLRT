# version 0.1.4
safer handling of posterior structural zi probability in 0/0 form


#  version 0.1.3
Defaults pseudo lambda LRT now uses use_cont


# version 0.1.2

Lots of changes. Primarily analysis is now centered around 
a pvlrt object that is created when the function pvlrt() is 
used. Various summary, plotting, printing methods for pvlrt 
are created.

# version 0.1.1

* Changes in lrt_zi_poisson():
 - safer null statistics evaluation via tryCatch()
 - better print and summary methods for pvlrt objects 
 - functions for extracting pvalues, lrstat, zi probabilities from pvlrt object. (Documentation needs to be added.) 

# version 0.1.0

* Changes in lrt_zi_poisson():

  - removes the argument skip_null_omega_estimation. Null omega estimation is not performed as it is not necessary.


# version 0.0.2.3

* Changes in lrt_zi_poisson():

  - used analytic expression for pseudo lrt statistic for lambda, and a reduced algebraic expression for llik of omega. Both speed up computations.

  - Replaced omega_est_vec by omega_vec. omega_est_vec is now deprecated: will throw a warning if used instead of omega_vec.

  - Added the argument test_drug_idx to specify which columns to test.




# version 0.0.2.2

- add a secondary optim() call in gamma smoothed omega estimation function



# version 0.0.2.1

- add in the 10 missing rows (with names containing string 'total') in statin


# version 0.0.2.0
* added options "test_omega", "pval_ineq_strict", "skip_null_omega_estimation", and "use_gamma_smooth_omega" to lrt_zi_poisson()
  - test_omega: performs pseudo LRT for the zi parameter omega
  - pval_ineq_strict: determines whether to use >= or > in the definition of the p-values of the test
  - skip_null_omega_estimation: determines whether or not to skip estimation of the omega on the null tables for lambda test. Defaults to TRUE, in which case the test effectively reduces to a parametric bootstrap pseudo-likelihood ratio test. NOTE: this changes the default behavior of the function. In previous versions the same null omegas were used in all null tables, which can currently be achieved by setting skip_null_omega = FALSE 
  - use_gamma_smooth_omega: determines whether or not to use a gamma prior (smoothing) while estimating the omegas. Defaults to FALSE. NOTE: this changes the default behavior of the function. In previous versions a gamma prior was always used, which in current version can be achieved by setting use_gamma_smooth_omega = FALSE 

* added 4 new data sets: statin (all pts), statin1491, statin46 and gbca

* added function r_contin_table_zip() to  generate random contingency tables for adverse effect (across rows) and drug (across columns) combinations given row and column marginal totals, embedded signals, and possibly with zero inflation



# version 0.1.1.2
* added option 'grouped_est_omega' (logical argument) in lrt_zi_poisson to specify whehter or not do grouped estimation of omega in extended LRT
* added a check in lrt_zi_poisson and lrt_vanilla_poisson to ensure that drug_class_idx has the same columns as the input contin_table



# version 0.0.1.1
* updated estimation of omega's in lrt_zi_poisson, to reflect grouping structure among drugs, if present, as provided through drug_class_idx.


# version 0.0.1
* Changed internals of lrt_zi_poisson and lrt_vanilla_poisson:
  - Instead of storing random null tables and LRT statistics through lapply, they are now computed inside a for loop, and p-value is updated on the fly.
  
* Removed pbapply from import; added progress


# version 0.0.0.9000
The first development version of the package


