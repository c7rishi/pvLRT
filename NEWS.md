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
