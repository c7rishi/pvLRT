# pvLRT v0.2.1

This is a new CRAN submission for R package pvLRT, which provides a suite of functions implementing likelihood ratio test-based approaches to pharmacovigilance. In this version we have made sure that:

* all examples take < 10s time to pass the CRAN check. One example failed this test during initial CRAN submission.


## Test environments

* local Windows 10 install, R 4.1.2
* Ubuntu 20.04.3 LTS, R 4.1.2
* win-builder (devel, release, oldrelease)

## R CMD Check Results

- There were no ERRORs or WARNINGs.
- There were no NOTEs in local Windows 10 and Ubuntu 20.04.3 checks
- There was one note in win-builder: 

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Saptarshi Chakraborty <chakra.saptarshi@gmail.com>'

New submission

Possibly mis-spelled words in DESCRIPTION:
  Pharmacovigilance (3:23)
  pharmacovigilance (8:97)

We believe this NOTE is a false alarm, as 'pharmacovigilance' refers to the field relating to the collection, detection, assessment, monitoring, and prevention of adverse effects with pharmaceutical products (Wikipedia; https://en.wikipedia.org/wiki/Pharmacovigilance). 
