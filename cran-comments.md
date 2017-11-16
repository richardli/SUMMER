## Test environments
* local OS X install, R 3.4.2
* ubuntu 12.04 (on travis-ci), R 3.4.2
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. There was 1 NOTE:

* Possibly mis-spelled words in DESCRIPTION:
  Spatio (3:8)
  spatio (9:63)

* Suggests or Enhances not in mainstream repositories:
  INLA
  
The word spatio is not misspelled. It is part of the hyphenated word spatio-temporal.
  
The INLA repository is listed under Additional_repositories.
