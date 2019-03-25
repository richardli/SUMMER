## Test environments
* local OS X install, R 3.4.2
* ubuntu 12.04 (on travis-ci), R 3.4.2
* win-builder (devel and release)

0 errors | 0 warnings | 0 notes

## Version 0.1.0 Checks

Result: NOTE 
    Package suggested but not available for checking: ‘INLA’ 
    
INLA is not held in a common repository so package so cannot be made available for checking.

checking re-building of vignette outputs ... [2s/3s] WARNING
Error in re-building vignettes:
  ...
Installing package into '/home/hornik/tmp/R.check/r-release-gcc/Work/build/Packages'
(as 'lib' is unspecified)

We added the vignette output to inst/doc. We also have switched to a static vignette, but with source code available and shown.

## Version 0.2.0 Update

* We expanded authorship of the package.
* We updated the links for URL and BugReports.
* We reworked our vignette.
* We added the capability to use a yearly model. Previously, only 5-year increments were possible.
* We changed the plotting function aesthetics and adjusted for the new yearly model. 
* We renamed our included data, added another data set, and clarified that the data do not represent any real country.
* We added the capability to do a simpler model that uses only binary outcome.
* We made minor updates to the documentation and examples.
* Add vignette to inst/doc to prevent WARNINGs.
* Switch to static vignette with source code available and shown so that INLA does not need to be installed.

## Version 0.2.1 Update
* We developed a new function for simulating spatial and temporal random effects.
* We fixed a problem with fitINLA() function that may arise when user-customized dataset with missing areas.
* We updated the INLA repository as requested.
* We removed packages that are no longer needed to remove NOTEs from https://cran.rstudio.com/web/checks/check_results_SUMMER.html.

## Version 0.2.2 Update
* Remove deprecated calls to cBind() and rBind().
* Update to vignette.
* Minor edits to documentation.
* Add in new data examples.
* Improve generalizability of some functions.

## Version 0.2.3 Update
* Improve generalizability of some functions.
