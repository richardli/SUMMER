## Test environments
* local OS X install, R 3.4.2
* ubuntu 12.04 (on travis-ci), R 3.4.2
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. There were 3 NOTEs:

* Maintainer: ‘Bryan D Martin <bmartin6@uw.edu>’

  New submission

  Suggests or Enhances not in mainstream repositories:
    INLA 
  Availability using Additional_repositories specification:
    INLA   yes   https://www.math.ntnu.no/inla/R/stable
  
**Response:** The first part does not seem to be an issue. For the latter, the INLA repository is included in Additional_repositories.

* Files ‘README.md’ or ‘NEWS.md’ cannot be checked without ‘pandoc’ being installed.

**Response:** I can remove this NOTE by adding README.md file to .Rbuildignore, but I don't think I should do this.

* Possibly mis-spelled words in DESCRIPTION:
  Spatio (3:8)
  spatio (9:63)
  
**Response:** The word spatio is not misspelled. It is part of the hyphenated word spatio-temporal.

## Round 1 Response

* Thanks, please omit the redundant 'in R' in your title.

**Response:** Done.

* Please add a reference for this method in the 'Description' field of your DESCRIPTION file in the form
  authors (year) <doi:...>
  authors (year) <arXiv:...>
  authors (year, ISBN:...)
  with no space after 'doi:', 'arXiv:' and angle brackets for auto-linking.
  
**Response:** The citation has been added. This added 1 NOTE:
checking DESCRIPTION meta-information ... NOTE
Malformed Description field: should contain one or more complete sentences.

This is due to the "." from "Mercer et al.". All sentences are complete.
  
* Please add more small executable examples in your Rd-files.

**Response:** All exported functions now include examples.

## Round 2 Response

* Thanks, but descriptions longer than one line should start with a space or tab. Otherwise a new description line is identified as a new DESCRIPTION field.

**Response:** Thanks for pointing this out. It has been fixed.

* Please also omit the blank in <doi: 10.1214/15-AOAS872>.

**Response:** We have searched for this blank in all files by checking every occurrence of "doi", and it does not occur. Just in case, we have manually re-typed  "doi:1" in each occurence of the reference.

## Version 0.2.0 Update

* We expanded authorship of the package.
* We updated the links for URL and BugReports.
* We reworked our vignette.
* We added the capability to use a yearly model. Previously, only 5-year increments were possible.
* We changed the plotting function aesthetics and adjusted for the new yearly model. 
* We renamed our included data, added another data set, and clarified that the data do not represent any real country.
* We added the capability to do a simpler model that uses only binary outcome.
* We made minor updates to the documentation and examples.

