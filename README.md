
# <img src="inst/logo.png" width="50" height="50"> synthpop with disclosure functions

### for more information about _synthpop_ you can visit **[www.synthpop.org.uk](www.synthpop.org.uk)**. Look under resources for materials from our courses on how to use the package.

### This is a new version (v1.9-0) of the synthpop package.

### The version on CRAN is 1.8-0 as of 6/2/2024.

### You can find it here **[_synthpop_ on CRAN](https://CRAN.R-project.org/package=synthpop)**.

The biggest change in the new version is the inclusion of two new functions to assess the disclosure risk of synthetic data, disclosure() and disclosure.summary().

Each of them can take as input either an oject of class synds that includes synthetic data set(s) created by synthpop or a synthetic data set or a list of synthetic data sets created by any method from the same original data.

These are described in the vignette disclosure.pdf.

To install the package from this repo in Rstudio use the command

### devtools::install_github("Gillian-Raab/synthpop",build_vignettes = TRUE) ###

If you then go to the package index you will see the documentation file 
disclosure.pdf in the vignettes directory,

## WARNING: 
If your computer does not have a program pdflatex to convert latex files to pdfs the installation is likely to fail. If so you will need to install it without vignettes using the command.

### devtools::install_github("Gillian-Raab/synthpop") ###

You will not find the new vignette (disclosure.pdf), that explaining the new disclosure functions, but you can find it in the vignettes directory on github. The other vignettes are unchanged from the version of synthpop on CRAN.
 
## Other changes include :
 
 Adding Gaussian DP noise in the DP versions of syn.catall and syn.ipf
 
 Extending syn.ipf to allow all 3-way or 4-way interactions as a group option
 
 New function (mergelevels.syn) for grouping small levels of factors
 
 New function (synorig.compare) that checks that the structure of the original and synthetic data
  are compatible. Called by several functions that evaluate the synthetic data when the synthetic 
  data are supplied as data set(s) that are not objects of class synds.
 



## Acknowledgements

Many thanks to Jason Klebes for all his help with github and his patience with my ignorance.
