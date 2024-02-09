\name{disclosure}
\alias{disclosure}
\alias{disclosure.synds}
\alias{disclosure.data.frame}
\alias{disclosure.list}
\alias{print.disclosure}
\title{Disclosure measures}
\description{
  Calculates disclosure measures for synthetic data. 
}
\usage{
\method{disclosure}{synds}(object, data, keys, target, denom_lim = 5,
  exclude_ov_denom_lim = FALSE, print.flag = TRUE,
  digits = 2, usetargetNA = TRUE, usekeysNA = TRUE,
  exclude.keys =NULL, exclude.keylevs = NULL, exclude.targetlevs = NULL,
  to.print =c("short"), \dots )

\method{disclosure}{data.frame}(object, data, keys, target, denom_lim = 5, 
  exclude_ov_denom_lim = FALSE, print.flag = TRUE,
  digits = 2, usetargetNA = TRUE, usekeysNA = TRUE,
  exclude.keys =NULL, exclude.keylevs = NULL, exclude.targetlevs = NULL,
  to.print =c("short"), \dots )

\method{disclosure}{list}(object, data, keys, target, denom_lim = 5, 
  exclude_ov_denom_lim = FALSE, print.flag = TRUE,
  digits = 2, usetargetNA = TRUE, usekeysNA = TRUE,
  exclude.keys =NULL, exclude.keylevs = NULL, exclude.targetlevs = NULL,
  to.print =c("short"), \dots )

\method{print}{disclosure}(x,  to.print =NULL, digits = NULL, \dots)}

\arguments{
  \item{object}{an object of class \code{synds}, which stands for 'synthesised 
    data set'. It is typically created by function \code{syn()} and it includes 
    \code{object$m} synthesised data set(s) as \code{object$syn}. This a single 
    data set when \code{object$m = 1} or a list of length \code{object$m} when 
    \code{object$m > 1}. Alternatively, when data are synthesised not using 
    \code{syn()}, it can be a data frame with a synthetic data set or a list 
    of data frames with synthetic data sets, all created from the same original 
    data with the same variables.}  
  \item{data}{the original (observed) data set.}
  \item{keys}{vector of variable names in data and synthetic data to act as 
       quasi-identifiers for identity or attribute disclosure.}
  \item{target}{name of target variable for attribute disclosure.}
  \item{denom_lim}{Size of the denominator from an apparent attribute
    disclosure risk in the synthetic data, that will trigger the 
    printing of information on the key and variable levels with 
    the maximum  denominators.}
  \item{exclude_ov_denom_lim}{logical to exclude disclosive records
    with denominators greater than \code{denom_lim} from DiSCO.}
  \item{print.flag}{logical value as to whether a line is printed as
    disclosure is calculated for each symthetic data set.}
  \item{digits}{number of digits to print for disclosure measures.} 
  \item{usetargetNA}{determines whether NA values in target are to be used 
            in checking for disclosure}
  \item{usekeysNA}{determines whether NA values in keys are to be used 
            in checking for disclosure.}
  \item{exclude.keys}{vector of names of keys that, with the next two items
     will define the target and key options to be excluded from the calculation 
         of disclosure measures.}
  \item{exclude.keylevs}{vector of the same length as exclude.keys that give 
      the levels  to be excluded for the corresponding key.}
  \item{exclude.targetlevs}{vector of the same length as exclude.keys that give 
      the levels  of target that will be excluded for each key and key level.}
  \item{to.print}{Vector to determine what aspects of an object of
                  class disclosure will be printed. Must consist  of 
                  one or more of the following "short", "ident",
                  "attrib","allCAPs", "denom_details","all".
                  Default is "short" giving a brief summary.}
  \item{\dots}{additional parameters}  
  \item{x}{an object of class \code{disclosure}.}
}

\details{ Calculates identity disclosure measures for a for a set of keys,
         (quasi identifiers) and attribute disclosure measures for  one
         variable from the same  set of keys considered as a target. The 
         function \code{\link{disclosure.summary}} calls this function and 
         summarises the attribute disclosure measures for a set of targets.
         See the vignette 
}

\value{An object of class \code{disclosure} which is a list with the following 
  components. 
  \item{call}{the call that created the object.}
  \item{ident}{Table of measures of identity disclosure
              one for each synthesis. Measures are "UiO","UiS","UiSiO" and
              "repU". See vignette disclosure.pdf for an explanation of these 
              and the following measures.}
  \item{attrib}{Table of measures of attribute disclosure  one for each 
  synthesis. THese include "DiO","DiS","DiSiO","DisCo" and "DiSCO_DiO. 
  The table also as gives the maximum and mean of the denominators for the
  "DiSCO" measure i.e. the distribution for every record that leads to a 
  correct disclosure of the number of observations with
  the same keys and the same correct target. Large denominators are often an
  indication that the discloure is something that might be expected from
  prior knowledge of relations.} 
  \item{allCAPs}{Table of the following measures of correct attribution
  probability: "baseCAPd","CAPd", "CAPs" , "DCAP" and  "TCAP"'} 
  \item{denom_details }{A list of length number of syntheses giving details 
  for each of the two-way combinations of target and keys where the 
  denominator exceeds \code{denom_lim}. If no denominators exceed
  denom_lim for one of the syntheses then the lisat element is NULL. 
  SUch disclosive combinations are often associated with a high prior 
  probability of   the target from just this level of one of the  keys 
  in the original data. The denominators  and proportions for this combination 
  in the original data are also given.
  Note these combinations will be identified  for data without 
  excluded or missing values of keys if there are any excluded combinations  
  or target if any of \code{usekeysNA} or \code{usetargetNA} are FALSE.} 
  \item{Nidentexclude}{Number of records excluded due to missing values in keys}
  \item{list_excludes}{A list of length number of syntheses with number of 
  records excluded from attribute measures for different reasons.}
  \item{keys}{as input}
  \item{target}{as input}
  \item{digits}{as input}
  \item{denom_lim}{as input}
  \item{usetargetNA}{as input}
  \item{usekeysNA}{as input}
  \item{Nout_miss}{Number of records excluded when any of \code{usetargetNA} or \code{usekeysNA} are FALSE.}
  \item{exclude.keys}{as input}
  \item{exclude.keylevs}{as input}
  \item{exclude.targetlevs}{as input}
  \item{Norig}{number of records in original data}
  \item{Nsyn}{Number of records in each synthetic data set.}
  \item{Nout_exclude}{Number of records excluded as identified by \code{exclude.keys}, \code{exclude.keylevs} and \code{exclude.targetlevs}.}
  \item{to.print}{as input}
  \item{call}{text of call to the disclosure function that created the object,}

}
\note{See package vignette disclosure.pdf for additional information including
  formal definitions of all quantities and worked examplesexamples.}

\references{
 See references in package vignette}
 
\seealso{  \code{\link{syn}} \code{\link{disclosure.summary}}
}

\examples{
library(synthpop)

ods <- SD2011[, c("sex", "age", "edu", "marital", "region", "income","ls",
              "trustfam","ymarr","smoke")]
odsF <- numtocat.syn(ods, numtocat = "income", catgroups = 7, cont.na = list(income = -8))
s1 <- syn(odsF$data, method = "ctree",seed = 75, m=5)
disc1 <- disclosure(s1, odsF$data, target = "income", 
         keys = c("sex", "age", "edu","marital", "region"))
print(disc1)
print(disc1, to.print = c("ident", "attrib","allCAPs"), digits = 3)

disc2 <- disclosure(s1$syn[[3]], odsF$data, target = "income", keys = c("sex", "age", "edu",
                                "marital", "region"), usetargetNA = FALSE)
print(disc2, to.print = "attrib")

disc3 <- disclosure(s1$syn[[2]], odsF$data, target = "smoke", keys = c("sex", "age", "edu",
         "marital", "region"), denom_lim =5)
print(disc3, to.print = "attrib")
print(disc3, to.print = "denom_details")
}
