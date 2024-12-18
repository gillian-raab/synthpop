\name{write.syn}
\alias{write.syn}
\title{Exporting synthetic data sets to external files}
\description{Exports synthetic data set(s) from synthesised data set 
(\code{synds}) object to external files of selected format. 
Currently supported file formats include:  SAS, csv, tab,
rda, RData and txt. Options to write SPSS or Stata are no longer
available. The filetype "SAS" writes two files,  a .SAS  file 
to read in to SAS, and a .txt file with the data read in by the .SAS file
In all cases information about the synthesis is written into a separate file.  
}
\usage{
write.syn(object, filename, 
filetype = c("csv", "tab", "txt",
                     "SAS", "rda", "RData"), 
 save.complete = TRUE, extended.info = TRUE, ...)
}
\arguments{
  \item{object}{an object of class \code{synds}, which stands
    for 'synthesised data set'. It is typically created by
    function \code{\link{syn}} and it includes \code{object$m} 
    synthesised data set(s).}
  \item{filename}{the name of the file (excluding extension) which the 
    synthetic data are to be written into. For multiple synthetic data sets
    it will be used as a prefix folowed respectively by \code{_1}, \code{_...},
     \code{_m}.}
  \item{filetype}{a desired format of the output files.}
  \item{save.complete}{a logical value indicating whether a complete 
    'synthesised data set' (\code{synds}) object should be saved into a file
    (\code{synobject_filename.RData}).}
  \item{extended.info}{a logical value indicating whether extended information
    should be saved into an information file.}
  \item{\dots}{additional parameters passed to write functions.}  
}
\value{File(s) with synthesised data set(s) and a text file with information 
  about synthesis are produced. Optionally a complete synthesised data set
  object is saved into \code{synobject_filename.RData} file.}



