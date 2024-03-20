synorig.compare <- function(syn,orig, print.flag = TRUE){
 needsfix = FALSE
 unchanged <- TRUE
  ## to convert any tibbles or matrices
  if (!is.data.frame(syn)) { syn <- data.frame(syn) ; unchanged = FALSE}
  if (!is.data.frame(orig)) { orig <- data.frame(orig); unchanged = FALSE}
  
  names_syn <- names(syn)
  names_orig <- names(orig)
  
  if (any(!(names_syn %in% names_orig))) {
  cat("Variables",names_syn[!(names_syn%in% names_orig)],"in synthetic but not in original: comparisons impossible
  ****************************************************************\n")
    needsfix <- TRUE
  }

	i.syn<- (1:length(names_syn))[names_syn %in% names_orig]
	i.orig<- (1:length(names_orig))[names_orig %in% names_syn]

  len.common <- length(i.syn)
  if (print.flag) cat(len.common, "variables in common out of" , length(names_syn), "in syn and out of", length(names_orig),"in orig\n")
  
 ##-------------- change common character variables to factors---------------------------
  nch_syn <- 0 ; nch_orig <- 0
  for( i in i.syn ) {
    if ( class( syn[,i] )[1] == "character") {
      syn[,i] <- factor( syn[,i] ) 
      nch_syn <- nch_syn +1
      unchanged = FALSE
    }
  }
  for( i in i.orig ) {
    if ( class( orig[,i] )[1] == "character") {
      orig[,i] <- factor( orig[,i] ) 
      nch_orig <- nch_orig +1
      unchanged = FALSE
    }
  }
  if (print.flag && nch_syn + nch_orig > 0) cat("\nCommon variables changed from character to factor",nch_syn,"in syn",nch_orig,"in orig\n\n")

  ##--------------- check data types in common variables----------------------------------

  for (i in 1:length(i.syn)){
    i1 <- i.syn[i]
    i2 <- i.orig[i]
    if ( length(class(syn[,i1])) !=  length(class(orig[,i2])) ||
         !all(class(syn[,i1])  ==  class(orig[,i2]) )	) {
      cat("\nDifferent classes for",names_syn[i1],"in syn:",class(syn[,i1]),"in orig:",class(orig[,i2]),
         "\n*************** must be changed before utility or disclosure functions *********\n")
      needsfix <- TRUE
    }
  
  } 
  ##--------------------------- compare missingness and levels for factors------------------------------
  for (i in 1:length(i.syn)){
    i1 <- i.syn[i]
    i2 <- i.orig[i]
    if (!any(is.na(orig[,i2])) & any(is.na(syn[,i1])) ) { 
    cat("Missing data for common variable", names(syn)[i1], "in syn but not in orig warning:
 ***********  comparison methods will fail need to sort out problem*************\n")
      needsfix = TRUE
    }
    if ( is.factor(syn[,i1])  & is.factor(orig[,i2]) ) {  
    if ( any(is.na(orig[,i2])) ) {
        syn[,i1] <- addNA(syn[,i1])
        orig[,i2] <- addNA(orig[,i2])
    }
      syn[,i1] <- factor(as.character(syn[,i1])) ### to get in right order
      orig[,i2] <- factor(as.character(orig[,i2])) ### maybe not needed
      lev1 <- levels(syn[,i1]) ;lev2 <- levels(orig[,i2])
      
      #cat("\n",(lev1),"\n",(lev2),"\n\n", lev1==lev2)
      if ( length(lev1) != length(lev2)  ||
           !all(lev1[!is.na(lev1)] == lev2[!is.na(lev2)]) ) {
        
        cat("Different levels of",names(orig)[i2],"in 1st and 2nd data set; now combined\n",
         "syn:"  , lev1,"\n","orig:",lev2,"\n\n not just due to absence of missings in syn, levels combined BUT please check  if looks OK",
  "\n************************************* checking highly recommended*******************\n\n")
        syn[,i1] <- factor(as.character(syn[,i1]), exclude = NULL, levels = unique(c(lev1,lev2)))
        orig[,i2] <- factor(as.character(orig[,i2]), exclude = NULL, levels = unique(c(lev1,lev2)))
        unchanged = FALSE
      }
    }
  }
  if (print.flag & needsfix) cat("More fixing needed before you can compare and evaluate these data sets.")
  res <- list(syn=syn, orig = orig, needsfix = needsfix, unchanged = unchanged)
}
