###----------------------- disclosure.summary------------------------------
disclosure.summary <- function(object, data, ...) UseMethod("disclosure.summary")

###-----disclosure.summary.default-----------------------------
disclosure.summary.default <- function(object, ...)
  stop("No disclosure.summary method associated with class ", class(object), call. = FALSE)

###-----disclosure.summary.data.frame---disclosure.list--------
disclosure.summary.data.frame <- disclosure.summary.list <- 
  function(object, data, keys , targets = NULL, denom_lim = 5, 
           exclude_ov_denom_lim = FALSE, print.flag = TRUE,  
           usetargetsNA = TRUE,  usekeysNA = TRUE, 
           ident.meas = "repU", attrib.meas = "DiSCO",
           digits = 2, plot = TRUE, print = TRUE,  ...)
    
  {
    if (is.null(object)) stop("Requires parameter 'object' to give name of the synthetic data.\n", call. = FALSE)   
    
    if (is.list(object) & !is.data.frame(object)) m <- length(object)
    else if (is.data.frame(object)) m <- 1
    else stop("object must be a data frame or a list of data frames.\n", call. = FALSE)
    
    object <- list(syn = object, m = m) 
    class(object) <- "synds"
    
    
    res <- disclosure.summary.synds(object, data, keys = keys , targets = targets, 
           denom_lim = denom_lim,  exclude_ov_denom_lim = exclude_ov_denom_lim, 
           print.flag = print.flag ,  usetargetsNA = usetargetsNA, usekeysNA = usekeysNA, 
           ident.meas = ident.meas, attrib.meas = attrib.meas, 
           digits = digits, plot = plot, print = print, ...) 
    res$call <- match.call()
    return(res)
  }


###-----disclosure.summary.synds-------------------------
disclosure.summary.synds <-     function(object, data, keys , targets = NULL, denom_lim = 5, 
                                         exclude_ov_denom_lim = FALSE, print.flag = TRUE,  
                                         usetargetsNA = TRUE,  usekeysNA = TRUE, 
                                         ident.meas = "repU", attrib.meas = "DiSCO",
                                         digits = 2, plot = TRUE, print = TRUE,  ...)
{
  ###----------------check input parameters ----
  if (!(is.data.frame(data)) )   stop("data  must be a data frame \n\n", call. = FALSE)
  data <- data.frame(data) ## in case it is table or a tibble
  if (!( inherits(object,"synds")))   stop(" object must be an object of class synds\n\n", call. = FALSE)
  
  if (object$m ==1) {
    names.syn <- names(object$syn)
  }
  else names.syn <- names(object$syn[[1]])

  # targets must be variables in  data and object$syn
  # keys must be a vector of variable names in data and in object$syn
  # target must not be in keys
  if (is.null(targets) ) targets <- names(data)[!names(data) %in% keys]

  if (!(all(keys %in% names(data)) && all(keys %in% names.syn) 
        && all(targets %in% names(data)) & all(targets %in% names.syn )))
    stop("keys and targets must be variables in data and synthetic data. \n", call. = FALSE)
  if ( any(targets %in% keys) )
    stop("no targets can be in keys \n", call. = FALSE)
  if (any(names(data)  == "target"))
    stop("your data has a variables called 'target' please rename in original and synthetic data.\n\n", call. = FALSE) 
  
  if (is.null(usekeysNA) ) usekeysNA <- TRUE
  
  # get keys in same order as in data
  keys <- names(data)[names(data) %in% keys]
  if (!(length(ident.meas) == 1  && length(attrib.meas) == 1 &&
        ident.meas %in% c("repU", "UiSiO") && attrib.meas %in% c("DiSCO", "DiSCOU","DCAP")) )
        stop('ident.meas and attrib.meas must be single values
             from c("repU", "UiSiO") and  c("DiSCO", "DiSCO_DiO", "DCAP","TCAP") respectively\n\n', call. = FALSE)

  if (length(usekeysNA) == 1) usekeysNA <- rep(usekeysNA, length(keys))
  if (length(usetargetsNA) == 1) usetargetsNA <- rep(usetargetsNA, length(targets))
  
  if (length(usekeysNA) != length(keys)) stop("usekeysNA must be same length as keys", call. = FALSE)
  if (length(usetargetsNA) != length(targets)) stop("usetargetsNA must be same length as targets", call. = FALSE)
  

###-------------------------- create output matrix--------------------
  result <- matrix(NA, length(targets),3)
  result[,3] <- FALSE
  dimnames(result)<- list(targets, c("DiO","","check denom"))
    dimnames(result)[[2]] [2] <- attrib.meas
  if (attrib.meas == "DCAP") dimnames(result)[[2]] [1] <- "CAPd"

  for (i in 1: length(targets)) {
    if (print.flag) cat("------------------",i,targets[i],"-------------------","\n")

    ttt <-disclosure(object, data, target = targets[i], keys = keys, 
            denom_lim = denom_lim, exclude_ov_denom_lim = exclude_ov_denom_lim,
            print.flag = print.flag, digits =digits,
            usekeysNA = usekeysNA, usetargetNA = usetargetsNA[i])
            
    if (ident.meas == "repU") ident.syn = mean(ttt$ident$repU[1])
    if (ident.meas == "UiSiO") ident.syn = mean(ttt$ident$UiSiO)
    ident.orig <- mean(ttt$ident$UiO)

    if (attrib.meas == "DiSCO" )
      result[i,1:2] <- c(mean(ttt$attrib$DiO), mean(ttt$attrib$DiSCO ))
    if (attrib.meas == "DiSCOU" )
      result[i,1:2] <- c(mean(ttt$attrib$DiO), mean(ttt$attrib$DiSCOU ))
    if (attrib.meas == "DCAP" )
      result[i,1:2] <- c(mean(ttt$allCAPs$CAPd), mean(ttt$allCAPs$DCAP ))
    if (length(ttt$denom_details) > 0) result[i,3] <- TRUE

  }
    
    ###-----ntoc---------------------------------------------------------------
    # to make labels for variables of constant length from an integer
    ntoc <- function(x)
    {
      nch <- nchar(max(x))
      res <- as.character(x)
      res <- str_pad(res, nch, "left", "0")
      return(res)
    }
    result <- data.frame(result)
    result <- result[order(result[,1]),]

    dimnames(result)[[1]] <- paste(ntoc(1:dim(result)[1]),dimnames(result)[[1]] )
    if (any(result[,3] >0)) dimnames(result)[[1]][result[,3] > 0] <- 
                      paste(dimnames(result)[[1]]," - check")[result[,3] >0] 

    nnn <- dimnames(result)[[1]]
    toplot <- result[rep( 1:dim(result)[[1]], rep(2,dim(result)[[1]]) ),]
    toplot$name <- rep(nnn, rep(2,length(nnn)))

    for ( i in 1:(dim(result)[1]) )       toplot[i*2,1] <- toplot[i*2,2]

    toplot<- as.data.frame(toplot[,-2])
    toplot$measure <- dimnames(result)[[2]][1]
    for ( i in 1:(dim(result)[1]) )  toplot$measure[i*2] <- dimnames(result)[[2]][2]

    names(toplot)[1] <- "VALUE"
    toplot <- toplot[,c(1,3,4)]

     attrib.plot <- ggplot(toplot) + 
     geom_point(data = toplot, size=5, aes(colour = .data$measure, x=.data$VALUE, y=.data$name)) +
       xlim(0,100)  +
       labs(x = "Disclosure measure", y = "", 
       title = "Comparison of attribute disclosure measures",
       subtitle= paste(dimnames(result)[[2]][2],"for synthetic data  to", dimnames(result)[[2]][1] ,"for original data.")) +
      geom_line(data = toplot, mapping = aes(x=.data$VALUE, y=.data$name), arrow = arrow(length=unit(0.30,"cm"), ends="first", type = "closed"))
     

     res <- list( attrib.table = result, attrib.plot = attrib.plot, keys = keys,
                       ident.orig = ident.orig, ident.syn = ident.syn,
                       denom_lim = denom_lim, exclude_ov_denom_lim = exclude_ov_denom_lim,
                       digits = digits, usetargetsNA = usetargetsNA, usekeysNA = usekeysNA, 
                       ident.meas = ident.meas, attrib.meas = attrib.meas, m = object$m,
                       plot = TRUE, print = TRUE)
     class(res) <- "disclosure.summary"
     return(res)
  }        


###-----print.disclosure.summary-----------------------------------------------
print.disclosure.summary <- function(x,  digits = NULL,  
                                 plot = NULL, print = NULL, ...) {
  
   if (is.null(digits)) digits <- x$digits
  if (is.null(plot)) plot <- x$plot
  if (is.null(print)) print <- x$print
  
  
  if (!print) {
    cat("\nFor more details use print = TRUE.\n")
  } else {
    if (x$m >1) cat("\nResults are averaged over" ,x$m,"syntheses")
    cat("\nIdentity disclosure measures\n")
    cat("from keys:", x$keys,"\n")
    cat("\nFor original, percent data unique from keys  ( UiO )", x$ident.orig,"\n")
    cat("For synthetic, replicated uniques from keys (", x$ident.meas,")", x$ident.syn,"\n")
    
    cat("\n\nTable of selected disclosure measures, same keys,\n")
    cat("Variables Ordered by original disclosure\n\n")
    print(x$attrib.table)
  }

  if (plot) {
    print(x$attrib.plot)
  }
  if (any(x$attrib.table$check >0)) cat("\n\nCheck denom_details in output from disclosure()\n",
                              "for targets with check >0\n")
  invisible(x)
}




