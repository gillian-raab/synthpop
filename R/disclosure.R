###-----new version with big denom checks---------------------
disclosure <- function(object, data, ...) UseMethod("disclosure")

###-----disclosure.Qdefault-----------------------------
disclosure.default <- function(object, ...)
  stop("No disclosure method associated with class ", class(object), call. = FALSE)

###-----disclosure.data.frame---disclosure.list--------
disclosure.data.frame <- disclosure.list <- 
  function(object, data, keys , target , denom_lim = 5,
           exclude_ov_denom_lim = FALSE, print.flag = TRUE, digits = 2, 
           usetargetNA = TRUE, usekeysNA = TRUE, 
           exclude.keys =NULL, exclude.keylevs = NULL, 
           exclude.targetlevs = NULL, to.print =c("short"), ...) 
    {
    if (is.null(object)) stop("Requires parameter 'object' to give name of the synthetic data.\n", call. = FALSE)   
    
    if (is.list(object) & !is.data.frame(object)) m <- length(object)
    else if (is.data.frame(object)) m <- 1
    else stop("object must be a data frame or a list of data frames.\n", call. = FALSE)
    
    object <- list(syn = object, m = m) 
    class(object) <- "synds"
    
    
    res <- disclosure(object, data, keys , target , denom_lim = denom_lim, 
                      exclude_ov_denom_lim = exclude_ov_denom_lim, 
                      print.flag = print.flag, digits = digits, 
                      usetargetNA = usetargetNA, usekeysNA = usekeysNA, 
                      exclude.keys = exclude.keys, exclude.keylevs =  exclude.keylevs,
                      exclude.targetlevs =  exclude.targetlevs,
                      to.print = to.print, ...) 
    
    res$call <- match.call()
    return(res)
  }
    
###-----disclosure.synds-------------------------
disclosure.synds <-  function(object, data, keys , target , denom_lim = 5,
           exclude_ov_denom_lim = FALSE, print.flag = TRUE, digits = 2, 
           usetargetNA = TRUE, usekeysNA = TRUE, 
           exclude.keys =NULL, exclude.keylevs = NULL, 
           exclude.targetlevs = NULL, to.print =c("short"), ...) 
  
  {

###----------------check input parameters ----
    if (!(is.data.frame(data)) )   stop("data  must be a data frame \n\n", call. = FALSE)
   data <- data.frame(data) ## in case it is table or a tibble
    if (!( inherits(object,"synds")))   stop(" object must be an object of class synds\n\n", call. = FALSE)

 if (object$m ==1) {
   names.syn <- names(object$syn)
 }
   else names.syn <- names(object$syn[[1]])

  # target must be a single variable in  data and object$syn
  # keys must be a vector of variable names in data and in s
  # target must not be in keys

if (!(all(keys %in% names(data) & keys %in% names.syn) & all(target %in% names(data)) & target %in% names.syn))
    stop("keys and target must be variables in data and synthetic data. \n", call. = FALSE)
   if (any(duplicated(keys)))
     stop("keys cannot include duplicated values. \n", call. = FALSE)
  if (!(length(target)==1))
    stop("target must be a single variable \n", call. = FALSE)
   if (target %in% keys)
     stop("target cannot be in keys \n", call. = FALSE)
   if (any(names(data)  == "target"))
     stop("your data have a variables called 'target' please rename in original and synthetic data.\n\n", call. = FALSE) 
   
   if (!(length(usetargetNA)==1))
     stop("usetargetNA must be a single logical value \n", call. = FALSE)
   if (length(usekeysNA) ==1 ) usekeysNA <- rep(usekeysNA, length(keys))
   if (!(length(usekeysNA)==length(keys)))
     stop("usekeysNA must be a logical value of same length as keys\n", call. = FALSE)
   
  # get keys and usekeysNA in same order as in data
   #cat(keys,"keys\n", usekeysNA,"usekeysNA\n")
   #cat((1:length(names(data)))[names(data) %in% keys],"(1:length(names(data)))[names(data) %in% keys]\n")
   oldkeys <- keys
   keys <- names(data)[(1:length(names(data)))[names(data) %in% keys]]
   usekeysNA <- usekeysNA[match(oldkeys,keys)]
   #cat(keys,"keys\n", usekeysNA,"line 79 usekeysNA\n")

     # check excluded combinations
  if (!is.null(exclude.keys)) {
    if (! length(exclude.keylevs) == length(exclude.keys) & 
        length(exclude.keylevs) == length(exclude.targetlevs) ) stop("All excludes must be the same length\n" , call. = FALSE)
    if (!all(exclude.keys %in% keys)) stop("exclude.keys must be the name of one of your keys", call.= FALSE)
  }
   
   if(!is.null(denom_lim)){
     if (!(round(denom_lim) == denom_lim && denom_lim > 0 )){
       cat(denom_lim ,"denom_lim\n")
       stop("\ndenom_lim must be an integer >0\n", call.= FALSE)
     }
   }
###---------------------- define output items-------------
  allCAPs         <- matrix(NA,object$m,5)
  attrib   <- matrix(NA,object$m,7)
  ident <-     Nidentexclude <- matrix(NA,object$m ,4)

  dimnames(allCAPs) <- list(1:object$m, c("baseCAPd","CAPd", "CAPs", "DCAP","TCAP"))
  dimnames(attrib) <- list(1:object$m,c("DiO","DiS","DiSiO","DiSCO", "DiSCO_DiO", "max_denom","mean_denom"))
  dimnames(ident) <- dimnames(Nidentexclude) <- list(1:object$m,c("UiO", "UiS","UiSiO", "repU"))
  list_denom <-list(1:object$m)
  list_excludes <- list(1:object$m)

  
  dd <- data ## rename target variable to target
  targetx =dd[, names(dd) == target]
  dd$target <- targetx
 # cat(names(dd),"names(dd)\n")
  dd <- as.data.frame(dd[,names(dd) %in% c( keys, "target")])


  if ( !is.factor(dd$target) )   dd$target <- factor(dd$target)
  for (i in 1:length(keys))   dd[,names(dd) == keys[i]] <- factor( dd[,names(dd) == keys[i]])
###------------------- loop over object$m syntheses---------------
 
  for ( jj in 1:object$m) {
###-------------------------------- PRINT WHERE AT ----------------------------------   
    if (print.flag) cat("Synthesis",jj,"\n") 
    
   if (object$m > 1 )    ss <- object$syn[[jj]]
    else    ss <- object$syn
    
###--------------- make new variable in synthetic data from the target variable called target--------------------------
    targetx =ss[, names(ss) == target]
    ss$target <- targetx
    ss <- ss[,names(ss) %in% c("keys","target", keys)]
    
    if ( !is.factor(ss$target) )   ss$target <- factor(ss$target)
    for (i in 1:length(keys))   ss[,names(ss) == keys[i]] <- factor( ss[,names(ss) == keys[i]])

###----------------------- form tables from keys and target-----------------     
###  replace Missing values with factor value of "Missing" to make tables easier
    tomissing <- function(x){
      if (!is.factor(x)) stop("x must be a factor\n",.call =FALSE)
      x <- as.character(x)
      x[is.na(x)] <- "Missing"
      x <- factor(x)
    }

    if ( any(is.na(dd$target)) ) dd$target <- tomissing(dd$target)
    if ( any(is.na(ss$target)) ) ss$target <- tomissing(ss$target)

    for ( i in 1:length(keys)) {
      if ( any(is.na(dd[,names(dd) == keys[i]])) ) dd[,names(dd) == keys[i]] <-
                                      tomissing(dd[,names(dd) == keys[i]])
      if ( any(is.na(ss[,names(ss) == keys[i]])) ) ss[,names(ss) == keys[i]] <-
          tomissing(ss[,names(ss) == keys[i]])

 }
    
     Nd <- dim(dd)[1]
     Ns <- dim(ss)[1]
  ###------------------------ make composite variable for keys 
if (length(keys) >1) {
   ss$keys <- apply(ss[, names(ss) %in% keys],1,function(x) paste(x, collapse = " | "))
   dd$keys <-   apply(dd[, names(dd)  %in% keys],1,function(x) paste(x, collapse = " | "))
} else {
  ss$keys <- ss[, names(ss) == keys]
  dd$keys <- dd[, names(dd) == keys]
}


  
###--------------------- make tables ---------------------------

  NKd <- length(table(dd$keys))
  NKs <- length(table(ss$keys))

  tab_kts <- table(ss$target,ss$keys)
  tab_ktd <- table(dd$target,dd$keys)   ## two way target and keys table orig
  
  ###  get keys in s and d
  #
  Kd <- names(table(dd$keys))
  Ks <- names(table(ss$keys))
  Kboth <- Kd[Kd %in% Ks]
  Kall <- c(Kd,Ks)
  Kall <- Kall[!duplicated(Kall)]

  ### same thing for target
  
  Td <- names(table(dd$target))
  Ts <- names(table(ss$target))
  Tboth <- Td[Td %in% Ts]
  Tall <- c(Td,Ts)
  Tall <- Tall[!duplicated(Tall)]

  ### augment keys tables to match

   if (!(all(Kd %in% Ks))) { ## some original not found in synthetic
    extraKd <- Kd[!(Kd %in% Ks) ]
    extra_tab <- matrix(0,dim(tab_kts)[1],length(extraKd))
    dimnames(extra_tab) <- list(dimnames(tab_kts)[[1]],extraKd)
    tab_kts <- cbind(tab_kts,extra_tab)
    tab_kts  <- tab_kts[, order(dimnames(tab_kts)[[2]])] 
  }

  if (!(all(Ks %in% Kd))) {  ## extra synthetic keys not in original
    extraKs <- Ks[!(Ks %in% Kd) ]
 
    extra_tab <- matrix(0,dim(tab_ktd)[1],length(extraKs))
    dimnames(extra_tab) <- list(dimnames(tab_ktd)[[1]],extraKs)
    tab_ktd <- cbind(tab_ktd,extra_tab)
    tab_ktd <- tab_ktd[,order(dimnames(tab_ktd)[[2]])]
  }

# same thing for target

  if (!(all(Td %in% Ts))) { ## some original target levels not found in synthetic
    extraTd <- Td[!(Td %in% Ts) ]
    extra_tab <- matrix(0,length(extraTd),dim(tab_kts)[2])
    dimnames(extra_tab) <- list(extraTd , dimnames(tab_kts)[[2]])
    tab_kts <- rbind(tab_kts,extra_tab)
    tab_kts  <- tab_kts[order(dimnames(tab_kts)[[1]]),] 
  }   else extraTd <- NULL
 
  if (!(all(Ts %in% Td))) {  ## extra synthetic target levels not in original  ############### edit this
    extraTs <- Ts[!(Ts %in% Td) ]
    extra_tab <- matrix(0,length(extraTs),dim(tab_ktd)[2],)
    dimnames(extra_tab) <- list(extraTs,dimnames(tab_ktd)[[2]])
    tab_ktd <- rbind(tab_ktd,extra_tab)
    tab_ktd <- tab_ktd[order(dimnames(tab_ktd)[[1]]),] 
  }   else extraTs <- NULL
 
###------------------------- calculate proportions and margins ---------#
  
  tab_ktd_p <- sweep(tab_ktd,2,apply(tab_ktd,2,sum),"/")
  tab_ktd_p[is.na(tab_ktd_p)] <- 0
  tab_kts_p <- sweep(tab_kts,2,apply(tab_kts,2,sum),"/")
  tab_kts_p[is.na(tab_kts_p)] <- 0

  tab_kd <- apply(tab_ktd,2,sum)
  tab_td<- apply(tab_ktd,1,sum)
  tab_ks <- apply(tab_kts,2,sum)
  tab_ts <- apply(tab_kts,1,sum)

  NKall <- length(tab_kd)
  NKboth <- length(Kboth)
  NTd <- length(Td)
  NTs <- length(Ts)
  Nboth <-  sum(tab_kd[names(tab_kd) %in% Kboth])
  Nd_ins <- sum(tab_kd[names(tab_kd) %in% Ks])
###----------------- get  CAP and DCAP measures---------------------

baseCAPd <- sum((tab_td/Nd)**2)*100 ## just from univariates
CAPd <-   sum(apply(tab_ktd_p^2,2,sum)*tab_kd)/Nd*100 ## from tables
CAPs <-   sum(apply(tab_kts_p^2,2,sum)*tab_ks)/Ns*100
DCAP <-   sum(apply(tab_kts_p*tab_ktd,2,sum))/Nd*100
##restrict to uniques in d
tab_kts_pU <- tab_kts_p
tab_kts_pU[,tab_kd != 1] <- 0
TCAP <-   sum(apply(tab_kts_pU*tab_ktd,2,sum))/Nd*100

allCAPs[jj,] <- c( baseCAPd,  CAPd,  CAPs, DCAP,TCAP)

###-----------------------get  identity disclosure measures -------------
tab_ks1 <- tab_ks[tab_ks == 1] 
tab_kd1 <- tab_kd[tab_kd == 1]
tab_ks1_d <- tab_ks1[names(tab_ks1) %in% Kd]
tab_ksd1 <-tab_kd[tab_ks == 1 & tab_kd == 1]

## exclude bits with missing levels of keys
Ndropident <- rep(0,4)

if (any(!usekeysNA)) {
  dropident <- rep(0,4)
  kk <- (1:(length(keys)))[!usekeysNA]
  
  drop_ks1 <- rep(FALSE,length(tab_ks1))
  drop_kd1 <- rep(FALSE,length(tab_kd1))
  drop_ks1_d <- rep(FALSE,length(tab_ks1_d))
  drop_ksd1 <- rep(FALSE,length(tab_ksd1))
  
  for ( i in kk) { 

    drop_ks1   <- drop_ks1   | word(names(tab_ks1),i, sep = fixed(" | ")) == "Missing"
    drop_kd1   <- drop_kd1   | word(names(tab_kd1),i, sep = fixed(" | ")) == "Missing"
    drop_ks1_d <- drop_ks1_d | word(names(tab_ks1_d),i, sep = fixed(" | ")) == "Missing"
    drop_ksd1  <- drop_ksd1  | word(names(tab_ksd1),i, sep = fixed(" | ")) == "Missing"
  }
    Ndropident[1] <-  sum(tab_kd1[drop_kd1])
    Ndropident[2] <-  sum(tab_ks1[drop_ks1])
    Ndropident[3] <-  sum(tab_ks1_d[drop_ks1_d])
    Ndropident[4] <-  sum(tab_ksd1[drop_ksd1])
    tab_ks1[drop_ks1] <- 0
    tab_kd1[drop_kd1] <- 0
    tab_ks1_d[drop_ks1_d] <- 0
    tab_ksd1[drop_ksd1] <- 0
}

ks_anonym1_pct <- sum(tab_ks1)/Ns*100
kd_anonym1_pct <- sum(tab_kd1)/Nd*100
ks_anonym1_ind_pct <- sum(tab_ks1_d)/Nd*100
rep_uniques_pct <- sum(tab_ksd1)/Nd*100

ident[jj,] <- c( kd_anonym1_pct,ks_anonym1_pct, ks_anonym1_ind_pct,rep_uniques_pct  )
Nidentexclude[jj,] <- Ndropident


###------------------------get tables for  calculating attribute disclosure measures-----------------------------------------------------

nots <- apply(tab_kts_p,2,function(x) !(any(x==1))) ## indicates not unique in the margins of tab_kts

did <- tab_ktd ; did[tab_ktd_p != 1] <- 0
dis <- tab_kts ; dis[tab_kts_p != 1] <- 0
dis_in_orig <- tab_ktd
dis_in_orig[,nots]<- 0
dis_and_correct <- tab_ktd; 
dis_and_correct[tab_kts_p != 1] <- 0 ## disclosive and in original and correct
dis_and_correct_did <- dis_and_correct  ; dis_and_correct_did[tab_ktd_p != 1] <- 0 ## disclosive in syn in original  correct and dis in orig


Nout <- rep(0, 4)
names(Nout) = c("missing target","missing in  keys", "set to exclude","over denom_lim")
Nexcludes <- matrix(0,5,4)
dimnames(Nexcludes) <- list(c("DiO","DiS","DiSiO", "DiSCO", "DiSCO_DiO"), c("missing target","missing in  keys", "set to exclude","over denom_lim"))
###----------------------------- now exclusions for  for DiSiO DiSCO and DiSCO----------------------------------
for (jjj in 1:5)  {

if (jjj == 1)  xx <- did
if (jjj == 2)  xx <- dis
if (jjj == 3)  xx <- dis_in_orig
else if(jjj == 4) xx <-  dis_and_correct
else if(jjj == 5) xx <-  dis_and_correct_did

###------------------------ items excluded fOr missing values -------------------------
#

#####   first missings excluded values from tables if set 

 if (!usetargetNA && any(dd$target == "Missing")) {
   Nout[1] <- sum(xx[dimnames(xx)[[1]] == "Missing",])
   xx[dimnames(xx)[[1]] == "Missing",] <- 0
 }

 for (i in 1:length(keys)){
   if (!usekeysNA[i]) {   ## do not use NA values for ith key
     key_levs <- dimnames(xx)[[2]]
     drop_d <-  word(key_levs,i, sep = fixed(" | ")) == "Missing"

     Nout[2] <- Nout[2] + sum(xx[,drop_d])
     xx[,drop_d] <- 0
  }
}
###-------------------------- remove any excluded two way  combinations----------------------------------

 if (!is.null(exclude.keys)) {

  if (!all(exclude.targetlevs %in% levels(dd$target))) stop("exclude.targetlevs must be one of levels of ",target,"\n", call. =FALSE)
  for (i in 1:length(exclude.keys)){
    vout <- (1:(dim(xx)[1]))[dimnames(xx)[[1]] == exclude.targetlevs[i]]
    klev <- levels(dd[, names(dd) == exclude.keys[i]])
    if (!all(exclude.keylevs[i] %in% klev)) stop("exclude.keylevs position ",i, " must be one of levels of ",keys[i],"\n", call. =FALSE)#
    kind <- (1:length(keys))[keys == exclude.keys[i]]
    wordk <- word(dimnames(xx)[[2]], start = rep(kind, dim(xx)[2]), sep = fixed(" | "))
    kout <- (1:dim(tab_ktd)[2])[wordk == exclude.keylevs[i]]
    Nout[3] <- Nout[3] + sum(xx[vout, kout])
    xx[rep(vout,length(kout)), kout] <- 0
   }
 }  

###------------------ exclude if over denom_lim ----------------------------------
 if (exclude_ov_denom_lim)  {
   Nout[4] <- sum(xx[xx > denom_lim ])
   xx[xx > denom_lim ] <- 0
 }

 Nexcludes[jjj,] <- Nout
 if (jjj == 1)  did <- xx
 else if (jjj == 2)  dis <- xx
 else if (jjj == 3)  dis_in_orig <- xx
 else if (jjj == 4)   dis_and_correct <- xx
 else if (jjj == 5)   dis_and_correct_did <- xx

 }
###--------------------- exclusions done-----------------------------------------


DiO <- sum(did)/Nd*100
DiS <- sum(dis)/Ns*100
DiSiO <- sum(dis_in_orig)/Nd*100
DiSCO <- sum(dis_and_correct)/Nd*100
DiSCO_DiO <- sum(dis_and_correct_did)/Nd*100

denom_DisCO <- as.vector(tab_ktd[dis_and_correct >0 ])

attrib[jj,] <- c( DiO,DiS,DiSiO,DiSCO,DiSCO_DiO,max(dis_and_correct),mean(dis_and_correct[dis_and_correct>0]))
list_excludes[[jj]] <- Nexcludes

###------------------------get details for large denominators-----------------
if (any(dis_and_correct > denom_lim))  {

xx <- dis_and_correct

denoms <- as.vector(xx[xx > denom_lim])
rows <- rep(1:dim(xx)[1],dim(xx)[2])
rows <- rows[xx > denom_lim]
cols <- rep(1:dim(xx)[2],rep(dim(xx)[1],dim(xx)[2]))
cols <- cols[xx > denom_lim]

target_levels = dimnames(xx)[[1]][rows]
key_levs = dimnames(xx)[[2]][cols]

props <- tots <- max_denom <- key_levels <-  rep(NA,length(keys))
maxprops <- maxpropkey <-maxpropkeylevs <-maxpropdenom <- maxproptargetlevs <-  rep(NA,length(denoms))


  for ( j in 1:length(denoms) ) {
    for ( i in 1:length(keys)) {
    key_levels[i] <- word(key_levs[j],i, sep = fixed(" | "))
    tab <- table(dd$target, dd[,names(dd) == keys[i]])

    ### get proportions for 2-way tables of each part of key with target
    
    line <- tab[,dimnames(tab)[[2]] == key_levels[i]]
    props[i] <- (line[ names(line) == target_levels[j] ])/(sum(line))
    tots[i] <- sum(line)
    }

  maxprops[j] <- max(props)
  maxpropkey[j] <- keys[order(props)][length(keys)]
  maxpropkeylevs[j] <- key_levels[order(props)][length(keys)]
  maxpropdenom[j] <- tots[order(props)][length(keys)]

  }

  det <- data.frame( key = maxpropkey,
                        key_levels  = maxpropkeylevs, 
                        target_levels  = target_levels, 
                        maxprops = maxprops,
                        maxpropdenom = maxpropdenom,
                        denoms = denoms)
  det$target_key <- paste(det$target_levels,det$key_levels,sep = fixed(" | "))
  
  ###aggregate details by key and target combinations
   det <- det[order(det$target_key),]
   total <- tapply(det$denoms,det$target_key,sum)
   max_denom <- tapply(det$denoms,det$target_key,max)
   key <- det$key[!duplicated(det$target_key)]
   target_key <- det$target_key[!duplicated(det$target_key)]
   prop <- det$maxprops[!duplicated(det$target_key)]
   denom_2way <- det$maxpropdenom[!duplicated(det$target_key) ] 
                                  
   details <- data.frame( total = total,max_denom = max_denom, key = key, 
     target_key = target_key,  proportion_2way = prop, denom_2way = denom_2way )
   details <- details[order(-details$total),]
                        
   dimnames(details)[[1]] <- 1:dim(details)[1] 

 } else  details <-   NULL

list_denom[[jj]] <- details
list_excludes[[jj]] <- Nexcludes
###-------------------------- end of jj loop -----------------------------------
  }

res <- list(ident = data.frame(ident), attrib = data.frame(attrib), 
            allCAPs = data.frame(allCAPs), denom_details = list_denom, 
            Nidentexclude = Nidentexclude, list_excludes = list_excludes,
             keys = keys, target = target, digits = digits,  denom_lim = denom_lim, 
            usetargetNA = usetargetNA, usekeysNA = usekeysNA, 
            exclude.keys = exclude.keys, exclude.keylevs = exclude.keylevs, 
            exclude.targetlevs = exclude.targetlevs,  
            Norig = Nd,Nsyn = Ns,to.print = to.print,
            call = match.call())
class(res) <- "disclosure"
return(res)
  
  }

###---------------------------print.disclosure---------------------------- 
print.disclosure <- function(x, to.print = NULL, digits = NULL,   ...)
{
  if (is.null(to.print)) to.print <- x$to.print
  if (is.null(digits)) digits <- x$digits

  if (!all(to.print %in% c("short","allCAPs","ident","attrib",
                           "denom_details","all")))   
  stop("to.print must be  choices from 'short','allCAPs',
        'ident','attrib', 'denom_details','all'\n", call. = FALSE)

  if (length(to.print) == 1 && to.print == "short") {

        cat("Identity  measures for keys", x$keys,
    "\nand disclosure measures for",x$target,"from the same keys\n\n" )
    short <- rbind(c(x$ident[1,1],x$attrib[1,1]),
                   cbind(x$ident[,4],x$attrib[,4]))
    dimnames(short) <- list(c("Original",paste("Synthesis", 1:dim(x$attrib)[1])),c("Identity (repU)","Attrib (DiSCO)"))
    print(round(short, digits))
    if (any(x$attrib[,6] > x$denom_lim)) cat("Some denominators for DiSCO exceed denom_lim of", x$denom_lim, 
                                             "\ncheck the component  $denom_details of the dieclosure object.\n")
  }
  
  if ( any( c("allCAPs","ident","attrib","denom_details","all") %in% to.print )) {
     cat("\nDisclosure measures from synthesis of original data with ", x$Norig, "records")

  if ("all" %in% to.print) {
    cat("\n\nOutput from the following call to function disclosure()\n")
    print(x$call)
  }
    
    cat("\nIdentity  measures for keys", x$keys,
        "\nand disclosure measures for",x$target,"from the same keys\n\n" )  
    
    if (any(c("ident","all") %in% to.print))  { 
      cat("\n\nIdentity disclosure measures for", dim(x$ident)[[1]] ,"synthetic data set(s) from keys:\n", x$keys,"\n\n")
      print(round(x$ident, digits))
      if (any(x$Nidentexclude >0)) {
        cat("\nNumber of records excluded for missing keys\n")
        print(x$Nidentexclude)
      }
    }

    if (any(c("attrib","all") %in% to.print)) {
      cat("\n\nAttribute disclosure measures for",x$target,"from keys:",x$keys,"\n")
       print(round(x$attrib, digits))
       
       ###-------------------------- print exclusions for each synth if any------------------------
       
       for ( i in 1:dim(x$attrib)[1]) {
         if (any(x$list_excludes[[i]] > 0)) {
           cat("\nSynthesis", i,"Potentially disclosive records excluded from attribute disclosure\n")
           print(x$list_excludes[[i]])
         }
       }  

    }

    if (any(x$attrib[,6] > x$denom_lim) & !("denom_details" %in% to.print) )  {
      cat("\n\nCheck data because", sum(x$attrib[,6] > x$denom_lim),"of",dim(x$attrib)[1],
          "syntheses  have disclosive key/target combinations\n",
          "with denominators for DiSCO that exceed the set denom_lim of",x$denom_lim,
          "\nInclude 'denom_details' in to.print to see details.\n\n")
      }
    }
  ###-------------------------- print denom details for each synth-------------------------------
    if (any(c("denom_details","all") %in% to.print)) {
      
      if (any(x$attrib[,6] > x$denom_lim)) {
      cat("\nTarget-key combinations with the greatest contributions to large denominators\n\n")
      cat("For each disclosive set of keys the key that best predicts the target from a \n")
      cat("two-way table is identified and its two-way prior given.\n\n")

    
       for (i in 1:dim(x$attrib)[1]){
         if (x$attrib[i,6] > x$denom_lim) {
           cat("\nSynthesis", i,"\n") 
           print(x$denom_details[[i]])
         }
       }
      }     else
    cat("\nNo disclosive key/target combinations for any synthesis \n",
              "have denominators that exceed the set denom_lim of",x$denom_lim,"\n\n")
    } 
  
   
    if (any(c("allCAPs","all") %in% to.print)) {
    cat("\nCAP measures for the target", x$target, "with keys\n")
      cat(x$keys,"\n\n")
      print(round(x$allCAPs, digits))
    }   
  
  invisible(x)
  }
            