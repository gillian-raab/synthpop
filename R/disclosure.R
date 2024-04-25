###-----new version with big denom checks---------------------
disclosure <- function(object, data, ...) UseMethod("disclosure")

###-----disclosure.default-----------------------------
disclosure.default <- function(object, ...)
  stop("No disclosure method associated with class ", class(object), call. = FALSE)

###-----disclosure.data.frame---disclosure.list--------
disclosure.data.frame <- disclosure.list <- 
  function(object, data, keys , target , denom_lim = 5,
           exclude_ov_denom_lim = FALSE, print.flag = TRUE, digits = 2, 
           usetargetNA = TRUE, usekeysNA = TRUE, not.targetlev = NULL,
           exclude.keys =NULL, exclude.keylevs = NULL, 
           exclude.targetlevs = NULL, 
           thresh_1way = c(50, 90),thresh_2way = c(5, 80),
           to.print =c("short"),  synorig.compare = FALSE, ...) 
    {
    if (is.null(object)) stop("Requires parameter 'object' to give name of the synthetic data.\n", call. = FALSE)   
    
    if (is.list(object) & !is.data.frame(object)) m <- length(object)
    else if (is.data.frame(object)) m <- 1
    else stop("object must be a data frame or a list of data frames.\n", call. = FALSE)
    
    if (synorig.compare) {
    
        if (m ==1) adjust.data <- synorig.compare(object,data, print.flag = FALSE) else
        if (m > 1) adjust.data <- synorig.compare(object[[1]],data, print.flag = FALSE)
        
        if (adjust.data$needsfix) stop("Synthetic data and/or original data needs more fixing before you can
      run the disclosure functions - see output. Use function synorig,compare() to check.", call. = FALSE)
        else if (!adjust.data$unchanged) {
          syn <- adjust.data$syn
          orig <- adjust.data$orig
          cat("Synthetic data or original or both adjusted with synorig.compare to try to make them comparable")
          if (m > 1) cat("only first element of the list has been adjusted and will be used here\n")
          m <- 1 }
        else if (print.flag) cat("Synthetic and original data checked with synorig.compare, no adjustment needed\n\n")
    }  

    object <- list(syn = object, m = m) 
    class(object) <- "synds"
    
    
    res <- disclosure(object, data, keys , target , denom_lim = denom_lim,
                      exclude_ov_denom_lim = exclude_ov_denom_lim, 
                      print.flag = print.flag, digits = digits, 
                      usetargetNA = usetargetNA, usekeysNA = usekeysNA,  
                      not.targetlev = not.targetlev,
                      exclude.keys = exclude.keys, exclude.keylevs =  exclude.keylevs,
                      exclude.targetlevs =  exclude.targetlevs,
                      thresh_1way = thresh_1way,thresh_2way = thresh_2way,
                      to.print = to.print, ...) 
    
    res$call <- match.call()
    return(res)
  }
    
###-----disclosure.synds-------------------------
disclosure.synds <-  function(object, data, keys , target , denom_lim = 5,
                              exclude_ov_denom_lim = FALSE, print.flag = TRUE, digits = 2, 
                              usetargetNA = TRUE, usekeysNA = TRUE, not.targetlev = NULL,
                              exclude.keys =NULL, exclude.keylevs = NULL, 
                              exclude.targetlevs = NULL, 
                              thresh_1way = c(50, 90),thresh_2way = c(5, 80),
                              to.print =c("short"), ...) 
  
  {
###-----------------------check input parameters ----
    if (!(is.data.frame(data)) )   stop("data  must be a data frame \n\n", call. = FALSE)
   data <- data.frame(data) ## in case it is table or a tibble
    if (!( inherits(object,"synds")))   stop(" object must be an object of class synds\n\n", call. = FALSE)
   
   if (is.numeric(keys) & !all(keys %in% 1:dim(data)[2])) stop("If keys are numeric they must be in range 1 to ",
                                      dim(data)[2] , call. = FALSE  )
   if (is.numeric(keys)) keys <- names(data)[keys]
   if (is.numeric(target) & !all(target %in% 1:dim(data)[2])) stop("If target is numeric it must be in range 1 to ",
                                                                   dim(data)[2] , call. = FALSE  )
   if (is.numeric(target)) target <- names(data)[target]
  
  Norig <- dim(data)[1]
 if (object$m ==1) {
   names.syn <- names(object$syn)
 }
   else names.syn <- names(object$syn[[1]])

  # target must be a single variable in  data and object$syn
  # keys must be a vector of variable names in data and in s
  # target must not be in keys


if (!(all(keys %in% names(data)) & all(keys %in% names.syn) & 
    all(target %in% names(data)) & all(target %in% names.syn)))
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
  ident  <- matrix(NA,object$m ,4)

  dimnames(allCAPs) <- list(1:object$m, c("baseCAPd","CAPd", "CAPs", "DCAP","TCAP"))
  dimnames(attrib) <- list(1:object$m,c("DiO","DiS","KDiSiO","DiSDiO", "DiSCO", "max_denom","mean_denom"))
  dimnames(ident) <- list(1:object$m,c("UiO", "UiS","UiOiS", "repU"))
  check_2way <-list(1:object$m)
  Nexclusions <- list(1:object$m)

  
  dd <- data ## rename target variable to target
  targetx =dd[, names(dd) == target]
  dd$target <- targetx
 # cat(names(dd),"names(dd)\n")
  dd <- as.data.frame(dd[,names(dd) %in% c( keys, "target")])


  if ( !is.factor(dd$target) )   dd$target <- factor(dd$target)
  for (i in 1:length(keys))   dd[,names(dd) == keys[i]] <- factor( dd[,names(dd) == keys[i]])
  
  check1way <- totalDisclosive <- PctDisAll <- totalLevel <- pctdisLevel <-  totLevAll <- LevpctAll <- rep(0, object$m)
  most_dis_lev <- rep("",object$m)
  
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

###----------------------- sort missings-----------------     
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
###------------------------ make composite variable for keys --------------------------
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
 
###------------------------- calculate proportions and margins ---------
  
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
###------------------------get tables for  calculating attribute disclosure measures-----------------------------------------------------

nots <- apply(tab_kts_p,2,function(x) !(any(x==1))) 
  
  ## indicates not unique in the margins of tab_kts


did <- tab_ktd ; did[tab_ktd_p != 1] <- 0
dis <- tab_kts ; dis[tab_kts_p != 1] <- 0
keys_orig <- apply(tab_ktd,2,sum)
keys_in_syn_in_orig <- dis 
keys_in_syn_in_orig[,keys_orig == 0] <- 0
dis_in_orig <- keys_in_syn_in_orig
dis_in_orig[tab_ktd == 0] <- 0 ## disclosive and in original 
dis_both_correct <- dis_in_orig  
dis_both_correct[tab_ktd_p != 1] <- 0 ## disclosive in syn in original  correct and dis in orig


Nout <- rep(0, 6)
names(Nout) = c("excluded target","missing target","missing in  keys", "set to exclude","over denom_lim","remaining")
Nexcludes <- matrix(0,7,6)
dimnames(Nexcludes) <- list(c("original","synthetic","DiO","DiS","KDiSiO", "DiSDiO", "DiSCO"),
                            c("excluded target","missing target","missing in  keys", "set to exclude","over denom_lim","remaining"))
###----------------------------- now exclusions----------------------------------

tab_exclude <- function(xx,col,Nexcludes) {
 total <- sum(xx)
 ###------------------- drop all target records with not.targetlev---------------------------------
 if (!is.null(not.targetlev) && !(not.targetlev == ""))
  {
   Nout[1] <- sum(xx[dimnames(xx)[[1]] %in% not.targetlev,])
   xx[dimnames(xx)[[1]] %in% not.targetlev,] <- 0
 }
  ###------------------------ items excluded fOr missing values -------------------------
  #
##### missings excluded values from tables if set

if (!usetargetNA && any(dd$target == "Missing")) {
  print(dimnames(xx)[[1]] )
  Nout[1] <- sum(xx[dimnames(xx)[[1]] == "Missing",])
  xx[dimnames(xx)[[1]] == "Missing",] <- 0
}

 for (i in 1:length(keys)){
   if (!usekeysNA[i]) {   ## do not use NA values for ith key
     key_levs <- dimnames(xx)[[2]]
     drop_d <-  word(key_levs,i, sep = fixed(" | ")) == "Missing"
     Nout[3] <- Nout[3] + sum(xx[,drop_d])
     xx[,drop_d] <- 0
  }
}
###-------------------------- remove any excluded two way  combinations----------------------------------

 if (!is.null(exclude.keys)) {

  if (!all(exclude.targetlevs %in% levels(dd$target))) stop("exclude.targetlevs must be one of levels of ",target,"\n", call. =FALSE)
  
   for (i in 1:length(exclude.keys)){
    vout <- (1:(dim(xx)[2]))[dimnames(xx)[[2]] == exclude.targetlevs[i]]
    klev <- levels(dd[, names(dd) == exclude.keys[i]])
    if (!all(exclude.keylevs[i] %in% klev)) stop("exclude.keylevs position ",i, " must be one of levels of ",keys[i],"\n", call. =FALSE)#
    kind <- (1:length(keys))[keys == exclude.keys[i]]
    wordk <- word(dimnames(xx)[[2]], start = rep(kind, dim(xx)[2]), sep = fixed(" | "))
    kout <- (1:dim(tab_ktd)[2])[wordk == exclude.keylevs[i]]
    Nout[4] <- Nout[4] + sum(xx[vout, kout])
    xx[rep(vout,length(kout)), kout] <- 0
   }
 }  

###------------------ exclude if over denom_lim ----------------------------------
 if (exclude_ov_denom_lim)  {
   Nout[5] <- sum(xx[xx > denom_lim ])
   xx[xx > denom_lim ] <- 0
 }
 Nout[6] <-  total -sum(Nout[1:5])
 
###---------------------------- copy Nout into rows of Nexcludes ------------------------------
 Nexcludes[col,] <- Nout
   return(list(tab = xx, Nexcludes = Nexcludes))
}
tab_ktdO <- tab_ktd
yy <- tab_exclude(tab_ktd,1,Nexcludes)
tab_ktd <- yy$tab; Nexcludes <- yy$Nexcludes
yy <- tab_exclude(tab_kts,2,Nexcludes)
tab_kts <- yy$tab; Nexcludes <- yy$Nexcludes
yy <- tab_exclude(did,3,Nexcludes)
did <- yy$tab; Nexcludes <- yy$Nexcludes
yy <- tab_exclude(dis,4,Nexcludes)
dis <- yy$tab; Nexcludes <- yy$Nexcludes
yy <- tab_exclude(keys_in_syn_in_orig,5,Nexcludes)
keys_in_syn_in_orig <- yy$tab; Nexcludes <- yy$Nexcludes
yy <- tab_exclude(dis_in_orig,6,Nexcludes)
dis_in_orig <- yy$tab; Nexcludes <- yy$Nexcludes
yy <- tab_exclude(dis_both_correct,7,Nexcludes)
dis_both_correct <- yy$tab; Nexcludes <- yy$Nexcludes
###--------------------- exclusions done-----------------------------------------
###-----------------------get  identity disclosure measures -------------
t1 <- apply(did,2,sum)
tab_ks <- apply(tab_kts,2,sum)
tab_kd <- apply(tab_ktd,2,sum)
tab_ks1 <- tab_ks
tab_ks1[tab_ks1>1] <- 0
tab_kd1 <- tab_kd 
tab_kd1[tab_kd1>1] <- 0
tab_kd1_s <- tab_kd1[names(tab_kd1) %in% Ks]
tab_ksd1 <-tab_kd[tab_ks == 1 & tab_kd == 1]
##print(tab_ktdO[,(t1-tab_kd1)<0][,1:10])
## exclude bits with missing levels of keys


UiS<- sum(tab_ks1)/Ns*100
UiO<- sum(tab_kd1)/Nd*100
UiOiS<- sum(tab_kd1_s)/Nd*100
repU <- sum(tab_ksd1)/Nd*100

ident[jj,] <- c( UiO,UiS, UiOiS,repU )
###----------------------------- attrib dis measures-------------------------
DiO <- sum(did)/Nd*100
DiS <- sum(dis)/Ns*100
KDiSiO <- sum(keys_in_syn_in_orig)/Nd*100
DiSDiO <- sum(dis_in_orig)/Nd*100
DiSCO <- sum(dis_both_correct)/Nd*100

denom_DisCO <- as.vector(tab_ktd[dis_both_correct])

attrib[jj,] <- c( DiO,DiS,KDiSiO,DiSDiO,DiSCO,max(dis_both_correct),mean(dis_both_correct[dis_both_correct>0]))
Nexclusions[[jj]] <- Nexcludes

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

###----------------------- checks for most_dis_lev 1 way ---------------------------
tab_target <- apply(tab_ktd,1,sum)
tab_dis_target <- apply(dis_in_orig,1,sum)
pctdisLev <- max(tab_dis_target)/sum(tab_dis_target)*100
 
totalDisclosive[jj] <- sum(tab_dis_target)
PctDisAll[jj] <- totalDisclosive[jj]/sum(tab_target)*100
tab_dis_target <- sort(tab_dis_target, decreasing = TRUE)
totalLevel[jj] <-tab_dis_target[1]


  most_dis_lev[jj] <- names(tab_dis_target)[1]
  totalLevel[jj] <-tab_dis_target[1]
  pctdisLevel[jj] <- tab_dis_target[1]/sum(dis_in_orig)*100

  if (sum(tab_dis_target) > thresh_1way[1] && pctdisLev >=  thresh_1way[2]) check1way[jj] <- 1

###------------------------get details for check_2way-----------------
### only implemented for DiSCO  though could be changed
  
  xx <- dis_both_correct
  if (any(xx > thresh_2way[1]))  {
    
   denoms <- as.vector(xx[xx > thresh_2way[1]])
  rows <- rep(1:dim(xx)[1],dim(xx)[2])
  rows <- rows[xx > thresh_2way[1]]
  cols <- rep(1:dim(xx)[2],rep(dim(xx)[1],dim(xx)[2]))
  cols <- cols[xx > thresh_2way[1]]
  

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
   pctall <- total*100/sum(xx)
   max_denom <- tapply(det$denoms,det$target_key,max)
   key <- det$key[!duplicated(det$target_key)]
   target_key <- det$target_key[!duplicated(det$target_key)]
   prop <- det$maxprops[!duplicated(det$target_key)]
   denom_2way <- det$maxpropdenom[!duplicated(det$target_key) ] 
                                  
   details <- data.frame( total = total,pctall = round(pctall,1), max_denom = max_denom, key = key, 
   target_key = target_key,  pct_2way = round(prop*100,1), denom_2way = denom_2way )
   if (all(details$pct_2way <= thresh_2way[2])) details <- NULL
   else {
     details <- details[details$pct_2way > thresh_2way[2],]
     details <- details[order(details$total, decreasing = TRUE),]
     dimnames(details)[[1]] <- 1:(dim(details)[1]) 
   }
   check_2way[[jj]] <- details
  }
  else check_2way[[jj]] <- NULL


###-------------------------- end of jj loop -----------------------------------
  }
###----------------------- checks for most_dis_lev 1 way ---------------------------
 
check_1way <- data.frame(check= check1way, Level = most_dis_lev, totalDisclosive =  totalDisclosive,
                        PctDisAll = PctDisAll, totalLevel = totalLevel,
                        pctdisLevel = pctdisLevel)
if (any (check_1way$check ==1 )) {
     if (object$m >1 ) {
       check_1way <- check_1way[check_1way$check ==1,]
       ulevs <- unique(check_1way$Level)
       nlevs <- length(ulevs)
       newcheck <- check_1way[1:nlevs,]
       newcheck$Level <- ulevs
       for (i in c(3:6)) newcheck[,i] <- tapply(check_1way[,i],check_1way[,2], mean)
       check_1way <- newcheck[order(-newcheck[,3]),-1]  
       }
     else check_1way <- check_1way
     checklev_1way <- paste(unique(check_1way$Level), sep = "|")
}
else {
  check_1way <- NULL
  checklev_1way <- ""
  }
  
if (length(check_2way) == 0) {
  check_2way <- NULL
  checklevs_2way <- ""
}
else {
  alllevs <- paste(check_2way[[1]]$target_key,"for key",check_2way[[1]]$key)
  if (object$m > 1) {
    for (i in 2:object$m){
      allevs <- c(alllevs,paste(check_2way[[i]]$target_key,"for key",check_2way[[i]]$key))
    }
  }
  alllevs <- unique(alllevs)
  checklevs_2way <- (alllevs)
}


res <- list(ident = data.frame(ident), attrib = data.frame(attrib), 
            allCAPs = data.frame(allCAPs), check_1way = check_1way, 
            checklev_1way = checklev_1way, check_2way = check_2way,
            checklevs_2way =checklevs_2way, Nexclusions = Nexclusions,
            keys = keys, target = target, digits = digits, Norig = Norig,
            to.print = to.print, call = match.call())
class(res) <- "disclosure"
return(res)
}
###---------------------------print.disclosure---------------------------- 
print.disclosure <- function(x, to.print = NULL, digits = NULL,   ...)
{
  if (is.null(to.print)) to.print <- x$to.print
  if (is.null(digits)) digits <- x$digits

  if (!all(to.print %in% c("short","allCAPs","ident","attrib",
                           "check_2way","check_1way","all","exclusions")))   
  stop("to.print must be  choices from 'short','allCAPs',
        'ident','attrib', 'check_2way','check_1way','all','exclusions'\n", call. = FALSE)
  
  if (length(to.print) >1 & ("short" %in% to.print)) stop("A 'Short' entry in to.print should not be combined with other options.\n", call. = FALSE)
  
  
  cat("Disclosure measures from synthesis for",x$Norig, "records in original data")

  
  nexcluded <- sapply(x$Nexclusions,function(x) sum(x[,-6]))
  if  (any(nexcluded>0)) {
    cat("\nSome records excluded from the evaluation of disclosure risks")
    cat("\nSee details by adding 'exclusions' to the parameter to.print of disclosure,
    or by printing the details with to,print including 'exclusions'.")
  }
  
  if (length(to.print) == 1 && to.print == "short") {
    cat("\nIdentity  measures for keys", x$keys,
        "\nand attribute measures for",x$target,"from the same keys" )
    short <- rbind(c(x$ident[1,1],x$attrib[1,1]),
                   cbind(x$ident[,4],x$attrib[,5]))
    dimnames(short) <- list(c("Original",paste("Synthesis", 1:dim(x$attrib)[1])),c("Identity (UiO/repU)","Attrib (DiO/DiSCO)"))
    print(round(short, digits))

  }
  if ("all" %in% to.print) {
    cat("\nOutput from the following call to function disclosure()")
    print(x$call)
  }
    
    if (any(c("ident","all") %in% to.print))  { 
      cat("\nIdentity disclosure measures for", dim(x$ident)[[1]] ,"synthetic data set(s) from keys:\n", x$keys,"\n\n")
      print(round(x$ident, digits))
      }
    
    if (any(c("attrib","all") %in% to.print)) {
      cat("\nAttribute disclosure measures for",x$target,"from keys:",x$keys,"\n")
       print(round(x$attrib, digits))
    }
if (any(x$Nidentexclude >0)) {
  cat("\nNumber of records excluded from identity disclosure for missing keys\n")
  print(x$Nidentexclude)  
}
###---------------- print CAP measures-----------------------------------------
  if (any(c("allCAPs","all") %in% to.print)) {
    cat("\nCAP measures for the target", x$target, "with keys\n")
    cat(x$keys,"\n\n")
    print(round(x$allCAPs, digits))
  }   
  
  ###------------------------ print check messages  
  if (x$checklev_1way != ""){
    cat("\nThe 1 way distributions of",x$target,"has a large contribution to disclosure from level",x$checklev_1way,"\n")
    cat("Please add 'check_1way' to to.print (e.g. print(disclosure_result, to.print = 'check_1way'),\n")
    cat("and look at original data to decide if backround knowledge would make this disclosure likely\n")
    cat("Consider excluding this level with the not.targetlev parameter to the disclosure function.\n\n")
  }
       
  if ("check_1way" %in% to.print){
      if (!is.null(x$check_1way)) {
      cat("\nDetails of target contributing disproportionately to disclosure\n")
      if (dim(x$check_1way)[1] >1) cat("averaged over all",dim(x$check_1way)[1],"syntheses\n")
      print(x$check_1way)
      }
      else {cat("No 1 way checks for your target are flagged at your settings for thresh_1way\n\n ")}
  }
  
  if (!(length(x$checklevs_2way) == 1 && x$checklevs_2way == "")){
    cat("\nDetails of target-key pairs contributing disproportionately to disclosure of",x$target,"\n")
    if (length(x$checklevs_2way) > 3) cat("Note only the first 3 printed here\n")
    for (i in 1:(min(length(x$checklevs_2way),3))) cat(x$checklevs_2way[i],"\n")
    cat("\nPlease examine component $check_2way of the disclosure object and look at original data.\n")
    cat("Consider excluding these key-target pairs with some the following parameters to disclosure:\n")
    cat("exclude_ov_denom_lim = TRUE or defining key-target combinations from exclude.targetlevs,\n")
    cat("exclude.keys and exclude.keylevs\n\n" )
  }
  
  if ("check_2way" %in% to.print){
      cat("\nDetails of target-key combinations contributing disproportionately to disclosure\n")
      if (!is.null(x$check_2way)) {
        cat("This is a list with one for each of",dim(x$check_2way),"syntheses\n")
        print(x$check_2way)
       }
      else {cat("No 2 way checks of key-target combinations are flagged at your settings for thresh_2way\n\n ")}
    }
  ###----------------------  print exclusions -----------------
  if ("exclusions" %in% to.print){
    for (i in 1:length(x$Nexclusions)){
      if (sum(x$Nexclusions[[i]]) > 0 ) {
        cat("\nRecords excluded from attribute exclusion measures for synthesis",i,"\n") 
        print(x$Nexclusions[[i]])}
      else cat("\nNo records excluded from attribute exclusion measures for synthesis",i,"\n")
    }
  }

  
  invisible(x)
    }
  