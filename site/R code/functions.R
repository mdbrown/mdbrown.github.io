
#function to perform stepwise forward regression for gee models
#this function assumes data is in the correct order for the gee function
# and that outcome is "Gleason.Progression"

stepforward_gee <- function(data, varlist, 
                            method = c("p-value", "QIC", "AUC"), 
                            threshold = 0.15, 
                            included.vars = NULL)
  {
    method = match.arg(method)
  
  
 
    keep.going = TRUE
  
    qic.last <- 1e6
    auc.last <- 0.5
    if(!is.null(included.vars) ) {
      tmp.formula <- as.formula(paste("Gleason.Progression", paste(included.vars, collapse = "+"), sep = "~"))
      #  cat(paste(tmp.formula)); cat("\n")
      tmpgee <- invisible(geeglm(tmp.formula, id = PASS.ID, family = binomial("logit"), 
                                 data = data ))
      risk.gee <- predict(tmpgee, type = "response")
      auc.last <-  roc(response = data$Gleason.Progression, predictor = risk.gee)$auc
      cat(auc.last); cat("\n")
    }
    while(keep.going){
      #keep adding variables to 'included.vars'
      #based on measure in metric list (p-value, delta AIC, delta AUC)
      metric.list <- rep(NA, length(varlist))
      i = 0
      for(tmpvar in varlist){
        i = i + 1
        
       tmp.formula <- as.formula(paste("Gleason.Progression", paste(append(included.vars, tmpvar), collapse = "+"), sep = "~"))
   #  cat(paste(tmp.formula)); cat("\n")
        tmpgee <- invisible(geeglm(tmp.formula, id = PASS.ID, family = binomial("logit"), 
                    data = data ))
      
        if(method == "p-value"){
          
          tmp.p <- coef(summary(tmpgee))[length(included.vars) + 2,4]
          metric.list[i] <- tmp.p
          
        }else if(method == "QIC"){
          metric.list[i] <- QIC(tmpgee)[1] 
        }else if(method == "AUC"){
         
          risk.gee <- predict(tmpgee, type = "response")
          metric.list[i] <- roc(response = data$Gleason.Progression, predictor = risk.gee)$auc - auc.last
        }
    }
  
    if(method == "QIC"){

      if(any(metric.list < qic.last)){
        
        chosen.ind <- which.min(metric.list)
        qic.last <- metric.list[chosen.ind]
        
        cat(paste("added:", varlist[chosen.ind], "\n"))
        included.vars <- append(included.vars, varlist[chosen.ind])
        varlist = varlist[-chosen.ind]
        
        
      }else{
        keep.going = FALSE
      }
      
    }else if(method == "p-value"){
      if(any(metric.list < threshold)){
        
         chosen.ind <- which.min(metric.list)
      

         cat(paste("added:", varlist[chosen.ind], "\n"))
         included.vars <- append(included.vars, varlist[chosen.ind])
         varlist = varlist[-chosen.ind]
   
      
      }else{
        keep.going = FALSE
      }
    }else if(method == "AUC"){
 
      if(any(metric.list > threshold)){
        
        chosen.ind <- which.max(metric.list)

        cat(paste("added:", varlist[chosen.ind], " AUC is:", round(auc.last + max(metric.list), 3), "\n"))
        auc.last <- metric.list[chosen.ind] + auc.last
        
        included.vars <- append(included.vars, varlist[chosen.ind])
        varlist = varlist[-chosen.ind]
        
        
      }else{
        keep.going = FALSE
      }
    }
  }
 
  return(included.vars)

}



###### 
# get folds for cross-validation
# stratify on outcome
# Y is binary outcome, V is number of folds. 
.cvFolds <- function(Y, V){
  # Create CV folds (stratify by outcome)   
  Y0 <- split(sample(which(Y==0)), rep(1:V, length=length(which(Y==0))))
  Y1 <- split(sample(which(Y==1)), rep(1:V, length=length(which(Y==1))))
  folds <- vector("list", length=V)
  for (v in seq(V)) {folds[[v]] <- c(Y0[[v]], Y1[[v]])}     
  return(folds)
}


#####
## output summary statistics to table
###

### calculate AUC in group A, B, and together. 
get.AB.AUC <- function(model, response, AB.grp, newdat = NULL){
  
  risk <- predict(model, newdata = newdat, type = 'response')
   
  #AUC 
  
  auc <- list( "both" = roc(predictor = c(risk), 
                            response = response), 
                "A" = roc(predictor = c(risk[AB.grp == "A"]), 
                          response = response[AB.grp == "A"]),
                "B" = roc(predictor = c(risk[AB.grp == "B"]), response = response[AB.grp == "B"]))

  list( "auc" = sapply(auc, FUN = function(x) x$auc), "full" = auc)

}
print.auc <- function(x, nround = 3){
  oo <- rbind((paste(round(x$auc[1], nround),
                     " (", round(x$auc.lower[1], nround), 
                     ",", round(x$auc.upper[1], nround), 
                     ")", sep = "")), 
              (paste( round(x$auc[2], nround),
                      " (", round(x$auc.lower[2], nround), 
                      ",", round(x$auc.upper[2], nround), 
                      ")", sep = "")),
              (paste( round(x$auc[3], nround),
                      " (", round(x$auc.lower[3], nround), 
                      ",", round(x$auc.upper[3], nround), 
                      ")", sep = "")))
  
  oo
}

get.AB.AUC_w.risk <- function(risk, response, AB.grp){
  
 # risk <- predict(model, newdata = newdat, type = 'response')
  
  #AUC 
  
  auc <- list( "both" = roc(predictor = c(risk), 
                            response = response), 
               "A" = roc(predictor = c(risk[AB.grp == "A"]), 
                         response = response[AB.grp == "A"]),
               "B" = roc(predictor = c(risk[AB.grp == "B"]), response = response[AB.grp == "B"]))
  
  list( "auc" = sapply(auc, FUN = function(x) x$auc), "full" = auc)
  
}


get.AB.AUC_withCI <- function(risk, response, AB.grp, bootstraps = 1000 ){
  
  auc.boot <- data.frame('both'= numeric(length = bootstraps), 
                         'A' = numeric(length = bootstraps),
                         'B' = numeric(length = bootstraps))
  
  for(b in 1:bootstraps){
    tmpind <- sample.int(length(risk), replace = TRUE)
    
    auc.boot[b,] <- get.AB.AUC_w.risk(risk[tmpind], response[tmpind], AB.grp[tmpind])$auc
    
  }
  auc.lower <- apply(auc.boot, 2,  quantile, c(0.025))
  auc.upper <- apply(auc.boot, 2, quantile, c(0.975))
 
  out <- get.AB.AUC_w.risk(risk, response, AB.grp)
  out$auc.lower <- auc.lower
  out$auc.upper <- auc.upper
  return(out)
  }
  
  

compare.AB.AUC_withCI <- function(risk1, risk2, response, AB.grp, bootstraps = 1000 ){
  
  auc.boot <- data.frame('both'= numeric(length = bootstraps), 
                         'A' = numeric(length = bootstraps),
                         'B' = numeric(length = bootstraps))

  for(b in 1:bootstraps){
    tmpind <- sample.int(length(risk1), replace = TRUE)
    
    auc.boot[b,] <- get.AB.AUC_w.risk(risk1[tmpind], response[tmpind], AB.grp[tmpind])$auc - get.AB.AUC_w.risk(risk2[tmpind], response[tmpind], AB.grp[tmpind])$auc
  }
  auc.lower <- apply(auc.boot, 2,  quantile, c(0.025))
  auc.upper <- apply(auc.boot, 2, quantile, c(0.975))
  
  est <- get.AB.AUC_w.risk(risk1, response, AB.grp)$auc - get.AB.AUC_w.risk(risk2, response, AB.grp)$auc
  out <- list(auc = est)
  out$auc.lower <- auc.lower
  out$auc.upper <- auc.upper
  return(out)
}



Summary.continue= function(colnm) 
{ dat.K$myvar<-dat.K[,colnm]

  junk = ddply(dat.K, .(pullSet, Gleason.Progression), summarise, round(mean(myvar,na.rm=TRUE),2), round(quantile(myvar, c(.25),na.rm=TRUE), 2),
               round(quantile(myvar, c(.75),na.rm=TRUE), 2),sum(is.na(myvar)))
  
  out = c(paste(junk[1,3],"(",junk[1,4],"-",junk[1,5],")",sep=''),
               paste(junk[2,3],"(",junk[2,4],"-",junk[2,5],")",sep=''),
               paste(junk[3,3],"(",junk[3,4],"-",junk[3,5],")",sep=''),
               paste(junk[4,3],"(",junk[4,4],"-",junk[4,5],")",sep=''))
  if(sum(is.na(dat.K$myvar))>0) {out = rbind(out,junk[,6])}
  out
  
}

Summary.category = function(colnm) {

  dat.K$myvar = dat.K[,colnm]
 
  myvar.val = eval(sort(unique(dat.K$myvar[!is.na(dat.K$myvar)])))
  np = length(myvar.val)
  out = matrix(0,np,4)
  ## this is clumsy but r will not work with change parameters: 
 if (np ==2) {
   out = t(ddply(dat.K, .(pullSet, Gleason.Progression),
                 summarise, 
                 sum(myvar==myvar.val[1], na.rm = TRUE), 
                 sum(myvar==myvar.val[2], na.rm = TRUE))[,3:4])
 }
 
  if(sum(is.na(dat.K$myvar))>0) {np = np+1; out = rbind(out,ddply(dat.K, .(pullSet, Gleason.Progression), summarise,sum(is.na(myvar)))[,3])}

 out.percent = round(out/VTM(apply(out,2,sum),np),2)
  out = round(out,2)
 junk = matrix(0,np,4)
  for (i in 1:(np)) {
   for (j in 1:4) {
     junk[i,j] = paste(out[i,j],"(",out.percent[i,j],")",sep="")
   }
  }
  junk 
}

VTM = function(vec,nr) {matrix(rep(vec,nr),nrow = nr)}

Summary.nbiopsy = function(colnm) {
  dat.K$myvar<-dat.K[,colnm]
  junk = ddply(dat.K, .(Biopsy.Number), summarise, round(mean(myvar, na.rm = TRUE),2), 
              round(quantile(myvar, c(.25), na.rm = TRUE), 2),round(quantile(myvar, c(.75), na.rm = TRUE), 2), 
              sum(is.na(myvar)), round(sum(is.na(myvar))/length(myvar),2))
  out=matrix(0,2,8)
  for (i in 1:8) {
    out[1,i] = paste(junk[i,2],'(', junk[i,3],'-', junk[i,4], ')',sep='')
    out[2,i] = paste(junk[i,5],'(', junk[i,6],')',sep='')
  }
  out  
}

Summary.nbiopsy.cat5 = function(colnm) {
  dat.K$myvar<-dat.K[,colnm]
  myvar.val = sort(unique(dat.K$myvar))
  if (sum(is.na(dat.K$myvar))>0) {
        junk = ddply(dat.K, .(Biopsy.Number), summarise, sum(myvar==myvar.val[1], na.rm = TRUE),
               sum(myvar==myvar.val[2], na.rm = TRUE),
               sum(myvar==myvar.val[3], na.rm = TRUE),
               sum(myvar==myvar.val[4], na.rm = TRUE),
               sum(myvar==myvar.val[5], na.rm = TRUE),
               sum(is.na(myvar)))  } else {
        junk = ddply(dat.K, .(Biopsy.Number), summarise, sum(myvar==myvar.val[1], na.rm = TRUE),
               sum(myvar==myvar.val[2], na.rm = TRUE),
               sum(myvar==myvar.val[3], na.rm = TRUE),
               sum(myvar==myvar.val[4], na.rm = TRUE),
               sum(myvar==myvar.val[5], na.rm = TRUE))
    }
    junk= junk[,-1]
    np = dim(junk)[2]
    junk.percent = round(junk/t(VTM(apply(junk,1,sum),np)),2)
    junk = t(junk)
    junk.percent = t(junk.percent)
  out=matrix(0,np,8)
  for (j in 1:8) {
    for (i in 1:np) {
      out[i,j] = paste(junk[i,j],"(",junk.percent[i,j],")",sep="")
    }
  }
  out  
}

#######
## leave one out cross-validation to get cross-validated estimates of risk
#######

leave.one.out.risk <- function(model, dat){
  
  risk <-numeric(nrow(dat))
  coef <- risk
  #for each individual (PASS.ID) I chuck them out, build a model, and estimate their risks
  ind.ids <- unique(dat$PASS.ID)
  
  for(id in ind.ids){
    #indices of id in dat
    dat.ind.myid <- which(dat$PASS.ID==id)
    
    tmp.fit <- lrm(model, x=TRUE, y=TRUE, data = dat[-dat.ind.myid,])
   # tmp.fit.adj <- robcov(tmp.fit, dat[-dat.ind.myid, "PASS.ID"]) #dont need var for this
               # cluster-adjusted Wald statistics
    # fastbw(g)         # cluster-adjusted backward elimination
    
 #   coef[dat.ind.myid] <- tmp.fit$coef[2]
    
    risk[dat.ind.myid] <- predict(tmp.fit, type = "fitted", newdata = dat[dat.ind.myid,])
    
  }
# browser()
 #tmp <- data.frame("coef" = coef, "outcome" = dat$Gleason.Progression)
 #ggplot(tmp, aes(coef, fill = outcome)) + geom_density(alpha = .5)
 
  risk
  
}

cv.fold.risk <- function(model, dat, nfold = 10){
  
  folds <- .cvFolds(Y = dat$Gleason.Progression, V = nfold)
  risk <-numeric(nrow(dat))
#browser()
#par(mfrow = c(2,5))
  for(f in 1:nfold){
    #indices of id in dat

    tmp.fit <- lrm(model, x=TRUE, y=TRUE, data = dat[-folds[[f]],])
    # tmp.fit.adj <- robcov(tmp.fit, dat[-dat.ind.myid, "PASS.ID"]) #dont need var for this
    # cluster-adjusted Wald statistics
    # fastbw(g)         # cluster-adjusted backward elimination
    
    #   coef[dat.ind.myid] <- tmp.fit$coef[2]
    
    risk[folds[[f]]] <- predict(tmp.fit, type = "fitted", newdata = dat[folds[[f]],])
    
   # plot(dat$Intact.PSA[folds[[f]]], risk[folds[[f]]])
  }
  # browser()
  #tmp <- data.frame("coef" = coef, "outcome" = dat$Gleason.Progression)
  #ggplot(tmp, aes(coef, fill = outcome)) + geom_density(alpha = .5)
  
  risk
  
}


