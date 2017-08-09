
myROC <- function(ftime, 
                  fstatus,
                  strata, 
                  times, 
                  data, 
                  span, 
                  strata.name){
  
  
  roc.dat <- NULL
  auc.dat <- NULL
  
  strata.levels <- unique(strata[!is.na(strata)])
  
  for(j in 1:length(strata.levels)){
    
    out <- NULL
    for(i in 1:length(times)){
      
      risk.ind <- which(names(data)==paste0("prob", times[i]))
      
      out2 = crCROC(Stime = ftime[is.element(strata, strata.levels[j])], 
             status = fstatus[is.element(strata, strata.levels[j])], 
             marker = data[is.element(strata, strata.levels[j]),risk.ind], 
             predict.time = times[i], 
             span = span)

      #out2$plot + ggtitle("2 year") + ylim(0, .1)
      out2$roc.dat <- data.frame(fpf = out2$FP, tpf = out2$TP[,1], year = paste(times[i], "year"))
      out2$roc.dat[[strata.name]] <- strata.levels[j]
      out2$FP <- NULL; out2$TP <- NULL; out2$cut.values = NULL; 
      out2$auc.dat <- data.frame(auc = out2$AUC[,1], year = paste(times[i], "year"))
      out2$auc.dat[[strata.name]] <- strata.levels[j]
      
      out[[i]] <- out2
    }
    
  
    roc.dat[[j]] <- bind_rows(lapply(out, function(x) x$roc.dat) )
    
    auc.dat[[j]] <- bind_rows(lapply(out, function(x) x$auc.dat))
    
  }
  
  roc.dat <- bind_rows(roc.dat)
  
  auc.dat <- bind_rows(auc.dat)
  
  
  roc.dat$year <- factor(roc.dat$year, levels = c("2 year", "5 year", "10 year"))
  auc.dat$year <- factor(auc.dat$year, levels = c("2 year", "5 year", "10 year"))
  
  return(list(roc.dat, auc.dat))
  
}

