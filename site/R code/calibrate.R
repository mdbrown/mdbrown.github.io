calibrationPlot <- function(ftime, fstatus, risk, groupNum = 10, time){
  
  require(cmprsk)


  risk.tiles <- c(0,quantile(risk, probs=c(1:groupNum/groupNum)[-groupNum], na.rm = TRUE), 1)
  
  risk.cut <- cut(risk, risk.tiles)
  
  groups = risk.cut
  
  tmpdat <- data_frame(ftime = ftime, fstatus = fstatus, groups = groups) 
  tmpdat <- tmpdat[complete.cases(tmpdat),]
  ftime = tmpdat$ftime; fstatus = tmpdat$fstatus; groups = tmpdat$groups
  #observed event rates from kaplan meier
  #cc <- CumIncidence(ftime = ftime, fstatus = fstatus, group = groups, t = c(2, 5, 10), level = .95)
  cuminc(ftime = ftime, fstatus = fstatus, cencode = 0)
  
  cc <- cuminc(ftime = ftime, fstatus = fstatus, group = groups, cencode = 0)
  obs.risk = timepoints(cc, times = time)
  
  level = .95
  z <- qnorm(1-(1-level)/2)
  lower <- obs.risk$est ^ exp(-z*sqrt(obs.risk$var)/(obs.risk$est*log(obs.risk$est)))
  upper <- obs.risk$est ^ exp(z*sqrt(obs.risk$var)/(obs.risk$est*log(obs.risk$est)))
  obs.risk$lower <- lower
  obs.risk$upper <- upper
  
  x <- lapply(obs.risk, c)
  obs.risk <- as_data_frame(x); rm(x)
  
  #average predicted risks by group
  
  obs.dat <- data.frame("observed" = obs.risk$est, 
                        "lower" = obs.risk$lower, 
                        "upper" = obs.risk$upper, 
                        "percentile" = c(c(1:groupNum/groupNum)[-groupNum], 1) - 1/(groupNum*2))[1:10,]
  
  risk = data.frame("risk" = t(t(risk)))
  risk$percentile = ecdf(risk$risk)(risk$risk)
  risk <- risk[order(risk$risk), ]
  risk.cc <- risk[complete.cases(risk),]
  
  p <- ggplot(obs.dat, aes(x=percentile, y=observed)) + 
    geom_point(data = obs.dat, size = 3) +  
    geom_errorbar(data = obs.dat, width = .025, aes(ymin = lower, ymax= upper)) + 
    geom_line(data = risk.cc, aes(y=risk))  + ylab("risk")
  
  
  p
  
  invisible(list("plot" =p, "risk" = risk.cc, "obs.dat" = obs.dat))
  
}

myCalibrate <- function(ftime, fstatus, times,strata, data,  strata.name){

obs.dat <- NULL
risk.dat <- NULL

strata.levels <- unique(strata[!is.na(strata)])

for(j in 1:length(strata.levels)){
  
  out <- NULL
  for(i in 1:length(times)){
   
    
    risk.ind <- which(names(data)==paste0("prob", times[i]))
    
    out2 <- calibrationPlot(ftime[is.element(strata, strata.levels[j])], 
                            fstatus[is.element(strata, strata.levels[j])], 
                            data[is.element(strata, strata.levels[j]),risk.ind],
                            10, times[i] )
    #out2$plot + ggtitle("2 year") + ylim(0, .1)
    out2$obs.dat <- cbind(out2$obs.dat, year = paste(times[i], "year"))
    out2$obs.dat[[strata.name]] <- strata.levels[j]
    out2$risk <- cbind(out2$risk, year = paste(times[i], "year"))
    out2$risk[[strata.name]] <- strata.levels[j]
    
    out[[i]] <- out2
  }

  obs.dat[[j]] <- bind_rows(lapply(out, function(x) x$obs.dat))
  
  risk.dat[[j]] <- bind_rows(lapply(out, function(x) x$risk))
  
}

obs.dat <- bind_rows(obs.dat)

risk.dat <- bind_rows(risk.dat)


obs.dat$year <- factor(obs.dat$year, levels = c("2 year", "5 year", "10 year"))
risk.dat$year <- factor(risk.dat$year, levels = c("2 year", "5 year", "10 year"))

return(list(obs.dat, risk.dat))

}
