# forecasting exercise citi Bike demand NY   ------------------------------

#Maastricht University
#School of Business and Economics
#Maastricht, October 2026
#Vyvey Declerck A., Vogelaar Q., van Hooy J.
#I6267153, I6344964, I6367300
#EBC4008 – Time series analytics
#Tutorial group 03
#29 September 2026
#Dr. Wilms I.
#Course term project
#a.vyveydeclerck@student.maastrichtuniversity.nl
#q.vogelaar@student.maastrichtuniversity.nl
#j.vanhooy@student.maastrichtuniversity.nl

# setup -------------------------------------------------------------------
rm(list=ls())
set.seed(299792458)
#packages
library(ggplot2)
library(forecast)
library(bootUR)
library(MCS)
library(smooth)
library(Matrix)
library(FoReco)
#data
load(file.choose())
data<-data_citibike
data1<-data

# functions ---------------------------------------------------------------
summary1<-function(x){
  st1<-round(mean(x),4)
  st2<-round(var(x),4)
  st3<-round(max(x),4)
  st5<-round(min(x),4)
  print("summary of the data")
  print(paste("the mean is:", st1))
  print(paste("the variance is:", st2))
  print(paste("the maximum is:", st3))
  print(paste("the minimum is:", st5))
}
# following functions determine optimal ARMA model based on different IC
# the functions adjust for changing sample size to make comparison possible
# the IC are stored in a matrix with P and Q as row and column
autolaglengthAIC<-function(p,q,y){
  m1<-max(p,q)
  y1<-if (m1 == 0) y else y[-(1:m1)]
  matrixaic<-matrix(nrow = p+1,ncol = q+1)
  for(i in 0:p){
    for(j in 0:q){
      M1<-Arima(y1,order = c(i,0,j))
      matrixaic[i+1,j+1]<-M1[["aicc"]]
    }
  }
  
  inds1 = which(matrixaic == min(matrixaic), arr.ind=TRUE)-1
  
  print("lag length according to AICc (col=q row =p:)")
  return(inds1)
} # built during quantitative techniques for fin. econ.
autolaglengthSC<-function(p,q,y){
  m1<-max(p,q)
  y1<-if (m1 == 0) y else y[-(1:m1)]
  t<-length(y1)
  matrixsc<-matrix(nrow = p+1,ncol = q+1)
  for(i in 0:p){
    for(j in 0:q){
      M1<-Arima(y1,order = c(i,0,j)) 
      a<-M1[["residuals"]]
      SC<-(i+j+1)*log(t)/t+log(mean(a^2))
      matrixsc[i+1,j+1]<-SC
    }
  }
  
  inds2 = which(matrixsc == min(matrixsc), arr.ind=TRUE)-1
  print("lag length according to SC (col=q row =p:)")
  
  return(inds2)
}# built during quantitative techniques for fin. econ.
autolaglengthAIC2<-function(p,q,y,xreg){
  xreg<-as.matrix(xreg)
  m1<-max(p,q)
  if (m1 == 0) {
    y1<-y
    xreg1<-xreg
  } else {
    y1<-y[-(1:m1)]
    xreg1<-xreg[-(1:m1), , drop=FALSE]
  }
  matrixaic<-matrix(nrow = p+1,ncol = q+1)
  for(i in 0:p){
    for(j in 0:q){
      M1<-Arima(y1,order = c(i,0,j),xreg = xreg1)
      matrixaic[i+1,j+1]<-M1[["aicc"]]
    }
  }
  
  inds1 = which(matrixaic == min(matrixaic), arr.ind=TRUE)-1
  
  print("lag length according to AICc (col=q row =p:)")
  return(inds1)
}# built during quantitative techniques for fin. econ.
autolaglengthSC2<-function(p,q,y,xreg){
  xreg<-as.matrix(xreg)
  m1<-max(p,q)
  if (m1 == 0) {
    y1<-y
    xreg1<-xreg
  } else {
    y1<-y[-(1:m1)]
    xreg1<-xreg[-(1:m1), , drop=FALSE]
  }
  t<-length(y1)
  matrixsc<-matrix(nrow = p+1,ncol = q+1)
  for(i in 0:p){
    for(j in 0:q){
      M1<-Arima(y1,order = c(i,0,j),xreg=xreg1) 
      a<-M1[["residuals"]]
      SC<-(i+j+1)*log(t)/t+log(mean(a^2))
      matrixsc[i+1,j+1]<-SC
    }
  }
  
  inds2 = which(matrixsc == min(matrixsc), arr.ind=TRUE)-1
  print("lag length according to SC (col=q row =p:)")
  
  return(inds2)
}# built during quantitative techniques for fin. econ.
# evaluation of forecasting errors
in.sample.criteria<-function(resid){
  resid<-resid[!is.na(resid)]
  x1<-(-sum(resid[resid < 0]))
  x2<-(-mean(resid[resid < 0]))
  x3<-sum(resid < 0)
  rmse<-(mean(resid^2))^0.5
  lad<-mean(abs(resid))
  print(paste("the total fitted overshoot is:", x1))
  print(paste("the mean of fitted overshoot is:", x2))
  print(paste("the number of fitted overshoot is:", x3))
  print(paste("the RMSE is:", rmse))
  print(paste("the MAE is:", lad))
}

# the models used in the forecast exercise
# arima() occasionally fails to converge on a single window; it is then re-estimated
# with exact maximum likelihood (method = "ML")
# (arima() is called directly inside each function, because predict() needs to find
# the xreg of the original call)
# each function gets: y = estimation sample (ts), xr = regressor in the sample,
# newxr = regressor over the forecast period, steps = number of steps ahead
hourly_models<-list(
arima_2_0_4_notrend=function(y,xr,newxr,steps){
    fit<-tryCatch(
      arima(y, order=c(2,0,4), include.mean=FALSE),
      error=function(e) arima(y, order=c(2,0,4), include.mean=FALSE, method="ML")
    )
    predict(fit, n.ahead=steps)$pred
  },
arima_3_0_5_t=function(y,xr,newxr,steps){
  fit<-tryCatch(
    arima(y,order=c(3,0,5),xreg=xr),
    error=function(e) arima(y,order=c(3,0,5),xreg=xr,method="ML")
  )
  predict(fit,n.ahead=steps,newxreg=newxr)$pred
},
arima_2_0_4_t=function(y,xr,newxr,steps){
    fit<-tryCatch(
      arima(y,order=c(2,0,4),xreg=xr),
      error=function(e) arima(y,order=c(2,0,4),xreg=xr,method="ML")
    )
    predict(fit,n.ahead=steps,newxreg=newxr)$pred
  },
# daily SARIMA(2,0,2)(1,0,1)[24] + trend (the auto.arima choice in-sample, c1):
# orders fixed, parameters re-estimated every day
sarima_2_0_2_101_t=function(y,xr,newxr,steps){
    fit<-tryCatch(
      arima(y,order=c(2,0,2),seasonal=list(order=c(1,0,1),period=frequency(y)),xreg=xr),
      error=function(e) arima(y,order=c(2,0,2),seasonal=list(order=c(1,0,1),period=frequency(y)),xreg=xr,method="ML")
    )
    predict(fit,n.ahead=steps,newxreg=newxr)$pred
  },
msarima_102_003_103 = function(y,xr,newxr,steps){
  
  fit <- smooth::msarima(
    y,
    orders=list(
      ar=c(1,0,1),
      i=c(0,0,0),
      ma=c(2,3,3)),
    lags=c(1,6,24),
    regressors="use",
    initial="backcasting")
  
  forecast(fit, h=steps)$mean
},
ets_MAM=function(y,xr,newxr,steps){
    forecast(ets(y,model="MAM"),h=steps)$mean
  },
# weekly SARIMA (smooth package, lags 1 and 168) with the orders selected in-sample
# (weekly_orders, from e2 in the weekly seasonality section): orders fixed, parameters
# re-estimated every day. constant = TRUE: drift (the trend after the weekly difference)
sarima_weekly=function(y,xr,newxr,steps){
    fit<-smooth::msarima(
      ts(as.numeric(y),frequency=168),
      orders=weekly_orders,
      lags=c(1,168),
      constant=TRUE,
      initial="backcasting")
    as.numeric(forecast(fit,h=steps)$mean)
  })

daily_models<-list(
arima_1_0_0=function(y,xr,newxr,steps){
    fit<-tryCatch(
      arima(y,order=c(1,0,0),xreg=xr),
      error=function(e) arima(y,order=c(1,0,0),xreg=xr,method="ML")
    )
    predict(fit,n.ahead=steps,newxreg=newxr)$pred
  },
hw_damped=function(y,xr,newxr,steps){
    hw(y,damped=TRUE,h=steps)$mean
  },
ets_ANA=function(y,xr,newxr,steps){
    forecast(ets(y,model="ANA"),h=steps)$mean
  })

# Rolling or expanding window forecasts for horizons h = 1,...,H
# - F.start: position of the first observation of the test sample (first forecast target)
# - every model is estimated ONCE per forecast origin and forecasts H steps ahead;
#   the j-th step of that forecast is the h = j forecast from that origin
#   (same forecasts as refitting per horizon, but H times fewer model estimations)
# - rolling vs expanding differs in ONE line: the start of the estimation sample
# - origin.step = 1 uses every hour of the test sample as forecast origin;
#   origin.step = 24 would use one origin per day (much faster)
# - a model that fails on a window gets NA for that origin; the evaluation drops those
window.forecast<-function(x,exo,freq.,F.start,H,models,window=c("rolling","expanding"),origin.step=1){
  window<-match.arg(window)
  n<-length(x)
  win.length<-F.start-1                      # rolling window keeps the training sample length
  origins<-seq(F.start,n,by=origin.step)     # first forecast target of each origin
  n.origins<-length(origins)
  
  # fc[[m]][i,j] = forecast of model m, made at origin i, for horizon j
  fc<-lapply(models,function(m) matrix(NA,nrow=n.origins,ncol=H))
  
  for(i in 1:n.origins){
    end.s<-origins[i]-1
    start.s<-if (window == "rolling") end.s-win.length+1 else 1
    steps<-min(H,n-end.s)
    
    e.sample<-ts(x[start.s:end.s],frequency = freq.)
    exo1<-exo[start.s:end.s]
    future.exo<-exo[(end.s+1):(end.s+steps)]
    
    for(m in names(models)){
      fc[[m]][i,1:steps]<-tryCatch(
        as.numeric(models[[m]](e.sample,exo1,future.exo,steps))[1:steps],
        error=function(e) rep(NA,steps)
      )
    }
    
    if (i %% 24 == 0) cat(window,"window: origin",i,"of",n.origins,"\n")
  }
  
  # reorganise per horizon: errors = actual - forecast
  out<-lapply(1:H,function(h){
    valid<-which(origins-1+h <= n)             # forecast target must lie in the data
    target<-origins[valid]-1+h
    forecasts<-data.frame(lapply(fc,function(M) M[valid,h]))
    errors<-x[target]-forecasts
    list(errors=errors,forecasts=forecasts,target=target)
  })
  names(out)<-paste0("horizon_",1:H)
  
  n.failed<-sapply(fc,function(M) sum(is.na(M[,1])))
  if (any(n.failed > 0)) {
    print("number of origins where a model could not be estimated:")
    print(n.failed)
  }
  
  out
}

# same criteria as in.sample.criteria(), but returned instead of printed,
# so they can be collected in a table per horizon and model
forecast.criteria<-function(resid){
  resid<-resid[!is.na(resid)]
  c(
    total_overshoot = -sum(resid[resid < 0]),
    mean_overshoot  = if (any(resid < 0)) -mean(resid[resid < 0]) else 0,
    n_overshoot     = sum(resid < 0),
    RMSE            = sqrt(mean(resid^2)),
    MAE             = mean(abs(resid))
  )
}

# long table: one row per horizon x model
criteria.by.horizon<-function(forecast_list){
  rows<-list()
  for(h in seq_along(forecast_list)){
    errors_h<-forecast_list[[h]][["errors"]]
    for(m in colnames(errors_h)){
      rows[[length(rows)+1]]<-data.frame(
        horizon=h,
        model=m,
        n_forecasts=sum(!is.na(errors_h[[m]])),
        t(forecast.criteria(errors_h[[m]]))
      )
    }
  }
  out<-do.call(rbind,rows)
  rownames(out)<-NULL
  out
}

# wide table of one criterion: rows = horizon, columns = models
metric.table<-function(crit,metric){
  wide<-reshape(
    crit[,c("horizon","model",metric)],
    idvar="horizon",
    timevar="model",
    direction="wide"
  )
  colnames(wide)<-sub(paste0(metric,"."),"",colnames(wide),fixed=T)
  rownames(wide)<-NULL
  wide
} 
plot.metric<-function(crit,metric,title){
  ggplot(crit,aes(x=horizon,y=.data[[metric]],linetype=model))+
    geom_line()+
    scale_linetype_discrete(labels=function(b) hourly_model_labels[b])+
    theme_minimal()+
    labs(
      title=title,
      x="Forecast horizon (hours)",
      y=metric,
      linetype="Model"
    )
}
# label names match the model names; MSARIMA orders corrected; weekly SARIMA added
hourly_model_labels<-c(arima_2_0_4_notrend="arima(2,0,4)_notrend",arima_3_0_5_t="ARIMA(3,0,5)",arima_2_0_4_t="ARIMA(2,0,4)",sarima_2_0_2_101_t="SARIMA(2,0,2)(1,0,1)[24]",msarima_102_003_103="MSARIMA(1,0,2)[1](0,0,3)[6](1,0,3)[24]",ets_MAM="ETS(MAM)",sarima_weekly="Weekly SARIMA [168]")

# DM p-values per horizon, every model against the benchmark
# two-sided test (H0: equal accuracy; the RMSE tables show which model is better),
# varestimator = "bartlett" keeps the long-run variance positive, and the overlap of the
# forecast errors depends on the distance between origins: with one origin per day
# (origin_step = 24) and h <= 24 the forecasts do not overlap, so dm.test gets h = 1
dm.pvalues<-function(forecast_list,benchmark){
  out<-data.frame(horizon=seq_along(forecast_list))
  others<-setdiff(colnames(forecast_list[[1]][["errors"]]),benchmark)
  for(m in others){
    out[[paste0(m,"_vs_",benchmark)]]<-sapply(seq_along(forecast_list),function(h){
      e<-forecast_list[[h]][["errors"]]
      ok<-complete.cases(e[,c(m,benchmark)])
      tryCatch(
        dm.test(e[[m]][ok],e[[benchmark]][ok],h=max(1,ceiling(h/origin_step)),power=2,alternative="two.sided",varestimator="bartlett")$p.value,
        error=function(err) NA
      )
    })
  }
  out
}

# squared-error loss matrix for the MCS, only origins where all models have a forecast
loss.matrix<-function(errors){
  as.matrix(na.omit(errors))^2
}

check.coherence <- function(f) {
  
  zone1.diff <- f$zone1 -
    (f$area1 + f$area2 + f$area3 + f$area4)
  
  zone2.diff <- f$zone2 -
    (f$area5 + f$area6)
  
  NYC.diff <- f$NYC -
    (f$zone1 + f$zone2)
  
  data.frame(
    zone1_max_diff = max(abs(zone1.diff), na.rm = TRUE),
    zone2_max_diff = max(abs(zone2.diff), na.rm = TRUE),
    NYC_max_diff   = max(abs(NYC.diff), na.rm = TRUE)
  )
}

# observed demand equilibrium observations -----------------------------------------
#the observed quantities demanded  for the different areas correspond to 
#the market equilibrium, where quantity supplied equals quantity demanded.
#this implies that a change in observed quantity demand corresponds to
# a possible shift in underlying real demand or real supply.A more accurate
# forecast can be probably obtained when one takes this into account.
#Unfortunately we don't have any additional regressors available.


# data exploration:moments, seasonality, stationarity, weak dependence --------------------------------------------------------
english.weekday<-function(date){c("Sunday","Monday","Tuesday","Wednesday","Thursday","Friday","Saturday")[as.POSIXlt(date)$wday+1]}

#setup
subset1<-data1
nrows<-nrow(data1)
colnames(subset1)<-c("year","month","day","hour","area1","area2" ,"area3","area4","area5","area6")
subset1$date<-as.Date(paste(subset1$year,subset1$month,subset1$day, sep="-"))
subset1$weekday <- english.weekday(subset1$date)
table(subset1$weekday)   # check: every weekday should appear
subset1$timeofday <- ifelse(subset1$hour >= 6 & subset1$hour < 12, "morning",ifelse(subset1$hour >= 12 & subset1$hour < 18, "afternoon",ifelse(subset1$hour >= 18, "evening", "night")))
subset1$hour0 <- as.integer(subset1$hour == 0) 
subset1$hour1 <- as.integer(subset1$hour == 1) 
subset1$hour2 <- as.integer(subset1$hour == 2) 
subset1$hour3 <- as.integer(subset1$hour == 3) 
subset1$hour4 <- as.integer(subset1$hour == 4) 
subset1$hour5 <- as.integer(subset1$hour == 5) 
subset1$hour6 <- as.integer(subset1$hour == 6) 
subset1$hour7 <- as.integer(subset1$hour == 7) 
subset1$hour8 <- as.integer(subset1$hour == 8) 
subset1$hour9 <- as.integer(subset1$hour == 9) 
subset1$hour10 <- as.integer(subset1$hour == 10) 
subset1$hour11 <- as.integer(subset1$hour == 11) 
subset1$hour12 <- as.integer(subset1$hour == 12) 
subset1$hour13 <- as.integer(subset1$hour == 13) 
subset1$hour14 <- as.integer(subset1$hour == 14) 
subset1$hour15 <- as.integer(subset1$hour == 15) 
subset1$hour16 <- as.integer(subset1$hour == 16) 
subset1$hour17 <- as.integer(subset1$hour == 17) 
subset1$hour18 <- as.integer(subset1$hour == 18) 
subset1$hour19 <- as.integer(subset1$hour == 19) 
subset1$hour20 <- as.integer(subset1$hour == 20) 
subset1$hour21 <- as.integer(subset1$hour == 21) 
subset1$hour22 <- as.integer(subset1$hour == 22) 
subset1$hour23 <- as.integer(subset1$hour == 23)
subset1$Mondaymorning    <- as.integer(subset1$weekday == "Monday"    & subset1$timeofday == "morning")
subset1$Mondayafternoon  <- as.integer(subset1$weekday == "Monday"    & subset1$timeofday == "afternoon")
subset1$Mondayevening    <- as.integer(subset1$weekday == "Monday"    & subset1$timeofday == "evening")
subset1$Mondaynight      <- as.integer(subset1$weekday == "Monday"    & subset1$timeofday == "night")
subset1$Tuesdaymorning   <- as.integer(subset1$weekday == "Tuesday"   & subset1$timeofday == "morning")
subset1$Tuesdayafternoon <- as.integer(subset1$weekday == "Tuesday"   & subset1$timeofday == "afternoon")
subset1$Tuesdayevening   <- as.integer(subset1$weekday == "Tuesday"   & subset1$timeofday == "evening")
subset1$Tuesdaynight     <- as.integer(subset1$weekday == "Tuesday"   & subset1$timeofday == "night")
subset1$Wednesdaymorning <- as.integer(subset1$weekday == "Wednesday" & subset1$timeofday == "morning") 
subset1$Wednesdayafternoon <- as.integer(subset1$weekday == "Wednesday" & subset1$timeofday == "afternoon") 
subset1$Wednesdayevening <- as.integer(subset1$weekday == "Wednesday" & subset1$timeofday == "evening") 
subset1$Wednesdaynight <- as.integer(subset1$weekday == "Wednesday" & subset1$timeofday == "night") 
subset1$Thursdaymorning <- as.integer(subset1$weekday == "Thursday" & subset1$timeofday == "morning") 
subset1$Thursdayafternoon <- as.integer(subset1$weekday == "Thursday" & subset1$timeofday == "afternoon") 
subset1$Thursdayevening <- as.integer(subset1$weekday == "Thursday" & subset1$timeofday == "evening") 
subset1$Thursdaynight <- as.integer(subset1$weekday == "Thursday" & subset1$timeofday == "night") 
subset1$Fridaymorning <- as.integer(subset1$weekday == "Friday" & subset1$timeofday == "morning") 
subset1$Fridayafternoon <- as.integer(subset1$weekday == "Friday" & subset1$timeofday == "afternoon") 
subset1$Fridayevening <- as.integer(subset1$weekday == "Friday" & subset1$timeofday == "evening") 
subset1$Fridaynight <- as.integer(subset1$weekday == "Friday" & subset1$timeofday == "night") 
subset1$Saturdaymorning <- as.integer(subset1$weekday == "Saturday" & subset1$timeofday == "morning") 
subset1$Saturdayafternoon <- as.integer(subset1$weekday == "Saturday" & subset1$timeofday == "afternoon") 
subset1$Saturdayevening <- as.integer(subset1$weekday == "Saturday" & subset1$timeofday == "evening") 
subset1$Saturdaynight <- as.integer(subset1$weekday == "Saturday" & subset1$timeofday == "night") 
subset1$Sundaymorning <- as.integer(subset1$weekday == "Sunday" & subset1$timeofday == "morning") 
subset1$Sundayafternoon <- as.integer(subset1$weekday == "Sunday" & subset1$timeofday == "afternoon") 
subset1$Sundayevening <- as.integer(subset1$weekday == "Sunday" & subset1$timeofday == "evening") 
subset1$Sundaynight <- as.integer(subset1$weekday == "Sunday" & subset1$timeofday == "night")
subset1$april <- as.integer(subset1$month == 4) 
subset1$may <- as.integer(subset1$month == 5) 
subset1$june <- as.integer(subset1$month == 6) 
subset1$july <- as.integer(subset1$month == 7) 
subset1$august <- as.integer(subset1$month == 8)
subset1$trend<-1:nrow(subset1)
subset1$outlier1<-as.integer(subset1$date >= as.Date("2023-04-28") &subset1$date <= as.Date("2023-05-01"))
subset1$outlier2<-as.integer(subset1$date >= as.Date("2023-05-19") &subset1$date <= as.Date("2023-05-21"))
subset1$outlier3 <-as.integer(subset1$date >= as.Date("2023-07-15") &subset1$date <= as.Date("2023-07-17"))
data1<-subset1

#monthly patterns
subset1a<-subset(subset1, month== 4)
subset1b<-subset(subset1, month== 5)
subset1c<-subset(subset1, month== 6)
subset1d<-subset(subset1, month== 7)
subset1e<-subset(subset1, month== 8)


# full sample daily
daily1 <- aggregate(area1 ~ date, data = subset1, FUN = sum)
daily2 <- aggregate(area2 ~ date, data = subset1, FUN = sum)
daily3 <- aggregate(area3 ~ date, data = subset1, FUN = sum)
daily4 <- aggregate(area4 ~ date, data = subset1, FUN = sum)
daily5 <- aggregate(area5 ~ date, data = subset1, FUN = sum)
daily6 <- aggregate(area6 ~ date, data = subset1, FUN = sum)
daily<-data.frame(daily1,daily2,daily3,daily4,daily5,daily6)
daily$weekday <- english.weekday(daily$date)
daily$month <- as.integer(format(daily$date,"%m"))
daily$dayofmonth <- as.integer(format(daily$date,"%d"))
crisis_dates <- data.frame(start = as.Date(c("2023-04-28","2023-05-19", "2023-07-15")),end= as.Date(c("2023-05-01","2023-05-21" ,"2023-07-17")))

#plot of daily data full sample
            ggplot(daily, aes(x=date)) +
              geom_line(aes(y=area1),color="steelblue",linewidth=0.5)+
              geom_line(aes(y=area2),color="green",linewidth=0.5)+
              geom_line(aes(y=area3),color="orange",linewidth=0.5)+
              geom_line(aes(y=area4),color="purple",linewidth=0.5)+
              geom_line(aes(y=area5),color="tomato",linewidth=0.5)+
              geom_line(aes(y=area6),color="black",linewidth=0.5)+
              theme_minimal()+
              labs(title="")+
              ylab("demand")+
              theme(legend.position="bottom")+
              geom_rect(data = crisis_dates,aes(xmin = start,xmax = end,ymin = -Inf,ymax = Inf),inherit.aes = FALSE,fill="red",alpha=0.2)

#distribution of demand across time of day and weekday
                ggplot(subset1, aes(x=weekday, y=area1,fill=timeofday)) +
                  geom_bar(stat="identity")+
                  theme_minimal()+
                  theme(legend.position="bottom")+
                  labs(title ="Area 1 demand by weekday and time of day")+
                  xlab("days of the week")+
                  ylab("demand")
                ggplot(subset1, aes(x=weekday, y=area2,fill=timeofday)) +
                  geom_bar(stat="identity")+
                  theme_minimal()+
                  theme(legend.position="bottom")+
                  labs(title ="Area 2 demand by weekday and time of day")+
                  xlab("days of the week")+
                  ylab("demand")
                ggplot(subset1, aes(x=weekday, y=area3,fill=timeofday)) +
                  geom_bar(stat="identity")+
                  theme_minimal()+
                  theme(legend.position="bottom")+
                  labs(title ="Area 3 demand by weekday and time of day")+
                  xlab("days of the week")+
                  ylab("demand")
                ggplot(subset1, aes(x=weekday, y=area4,fill=timeofday)) +
                  geom_bar(stat="identity")+
                  theme_minimal()+
                  theme(legend.position="bottom")+
                  labs(title ="Area 4 demand by weekday and time of day")+
                  xlab("days of the week")+
                  ylab("demand")
                ggplot(subset1, aes(x=weekday, y=area5,fill=timeofday)) +
                  geom_bar(stat="identity")+
                  theme_minimal()+
                  theme(legend.position="bottom")+
                  labs(title ="Area 5 demand by weekday and time of day")+
                  xlab("days of the week")+
                  ylab("demand")
                ggplot(subset1, aes(x=weekday, y=area6,fill=timeofday)) +
                  geom_bar(stat="identity")+
                  theme_minimal()+
                  theme(legend.position="bottom")+
                  labs(title ="Area 6 demand by weekday and time of day")+
                  xlab("days of the week")+
                  ylab("demand")

# distribution of demand per area: does CLT apply?-> does I.I.D apply? 
                ggplot(subset1, aes(x=area1))+
                  geom_histogram(color="black", fill="steelblue")+
                  xlab("")+
                  ylab("")+
                  labs(title ="Distribution of Area 1 hourly demand")+
                  theme_minimal()
                ggplot(subset1, aes(x=area2))+
                  geom_histogram(color="black", fill="steelblue")+
                  xlab("")+
                  ylab("")+
                  labs(title ="Distribution of Area 2 hourly demand")+
                  theme_minimal()
                ggplot(subset1, aes(x=area3))+
                  geom_histogram(color="black", fill="steelblue")+
                  xlab("")+
                  ylab("")+
                  labs(title ="Distribution of Area 3 hourly demand")+
                  theme_minimal()
                ggplot(subset1, aes(x=area4))+
                  geom_histogram(color="black", fill="steelblue")+
                  xlab("")+
                  ylab("")+
                  labs(title ="Distribution of Area 4 hourly demand")+
                  theme_minimal()
                ggplot(subset1, aes(x=area5))+
                  geom_histogram(color="black", fill="steelblue")+
                  xlab("")+
                  ylab("")+
                  labs(title ="Distribution of Area 5 hourly demand")+
                  theme_minimal()
                ggplot(subset1, aes(x=area6))+
                  geom_histogram(color="black", fill="steelblue")+
                  xlab("")+
                  ylab("")+
                  labs(title ="Distribution of Area 6 hourly demand")+
                  theme_minimal()

# seasonality plot raw data: -> indication of seasonality (informal eyeballing)
              s1<-ts(subset1$area6, frequency = 24)
              seasonplot(s1, main = "Area 6 hourly seasonality")
              #conclusion: seems to be two seasonal patterns? -> not the seasonality of the data
              #(this follows from discussion first presentation, initially we disregarded the less strong pattern for noise)
              s2<-ts(subset1$area6, frequency = 6)
              seasonplot(s2, main = "Area 6 time of day seasonality")
              #conclusion: seems to be at least 4 seasonal patterns? -> not the seasonality of the data
              s3<-ts(subset1$area6, frequency = 24*7)
              seasonplot(s3, main = "Area 6 weekly seasonality")
              #conclusion: seems to be one seasonal pattern (although some noise)-> clear pattern for weekends (first two peaks) and weekdays (five later peaks)
              #=> we seem to have at least a weekly seasonal pattern in the data, (although it is possible there is also additional daily seasonality)
              

# seasonality plot daily data: -> no indication of seasonality (informal eyeballing)
              ggplot(daily, aes(x=dayofmonth, y=area6, colour=factor(month.name[month], levels=month.name))) +
                geom_line()+theme_minimal()+labs(title="Area 6 daily demand per month",x="day of the month",y="demand",colour="month")

# weak dependence?-> does autocorrelation die out quickly? (informal eyeballing)
              acf(daily$area1, main = "ACF of Area 1 daily demand")
              acf(daily$area2, main = "ACF of Area 2 daily demand")
              acf(daily$area3, main = "ACF of Area 3 daily demand")
              acf(daily$area4, main = "ACF of Area 4 daily demand")
              acf(daily$area5, main = "ACF of Area 5 daily demand")
              acf(daily$area6, main = "ACF of Area 6 daily demand")
              #conclusion: no weak dependence due to seasonality at least

# investigation seasonality for areas by regression-> deseasonalizing
# we perform a simplistic deseasonalizing strategy for expositional purpose only
                a1<-lm(area1~hour1+hour2+hour3+hour4+hour5+hour6+hour7+hour8+hour9+hour10+hour11+hour12+hour13+hour14+hour15+hour16+hour17+hour18+hour19+hour20+hour21+hour22+hour23+Tuesdaymorning+Tuesdayafternoon+Tuesdayevening+Tuesdaynight+Wednesdaymorning+Wednesdayafternoon+Wednesdayevening+Wednesdaynight+
                         Thursdaymorning+Thursdayafternoon+Thursdayevening+Thursdaynight+Fridaymorning+Fridayafternoon+Fridayevening+Fridaynight+Saturdaymorning+Saturdayafternoon+Saturdayevening+Saturdaynight+
                         Sundaymorning+Sundayafternoon+Sundayevening+Sundaynight+may+june+july+august,subset1)
                
                a2<-lm(area2~hour1+hour2+hour3+hour4+hour5+hour6+hour7+hour8+hour9+hour10+hour11+hour12+hour13+hour14+hour15+hour16+hour17+hour18+hour19+hour20+hour21+hour22+hour23+Tuesdaymorning+Tuesdayafternoon+Tuesdayevening+Tuesdaynight+Wednesdaymorning+Wednesdayafternoon+Wednesdayevening+Wednesdaynight+
                         Thursdaymorning+Thursdayafternoon+Thursdayevening+Thursdaynight+Fridaymorning+Fridayafternoon+Fridayevening+Fridaynight+Saturdaymorning+Saturdayafternoon+Saturdayevening+Saturdaynight+
                         Sundaymorning+Sundayafternoon+Sundayevening+Sundaynight+may+june+july+august,subset1)
                
                a3<-lm(area3~hour1+hour2+hour3+hour4+hour5+hour6+hour7+hour8+hour9+hour10+hour11+hour12+hour13+hour14+hour15+hour16+hour17+hour18+hour19+hour20+hour21+hour22+hour23+Tuesdaymorning+Tuesdayafternoon+Tuesdayevening+Tuesdaynight+Wednesdaymorning+Wednesdayafternoon+Wednesdayevening+Wednesdaynight+
                         Thursdaymorning+Thursdayafternoon+Thursdayevening+Thursdaynight+Fridaymorning+Fridayafternoon+Fridayevening+Fridaynight+Saturdaymorning+Saturdayafternoon+Saturdayevening+Saturdaynight+
                         Sundaymorning+Sundayafternoon+Sundayevening+Sundaynight+may+june+july+august,subset1)
                
                a4<-lm(area4~hour1+hour2+hour3+hour4+hour5+hour6+hour7+hour8+hour9+hour10+hour11+hour12+hour13+hour14+hour15+hour16+hour17+hour18+hour19+hour20+hour21+hour22+hour23+Tuesdaymorning+Tuesdayafternoon+Tuesdayevening+Tuesdaynight+Wednesdaymorning+Wednesdayafternoon+Wednesdayevening+Wednesdaynight+
                         Thursdaymorning+Thursdayafternoon+Thursdayevening+Thursdaynight+Fridaymorning+Fridayafternoon+Fridayevening+Fridaynight+Saturdaymorning+Saturdayafternoon+Saturdayevening+Saturdaynight+
                         Sundaymorning+Sundayafternoon+Sundayevening+Sundaynight+may+june+july+august,subset1)
                
                a5<-lm(area5~hour1+hour2+hour3+hour4+hour5+hour6+hour7+hour8+hour9+hour10+hour11+hour12+hour13+hour14+hour15+hour16+hour17+hour18+hour19+hour20+hour21+hour22+hour23+Tuesdaymorning+Tuesdayafternoon+Tuesdayevening+Tuesdaynight+Wednesdaymorning+Wednesdayafternoon+Wednesdayevening+Wednesdaynight+
                         Thursdaymorning+Thursdayafternoon+Thursdayevening+Thursdaynight+Fridaymorning+Fridayafternoon+Fridayevening+Fridaynight+Saturdaymorning+Saturdayafternoon+Saturdayevening+Saturdaynight+
                         Sundaymorning+Sundayafternoon+Sundayevening+Sundaynight+may+june+july+august,subset1)
                
                a6<-lm(area6~hour1+hour2+hour3+hour4+hour5+hour6+hour7+hour8+hour9+hour10+hour11+hour12+hour13+hour14+hour15+hour16+hour17+hour18+hour19+hour20+hour21+hour22+hour23+Tuesdaymorning+Tuesdayafternoon+Tuesdayevening+Tuesdaynight+Wednesdaymorning+Wednesdayafternoon+Wednesdayevening+Wednesdaynight+
                         Thursdaymorning+Thursdayafternoon+Thursdayevening+Thursdaynight+Fridaymorning+Fridayafternoon+Fridayevening+Fridaynight+Saturdaymorning+Saturdayafternoon+Saturdayevening+Saturdaynight+
                         Sundaymorning+Sundayafternoon+Sundayevening+Sundaynight+may+june+july+august,subset1)

                    summary(a1)# conclusion: significant seasonality
                    summary(a2)# conclusion: significant seasonality
                    summary(a3)# conclusion: significant seasonality
                    summary(a4)# conclusion: significant seasonality
                    summary(a5)# conclusion: significant seasonality
                    summary(a6)# conclusion: significant seasonality

deseasonalyarea1<-residuals(a1)
deseasonalyarea2<-residuals(a2)
deseasonalyarea3<-residuals(a3)
deseasonalyarea4<-residuals(a4)
deseasonalyarea5<-residuals(a5)
deseasonalyarea6<-residuals(a6)

                    acf(subset1$area1, main = "ACF of Area 1 hourly demand")
                    acf(deseasonalyarea1, main = "ACF of deseasonalized Area 1 hourly demand") # autocorrelation dies out-> weak dependence better
                    acf(subset1$area2, main = "ACF of Area 2 hourly demand")
                    acf(deseasonalyarea2, main = "ACF of deseasonalized Area 2 hourly demand") # autocorrelation dies out-> weak dependence better
                    acf(subset1$area3, main = "ACF of Area 3 hourly demand")
                    acf(deseasonalyarea3, main = "ACF of deseasonalized Area 3 hourly demand") # autocorrelation dies out-> weak dependence better
                    acf(subset1$area4, main = "ACF of Area 4 hourly demand")
                    acf(deseasonalyarea4, main = "ACF of deseasonalized Area 4 hourly demand") # autocorrelation dies out-> weak dependence better
                    acf(subset1$area5, main = "ACF of Area 5 hourly demand")
                    acf(deseasonalyarea5, main = "ACF of deseasonalized Area 5 hourly demand") # autocorrelation dies out-> weak dependence better
                    acf(subset1$area6, main = "ACF of Area 6 hourly demand")
                    acf(deseasonalyarea6, main = "ACF of deseasonalized Area 6 hourly demand") # autocorrelation dies out-> weak dependence better
                    
                    s4<-ts(deseasonalyarea6, frequency = 24)
                    seasonplot(s4, main = "Area 6 hourly deseasonalized") # seasonality seems taken care off
                    s5<-ts(deseasonalyarea6, frequency = 24*7)
                    seasonplot(s5, main = "Area 6 hourly deseasonalized") # remains some seasonality due to not correctly accounting for everything

#summary statistics for different areas:
                    summary1(subset1$area1)
                    summary1(subset1$area2)
                    summary1(subset1$area3)
                    summary1(subset1$area4)
                    summary1(subset1$area5)
                    summary1(subset1$area6) # clearly larger demand comming from here

#conclusion:
# non weak stationarity due to seasonality and probable trend
# potential influential outliers
# assumption 1 and 2 for weak stationarity no evidence of problem
# assumption 3 problem due tot seasonality
# weak dependency no problem after accounting seasonality
# weak stationarity no problem after accounting seasonality
# seasonality: weekly, hourly time of day
# data limitations to capture true demand-> equilibrium observations


# unit root test area 6 hourly data  ---------------------------------------------------------
# since we have seasonal data we should test for a unit-root based on deseasonalised data
# deseasonalizing before running test induces positive bias on finite sample: h0 less likely to be rejected
# Pantula principle: start at the highest order of integration considered (d = 2)
# step 1: H0 I(2) vs H1 I(1) -> ADF on the first difference
# step 2: only if step 1 rejects: H0 I(1) vs H1 I(0) -> ADF on the level
# bootstrapping is always good practise although the sample is fairly large
                    
#setup
# the unit root test uses the training sample only (April-July); August is test data.
#          The seasonality regression is re-estimated on the training sample (without the
#          august dummy, which is zero for all training observations).
training_rows<-which(subset1$month >= 4 & subset1$month <= 7)
a6_train<-lm(update(formula(a6), . ~ . - august), data = subset1[training_rows,])
deseasonalyarea6_train<-residuals(a6_train)

subset2a<-subset1
subset2a$deseasonalarea6<-NA
subset2a$deseasonalarea6[training_rows]<-deseasonalyarea6_train

# pantula algorithm
            boot_adf(data = diff(deseasonalyarea6_train),deterministics = "intercept")
            boot_adf(data = deseasonalyarea6_train,deterministics = "trend")
            # conclusion: reject that we have I(1) data (even with positive bias on ho), but do we have a trend?
# deterministic trend investigation
            summary(lm(area6~trend,subset1[training_rows,]))
            # conclusion: our data is trend stationary at 1% significance, although this might be caused by the large sample size
            #             we should look at the forecast performance of both models with and without trend

# unit root test area 6 daily data ----------------------------------------
#setup
subset2b<-daily
subset2b$weekday <- english.weekday(daily$date)
subset2b$Monday<- as.integer(subset2b$weekday == "Monday")
subset2b$Tuesday<-as.integer(subset2b$weekday == "Tuesday")
subset2b$Wednesday<-as.integer(subset2b$weekday == "Wednesday")
subset2b$Thursday<-as.integer(subset2b$weekday == "Thursday")
subset2b$Friday<-as.integer(subset2b$weekday == "Friday")
subset2b$Saturday<-as.integer(subset2b$weekday == "Saturday")
subset2b$Sunday<-as.integer(subset2b$weekday == "Sunday")
nrows2<-nrow(daily)
subset2b$trend <-1:nrows2

# training sample only (April-July)
daily_training_rows<-which(as.Date(subset2b$date) < as.Date("2023-08-01"))
a6d<-lm(area6~Tuesday+Wednesday+Thursday+Friday+Saturday+Sunday,subset2b[daily_training_rows,])
deseasonalyarea6daily<-residuals(a6d)
subset2b$deseasonalarea6daily<-NA
subset2b$deseasonalarea6daily[daily_training_rows]<-deseasonalyarea6daily

acf(daily$area6, main = "ACF of Area 6 daily demand")
acf(deseasonalyarea6daily, main = "ACF of deseasonalized Area 6 daily demand")

#we consider bootstrap because of small sample size (Pantula principle as above)
            boot_adf(data = diff(deseasonalyarea6daily),deterministics = "intercept")
            boot_adf(data = deseasonalyarea6daily,deterministics = "trend")
            #conclusion: reject Ho of I(1)
            summary(lm(area6~trend,subset2b[daily_training_rows,]))
            #no trend indication
#conclusion it is stationary due to twice rejecting the null


# ARMA choice hourly data -------------------------------------------------
#setup
subset3<-subset2a
subset3a<-subset3[subset3$month >= 4 & subset3$month <= 7, ]   # training sample: April - July

#based on ACF and PACF (with outlier)
          acf(subset3a$area6, main = "ACF of Area 6 training data")
          pacf(subset3a$area6, main = "PACF of Area 6 training data")
          #conclusion: impossible to determine

#based on information criteria (without considering impact outlier)
          auto.arima(subset3a$area6)
          #conclusion: arima (4,1,2) without considering the deterministic trend we should difference the data?-> contradiction unit root test
          #            informal indication that we should account for a deterministic trend
          
          autolaglengthAIC2(5,5,subset3a$area6,subset3a$trend)
          autolaglengthSC2(5,5,subset3a$area6,subset3a$trend)
          auto.arima(subset3a$area6,xreg = subset3a$trend)
          # conclusion: either ARMA(3,0,5) (own functions) or ARMA(2,0,4) (auto.arima)
          #             difference may come from the use of AICc instead of AIC or other manipulations used in auto.arima
          #             our functions recreates the selection procedure done in auto.arima without considering these manipulations
          
          auto.arima(subset3a$deseasonalarea6,xreg = subset3a$trend)
          # conclusion; incorrectly deseasonalising induces higher model complexity
          
#diagnostic residual autocorrelation (without considering impact outlier)
                        b<-arima(x = subset3a$area6,c(2,0,4))
                        b1<-arima(x = subset3a$area6,c(2,0,4),xreg = subset3a$trend)
                        b2<-arima(x = subset3a$area6,order = c(3,0,5),xreg = subset3a$trend,method = "ML")
                        
                        acf(residuals(b1), main = "ACF of ARIMA(2,0,4) residuals")
                        acf(residuals(b2), main = "ACF of ARIMA(3,0,5) residuals")
                        #conclusion: still significant autocorrelation?
                        #             we don't account for seasonality
                        plot(residuals(b1), main = "Residuals of ARIMA(2,0,4)")
                        plot(residuals(b2), main = "Residuals of ARIMA(3,0,5)")
                        #conclusion: still significant seasonality
                        #             we don't have weak stationarity                      
                        adf(residuals(b1),deterministics = "intercept")
                        adf(residuals(b2),deterministics = "intercept")
                        #conclusion: no unit root

#based on information criteria (controlling for the outliers)
#setup
exo<-cbind(outlier1=subset3a$outlier1,outlier2=subset3a$outlier2,outlier3=subset3a$outlier3,trend=subset3a$trend)
#based on information criteria
                        autolaglengthAIC2(5,5,subset3a$area6,exo)
                        autolaglengthSC2(5,5,subset3a$area6,exo)
                        auto.arima(subset3a$area6,xreg = exo)
                        # conclusion: outliers have no effect on previous choice


#diagnostic residual autocorrelation (controlling for the outliers)
                      b4<-arima(x = subset3a$area6,c(2,0,4),xreg = exo)
                      b5<-arima(x = subset3a$area6,order = c(3,0,5),xreg = exo,method = "ML")
                      
                      acf(residuals(b4), main = "ACF of ARIMA(2,0,4) residuals, outlier dummies included")
                      acf(residuals(b5), main = "ACF of ARIMA(3,0,5) residuals, outlier dummies included")
                      
                      plot(residuals(b4), main = "Residuals of ARIMA(2,0,4), outlier dummies included")
                      plot(residuals(b5), main = "Residuals of ARIMA(3,0,5), outlier dummies included")
  
                      adf(residuals(b4),deterministics = "intercept")
                      adf(residuals(b5),deterministics = "intercept")
                      #overall no difference in conclusion, so we don't care about outliers

# ARMA choice daily data --------------------------------------------------
subset3b<-subset2b
subset3b<- subset2b[as.Date(subset2b$date) < as.Date("2023-08-01"),]

#based on information criteria (with outlier)
            autolaglengthAIC2(5,5,subset3b$area6,subset3b$trend)
            autolaglengthSC2(5,5,subset3b$area6,subset3b$trend)
            auto.arima(subset3b$area6,xreg=subset3b$trend)
            # conclusion all are consistent ARMA(1,0,0)
# diagnostics
            b22<-arima(subset3b$area6,c(1,0,0),xreg=subset3b$trend)

            acf(residuals(b22), main = "ACF of daily ARIMA(1,0,0) residuals")
            # no indication of leftover serial correlation


# SARIMA model hourly data -------------------------------------------------
# SARIMA models allows us to model seasonality in the data by specifying an ARIMA structure on it in tandem with the non-seasonal part of the data.
# a SARIMA (1,0,0)[1,0,0][24] models the current hourly [24] observation as a function of the previous hour (1,0,0) and the previous hour one day before [1,0,0].
#setup
            subset4a<-subset3a
            c1a <- ts(subset4a$area6, frequency = 24) # freq seasonality > freq. of data = hourly (incorrectly)
            c1as <- ts(subset4a$area6, frequency = 24*7) # freq = seasonal relationship (correctly)
            c1ass<-ts(subset4a$area6, frequency = 6) # freq seasonality > freq. of data = time of day (incorrectly)
            c1t<-1:length(c1as)

# information critera
#we investigate the need of the linear trend since it was significant in the previous analysis but
#  it may simply try to capture some of the seasonality.
auto.arima(c1a)
auto.arima(c1as)
auto.arima(c1a,xreg=subset4a$trend)
auto.arima(c1as,xreg=subset4a$trend)
auto.arima(c1ass,xreg=subset4a$trend)
#observations:when considering the freq. of the data to be 24*7 and trend = T we get SARIMA (4,0,1) (case 1)
#             when considering the freq. of the data to be 24*7 and trend = F we get SARIMA (2,0,1)(0,1,0)[168] (case 2)
#             when considering the freq. of the data to be 24 and trend = T we get SARIMA  (2,0,2)[1,0,1][24] (case 3)
#             when considering the freq. of the data to be 24 and trend = F we get SARIMA  (4,0,0)(2,1,0)[24] (case 3)
#             when considering the freq. of the data to be 6 and trend = F we get SARIMA  (5,0,5)(2,0,1)[6] (case 4)
#            -> in the first case we consider the freq data = the actual seasonal pattern + account for trend  -> no seasonality detected 
#            -> in the second case we consider the freq data = the actual seasonal pattern  -> seasonality detected 
#            => we have seasonality of pattern length 24*7 (clear based on season plot)  for which we find it can be described by (0,1,0) 
#            => the trending behavior detected so far tried probably to capture the seasonal variation since SARIMA (4,0,1) is dominated by our ARMA models (previous section)
#            -> in the third and fourth case we consider the freq data < the actual seasonal pattern  -> seasonality detected (incorrectly)
# conclusion:   we go for SARIMA(2,0,1)(0,1,0)[168] interpretation: we should look at weekly differences
#             ! Since the [168] is big auto.arima might not be able to find the optimal model see section: "weekly seasonality"

e1<-auto.arima(c1as,xreg=subset4a$trend)   # stored as e1: used again in the weekly seasonality section
e1

# information critera 2
# auto.arima can only grasp one potential seasonal relationship. are there more? -> we test for hourly, time of day patterns
# weekly patterns are not included here, because auto.msarima with a lag of 168 next to
#          6 and 24 is very slow; the weekly seasonality is modelled separately in the section
#          "weekly seasonality" below (the season plot with frequency 168 shows a clear weekly pattern)
# including a deterministic trend causes perfect colinearity, hence is left out.
          auto.msarima(c1a, orders = list(ar = c(3, 3,3), i = c(0,0,0), ma = c(4,4,4)),lags = c(6,24), initial = c("backcasting"), ic = c("AICc"), loss = c("likelihood","MSE", "MAE", "HAM", "MSEh", "TMSE", "GTMSE", "MSCE"), h = 10,holdout = FALSE, bounds = c("usual", "admissible", "none"), silent = TRUE)
          # conclusion: disregarding the possibility of weekly relations, we have both hourly seasonality and time-of-day seasonality captured by SARIMA(1,0,2)[1](0,0,3)[6](1,0,3)[24]
          #              the lag-6 terms act as extra short-term dynamics within the day, not as a separate 6-hour season
# diagnostics
          c1<-arima(c1a,order = c(2, 0, 2),seasonal = list(order = c(1, 0, 1),period = 24),xreg=subset4a$trend)
          c2<-arima(c1ass,order = c(5, 0, 5),seasonal = list(order = c(2, 0, 1),period = 6),xreg=subset4a$trend)
          c3<-msarima(c1a, orders = list(ar = c(1, 0,1), i = c(0,0,0), ma = c(2,3,3)),lags = c(frequency(c1ass),frequency(c1a)),
                  constant = F,initial = c("backcasting"), ic = c("AICc"), loss = c("likelihood", "MSE", "MAE", "HAM", "MSEh", "TMSE","GTMSE", "MSCE"), h = 10, holdout = FALSE, bounds = c("usual","admissible", "none"), silent = TRUE)
          
          acf(residuals(c1), main = "ACF of hourly SARIMA residuals")
          plot(residuals(c1), main = "Residuals of hourly SARIMA model")
          adf(residuals(c1),deterministics = "intercept")
          #conclusion: improvement over the previous arma models
          acf(residuals(c2), main = "ACF of hourly SARIMA residuals")
          plot(residuals(c2), main = "Residuals of hourly SARIMA model")
          adf(residuals(c2),deterministics = "intercept")
          #conclusion: improvement over the previous arma models same as hourly
          acf(residuals(c3), main = "ACF of hourly SARIMA residuals")
          plot(residuals(c3), main = "Residuals of hourly SARIMA model")
          adf(residuals(c3),deterministics = "intercept")
          #conclusion: improvement over the previous arma models same as hourly

# brute force deseasonalizing
            c2a<-ts(subset4a$deseasonalarea6, frequency = 24)
            c2a<-auto.arima(c2a,xreg=subset4a$trend)
            # conclusion: not correctly deleting the seasonality induces higher model complexity
# diagnostics            
            acf(residuals(c2a), main = "ACF of deseasonalized SARIMA residuals")
            plot(residuals(c2a), main = "Residuals of deseasonalized SARIMA model")
            adf(residuals(c2a),deterministics = "intercept")
            # conclusion: no improvement


# SARIMA model daily data --------------------------------------------------
#setup
subset4b<-subset3b
#information criteria based selection
          c1b <- ts(subset4b$area6, frequency = 7)
          auto.arima(c1b,xreg=subset4b$trend)
          c11<-arima(subset4b$area6,c(1,0,0),xreg=subset4b$trend)
          # conclusion: we did not have any seasonality so no improvement over ARMA
#diagnostics
          acf(residuals(c11), main = "ACF of daily SARIMA residuals")
          plot(residuals(c11), main = "Residuals of daily SARIMA model")
          adf(residuals(c11))
          #conclusion: no improvements


# es model hourly data ----------------------------------------------------
# Exponentional smoothening models are state space models that allow for additive or multiplicative seasonality, trending and errors
# some linear exponential models are special cases of ARIMA models
# ets models can allow for dampening of the effect of a trend over time
# innovation errors are used as predictors of the model so we need to look at "forecast errors"

#setup
subset5a<-subset4a
d1a<-ts(subset4a$area6, frequency = 24)

# information criteria based selection
          d2<-ets(d1a)
          # conclusion: we have an MAM model implying multiplicative error and seasonality while having additive trend.

# diagnostics
          acf(residuals(d2), main = "ACF of ETS(MAM) residuals")
          plot(residuals(d2), main = "Residuals of ETS(MAM)")
          adf(residuals(d2),deterministics = "intercept")
          #conclusion: ACF seems worse than what we found for SARIMA model


# weekly seasonality: auto.arima and ETS ------------------------------------
# The season plot with frequency 168 shows a clear weekday/weekend pattern,
# so we try to model the weekly seasonality (168 hours) directly.
#setup
c1w<-c1as   # area 6 training data as a time series with frequency 168

# 1. auto.arima with period 168 (e1, estimated in the SARIMA section above)
          arimaorder(e1)
          # conclusion: only p, d, q are shown -> NO weekly part. Seasonal ARIMA with a period
          #             of 168 is too heavy to estimate in forecast::Arima, so auto.arima silently
          #             drops the seasonal candidates (Hyndman & Athanasopoulos: seasonal ARIMA
          #             is designed for short seasonal periods)

# 2. weekly SARIMA estimated with the smooth package (an MSARIMA with one seasonal lag
#    is a SARIMA, and its estimation is much faster). Orders chosen by AICc.
#    No normal differencing (unit root test), weekly difference allowed.
          e2<-auto.msarima(c1w,orders=list(ar=c(5,5),i=c(0,1),ma=c(5,5)),lags=c(1,168),initial="backcasting")
          print(e2)
          # store the selected orders for the out-of-sample forecast exercise (sarima_weekly).
          #          smooth::orders() shortens an all-zero vector to a single 0, so pad to 2 lags.
          pad2<-function(x) c(as.numeric(x),rep(0,2-length(x)))[1:2]
          weekly_orders<-list(ar=pad2(smooth::orders(e2)$ar),i=pad2(smooth::orders(e2)$i),ma=pad2(smooth::orders(e2)$ma))
          print(weekly_orders)
          # the model written in SARIMA notation
          cat(sprintf("SARIMA(%d,%d,%d)(%d,%d,%d)[168]\n",
                      weekly_orders$ar[1], weekly_orders$i[1], weekly_orders$ma[1],
                      weekly_orders$ar[2], weekly_orders$i[2], weekly_orders$ma[2]))
          acf(as.numeric(c1w)-as.numeric(fitted(e2)), main="ACF of weekly SARIMA residuals")
          # conclusion: alternative seasonal model: SARIMA(5,0,4)(4,1,3)[168]

# 3. ETS with period 168 (forecast package)
          e3<-tryCatch(ets(c1w),error=function(e){print(conditionMessage(e)); NULL})
          if (!is.null(e3)) print(e3$method)
          # conclusion: ets() cannot handle seasonal periods above 24: it gives a warning and
          #             drops the seasonality (it would need 167 initial seasonal states)

# 4. ETS with weekly seasonality, smooth package (es); model chosen by AICc.
#    backcasting avoids estimating the 168 initial seasonal states
          e4<-es(c1w,initial="backcasting")
          print(e4)
          acf(as.numeric(c1w)-as.numeric(fitted(e4)), main="ACF of weekly ETS residuals")

# 5. forced Auto-SARIMA with period 168 (forecast package, not smooth).
#    The free auto.arima (e1) drops the weekly part; here one weekly difference is forced
#    (D = 1) and auto.arima chooses the rest. d = 0: unit root test. max.P = max.Q = 1:
#    otherwise far too slow. No trend regressor: after the weekly difference the trend is
#    a constant, included as drift if AICc prefers it. trace = TRUE shows every model tried.
#    SLOW: 30 minutes to a few hours. Saved to disk, so it never has to be estimated twice.
          if (file.exists("e5_auto_sarima168.rds")) {
            e5<-readRDS("e5_auto_sarima168.rds")
          } else {
            e5<-auto.arima(c1w,d=0,D=1,max.P=1,max.Q=1,trace=TRUE)
            saveRDS(e5,"e5_auto_sarima168.rds")
          }
          arimaorder(e5)
          summary(e5)
          acf(residuals(e5),lag.max=400,main="ACF of forced Auto-SARIMA [168] residuals")
          o5<-arimaorder(e5)
          e5_name<-sprintf("Auto-SARIMA(%d,%d,%d)(%d,%d,%d)[168]%s (forecast pkg, D forced)",
                           o5["p"],o5["d"],o5["q"],o5["P"],o5["D"],o5["Q"],
                           if ("drift" %in% names(coef(e5))) " + drift" else "")
          # check of the start-up hours: RMSE without the first week
          in.sample.criteria(residuals(e5)[169:length(residuals(e5))])


# es model daily data -----------------------------------------------------
#setup
subset5b<-subset4b
d1b<-ts(subset5b$area6, frequency = 7)
# information criteria selection
        d22<-ets(d1b)
        #conclusion: a ANA model is the best
# diagnostics
        acf(residuals(d22), main = "ACF of daily ETS residuals")
        plot(residuals(d22), main = "Residuals of daily ETS")
        # conclusion: seems fine
        
# in sample selection hourly -----------------------------------------------------
# the residuals of all models except ARMA models seem to satisfy weak stationarity->consistency
# the first 100 residuals are dropped for every model (start-up period), so all models are compared on exactly the same hours

#setup
arima204_nt<-residuals(b)[100:length(residuals(b))]      
arima204_t<-residuals(b1)[100:length(residuals(b1))]
arima305_t<-residuals(b2)[100:length(residuals(b2))]
sarima202101_t<-residuals(c1)[100:length(residuals(c1))]
sarima505201_t<-residuals(c2)[100:length(residuals(c2))]
sarima202003103<-residuals(c3)[100:length(residuals(c3))]
ets2<-residuals(d2,type = "response")[100:length(residuals(d2))]

# in sample criteria: (=not accounting for model complexity)
            in.sample.criteria(arima204_nt)
            in.sample.criteria(arima204_t)
            in.sample.criteria(arima305_t)
            in.sample.criteria(sarima202101_t)
            in.sample.criteria(sarima505201_t)
            in.sample.criteria(sarima202003103)
            in.sample.criteria(ets2)
            #conclusion: SARIMA performs the best from the models in sample

#information criteria (=accounting for model complexity)
#all models are not directly comparable due different classes   
# given we are performing a forecasting exercise we put more weight on AIC then BIC in evaluation (although we are in sample)
            # AIC arima models
              AIC(b)     # ARIMA(2,0,4)_nt
              AIC(b1)    # ARIMA(2,0,4)_t
              AIC(b2)    # ARIMA(3,0,5)_t
              AIC(c1)    # SARIMA(2,0,2)(1,0,1)[24]_t
              AIC(c2)    # SARIMA(5,0,5)(2,0,1)[6]_t
              AIC(c3)    # MSARIMA(1,0,2)[1](0,0,3)[6](1,0,3)[24]              
            # AIC ets models               
              d2$aic     # ETS(MAM)
              # AIC is only comparable within the same model class and differencing:
              #          compare the stats::arima models (b, b1, b2, c1, c2) among each other; c3 (smooth)
              #          and d2 (ETS) use different likelihood set-ups. Across classes we compare RMSE/MAE
              #          in the in-sample comparison table below.
              
            # BIC arima
              BIC(b)    # ARIMA(2,0,4)_nt
              BIC(b1)   # ARIMA(2,0,4)_t
              BIC(b2)   # ARIMA(3,0,5)_t
              BIC(c1)    # SARIMA(2,0,2)(1,0,1)[24]_t
              BIC(c2)    # SARIMA(5,0,5)(2,0,1)[6]_t
              BIC(c3)    # MSARIMA(1,0,2)[1](0,0,3)[6](1,0,3)[24] 
            # BIC ets models 
              d2$bic    # ETS(MAM)
              # same as for AIC: only compare within a model class

# in-sample comparison table ----------------------------------------------
# one row per model: AIC/BIC, RMSE/MAE, number and mean of overshoots
# residuals = actual - fitted, first 99 hours dropped (same hours for every model)
# NOTE: AIC/BIC are only comparable within the same model class (stats::arima,
#       smooth, forecast::ets) and with the same differencing; compare all models on RMSE/MAE
y_train<-as.numeric(subset4a$area6)

table.row<-function(name,resid,aic,bic){
  r<-as.numeric(resid)
  r<-r[100:length(r)]
  r<-r[!is.na(r)]
  data.frame(
    Test=name,
    AIC=round(as.numeric(aic),1),
    BIC=round(as.numeric(bic),1),
    RMSE=round(sqrt(mean(r^2)),1),
    MAE=round(mean(abs(r)),1),
    n_overshoot=sum(r<0),
    mean_overshoot=round(if (any(r<0)) -mean(r[r<0]) else 0,1)
  )
}

rows<-list(
  table.row("ARIMA(2,0,4)",residuals(b),AIC(b),BIC(b)),
  table.row("ARIMA(2,0,4) + trend",residuals(b1),AIC(b1),BIC(b1)),
  table.row("ARIMA(3,0,5) + trend",residuals(b2),AIC(b2),BIC(b2)),
  table.row("SARIMA(2,0,2)(1,0,1)[24] + trend",residuals(c1),AIC(c1),BIC(c1)),
  table.row("SARIMA(5,0,5)(2,0,1)[6] + trend",residuals(c2),AIC(c2),BIC(c2)),
  table.row("MSARIMA(1,0,2)[1](0,0,3)[6](1,0,3)[24]",y_train-as.numeric(fitted(c3)),AIC(c3),BIC(c3)),
  table.row("ETS(MAM) [24]",residuals(d2,type="response"),d2$aic,d2$bic),
  table.row("auto.arima [168] (no weekly part)",residuals(e1),AIC(e1),BIC(e1)),
  table.row("Weekly SARIMA [168] (smooth)",y_train-as.numeric(fitted(e2)),AIC(e2),BIC(e2)),
  table.row("Weekly ETS [168] (smooth)",y_train-as.numeric(fitted(e4)),AIC(e4),BIC(e4)),
  table.row(e5_name,residuals(e5),AIC(e5),BIC(e5))
)
if (!is.null(e3)) {
  rows[[length(rows)+1]]<-table.row("ets() [168] (seasonality dropped)",
                                    residuals(e3,type="response"),e3$aic,e3$bic)
}

insample_table<-do.call(rbind,rows)
insample_table<-insample_table[order(insample_table$RMSE),]
print(insample_table,row.names=FALSE)

# same table in the layout of the slide: AIC/BIC and RMSE/MAE combined
slide_table<-data.frame(
  Test=insample_table$Test,
  `AIC/BIC`=paste(insample_table$AIC,"/",insample_table$BIC),
  `RMSE/MAE`=paste(insample_table$RMSE,"/",insample_table$MAE),
  `# overshoot`=insample_table$n_overshoot,
  `mean overshoot`=insample_table$mean_overshoot,
  check.names=FALSE
)
print(slide_table,row.names=FALSE)
# conclusion:
#RMSE, number and mean of overshoots ->SARIMA(5,0,4)(4,1,3)[168] (but don't account for model complexity)
#AIC/BIC ->SARIMA(5,0,4)(4,1,3)[168] (MSARIMA class)
#AIC/BIC -> the SARIMA(2,0,2)(1,0,1)[24] + trend (ARIMA class)
#MAE -> weekly ETS (but don't account for model complexity)
# => no clear results found but also not relevant for forecast exercise.
# => to limit running time we only use SARIMA(5,0,4)(4,1,3)[168] and not weekly ETS in out-of-sample comparison

# in sample selection daily -----------------------------------------------
burnin.daily<-20
arima100<-residuals(b22)[burnin.daily:length(residuals(b22))]
ets11<-residuals(d22,type = "response")[burnin.daily:length(residuals(d22))]

in.sample.criteria(arima100)
in.sample.criteria(ets11)

#conclusion: ets performs the best from the models

#information criteria 
# ARIMA class
AIC(b22)
BIC(b22)

# ETS class

d22$aic         # ETS(ANA)
d22$bic         # ETS(ANA)


# rolling window forecast daily, hourly and evaluation ----------------------------
# setup daily
rolling_daily_forecast <- window.forecast(subset2b$area6,subset2b$trend,freq. = 7,F.start = 123,H = 1,
                models = daily_models,
                window = "rolling",
                origin.step = 1)[["horizon_1"]]
#evaluation             
              in.sample.criteria(rolling_daily_forecast[["errors"]][["arima_1_0_0"]])
              in.sample.criteria(rolling_daily_forecast[["errors"]][["ets_ANA"]])
              #conclusion: based on different evaluation criterium different model is better.
              #            for absolute statistics the ets model is better
              #            for RMSE and mean fitted overshoot the arma model is better model is better
              
# setup hourly 
# test sample = all of August (as in the case):
# training sample: observations 1-2928 (April-July 2023)
# test sample: observations 2929-3672 (August 2023)
# every day at midnight (origin_step = 24) all parameters are re-estimated on the window
# and the models forecast the next day, h = 1,...,24: 31 origins, so 31 forecasts per horizon.
# (origin_step = 1 would make a new forecast every hour: about 24 times slower)
origin_step <- 24
hourly_forecast_data <- subset3[1:(2928 + 744), ]
rolling_hourly_forecasts <- window.forecast(hourly_forecast_data$area6,hourly_forecast_data$trend,freq. = 24,F.start = 2929,H = 24,
  models = hourly_models,
  window = "rolling",
  origin.step = origin_step)

#evaulation
                      rolling_hourly_criteria<-criteria.by.horizon(rolling_hourly_forecasts)
                      
                      # h = 1
                      print(rolling_hourly_criteria[rolling_hourly_criteria$horizon == 1, ])
                      #conclusion: RMSE/MAE, number and mean of overshoots SARIMA(5,0,4)(4,1,3)[168] (smooth) is best
                      
                      # per criterion over all horizons
                      metric.table(rolling_hourly_criteria,"RMSE")
                      #conclusion:  weekly SARIMA: "SARIMA(5,0,4)(4,1,3)[168] (smooth)" is the best
                      metric.table(rolling_hourly_criteria,"MAE")
                      #conclusion: weekly SARIMA: "SARIMA(5,0,4)(4,1,3)[168] (smooth)" is the best
                      metric.table(rolling_hourly_criteria,"total_overshoot")
                      #conclusion: weekly SARIMA: "SARIMA(5,0,4)(4,1,3)[168] (smooth)" is the best
                      metric.table(rolling_hourly_criteria,"mean_overshoot")
                      #conclusion: weekly SARIMA: "SARIMA(5,0,4)(4,1,3)[168] (smooth)" is the best
                      metric.table(rolling_hourly_criteria,"n_overshoot")
                      #conclusion: weekly SARIMA: "SARIMA(5,0,4)(4,1,3)[168] (smooth)" is the best
                      
                      plot.metric(rolling_hourly_criteria,"RMSE","Rolling-window RMSE by forecast horizon")
                      plot.metric(rolling_hourly_criteria,"total_overshoot","Rolling-window total overshoot by forecast horizon")
                      #why mean_overshoot,n_overshoot? -> taking inventory but no demand more expensive than missing demand

# expanding window forecast daily, hourly  and evaluation --------------------------
#setup daily
expanding_daily_forecast <- window.forecast(
  subset2b$area6,
  subset2b$trend,
  freq. = 7,
  F.start = 123,
  H = 1,
  models = daily_models,
  window = "expanding",
  origin.step = 1
)[["horizon_1"]]
#evaulation
                  in.sample.criteria(expanding_daily_forecast[["errors"]][["arima_1_0_0"]])
                  in.sample.criteria(expanding_daily_forecast[["errors"]][["ets_ANA"]])
                  #conclusion ets performs best for most criteria

expanding_hourly_forecasts <- window.forecast(
  hourly_forecast_data$area6,
  hourly_forecast_data$trend,
  freq. = 24,
  F.start = 2929,
  H = 24,
  models = hourly_models,
  window = "expanding",
  origin.step = origin_step
)

#evaluation
expanding_hourly_criteria<-criteria.by.horizon(expanding_hourly_forecasts)
expanding_hourly_criteria[expanding_hourly_criteria$horizon == 1, ]  
#conclusion: weekly SARIMA: "SARIMA(5,0,4)(4,1,3)[168] (smooth)" is the best

metric.table(expanding_hourly_criteria,"RMSE")
#conclusion: weekly SARIMA: "SARIMA(5,0,4)(4,1,3)[168] (smooth)" is the best
metric.table(expanding_hourly_criteria,"MAE")
#conclusion: weekly SARIMA: "SARIMA(5,0,4)(4,1,3)[168] (smooth)" is the best
metric.table(expanding_hourly_criteria,"total_overshoot")
#conclusion: weekly SARIMA: "SARIMA(5,0,4)(4,1,3)[168] (smooth)" is the best
metric.table(expanding_hourly_criteria,"mean_overshoot")
#conclusion: weekly SARIMA: "SARIMA(5,0,4)(4,1,3)[168] (smooth)" is the best
metric.table(expanding_hourly_criteria,"n_overshoot")
#conclusion: weekly SARIMA: "SARIMA(5,0,4)(4,1,3)[168] (smooth)" is the best

plot.metric(expanding_hourly_criteria,"RMSE","Expanding-window RMSE by forecast horizon")
plot.metric(expanding_hourly_criteria,"total_overshoot","Expanding-window total overshoot by forecast horizon")

# forecast model comparison hourly ----------------------------------------
# the diebold mariano test pairwise evaluates whether the forecast error are significantly different from each other
# DM tests for each forecast horizon (two-sided, H0: both models equally accurate)
# NA = dm.test could not be computed for that horizon
#set up hourly
rolling_hourly_dm<-dm.pvalues(rolling_hourly_forecasts,"msarima_102_003_103")
expanding_hourly_dm<-dm.pvalues(expanding_hourly_forecasts,"msarima_102_003_103")
trend_or_no_trend_R<-dm.pvalues(rolling_hourly_forecasts,"arima_2_0_4_notrend")
trend_or_no_trend_E<-dm.pvalues(expanding_hourly_forecasts,"arima_2_0_4_notrend")

#evaluation
                round(rolling_hourly_dm,4)
                round(expanding_hourly_dm,4)

                round(trend_or_no_trend_R,4)
                round(trend_or_no_trend_E,4)
                # accounting for a deterministic trend reduce forecast error significantly

# MCS for every forecast horizon
#The MCS procedure is designed to identify a set of models that are statistically consistent with being the best models.
# hence, it does not involve pairwise tesing as for the DM-test
#setup
rolling_hourly_mcs_results<-list()
expanding_hourly_mcs_results<-list()
for(h in 1:24){
                
                rolling_hourly_mcs_results[[h]]<-MCSprocedure(
                  Loss=loss.matrix(rolling_hourly_forecasts[[h]][["errors"]]),
                  alpha=0.05,
                  B=5000,
                  statistic="Tmax",
                  verbose=FALSE
                )
                
                expanding_hourly_mcs_results[[h]]<-MCSprocedure(
                  Loss=loss.matrix(expanding_hourly_forecasts[[h]][["errors"]]),
                  alpha=0.05,
                  B=5000,
                  statistic="Tmax",
                  verbose=FALSE
                )
              }
#evaluation
            # rolling window             
              for(h in 1:24){
                cat("\n==== Rolling window MCS, horizon", h, "====\n")
                print(rolling_hourly_mcs_results[[h]])
              }

             # expanding 
              for(h in 1:24){
                cat("\n==== Expanding window MCS, horizon", h, "====\n")
                print(expanding_hourly_mcs_results[[h]])
              }

# extra evaluation of the hourly forecasts --------------------------

# 1. average RMSE over all 24 horizons, per model
avg_rmse<-data.frame(
  rolling=round(colMeans(metric.table(rolling_hourly_criteria,"RMSE")[,-1]),1),
  expanding=round(colMeans(metric.table(expanding_hourly_criteria,"RMSE")[,-1]),1)
)
avg_rmse<-avg_rmse[order(avg_rmse$expanding),]
print(avg_rmse)

# 2. MCS: in how many of the 24 horizons is each model in the superior set?
#    (read from the printed output: the "Super Set of Models is composed by" line)
mcs_members<-function(x){
  out<-capture.output(print(x))
  line<-grep("Super Set of Models is composed by",out,value=TRUE)
  strsplit(trimws(sub(".*:\\s*","",line)),"\\s+")[[1]]
}
mcs_count<-function(res){
  models<-colnames(rolling_hourly_forecasts[[1]][["errors"]])
  table(factor(unlist(lapply(res,mcs_members)),levels=models))
}
mcs_counts<-data.frame(
  rolling=as.numeric(mcs_count(rolling_hourly_mcs_results)),
  expanding=as.numeric(mcs_count(expanding_hourly_mcs_results)),
  row.names=colnames(rolling_hourly_forecasts[[1]][["errors"]])
)
print(mcs_counts)

#conclusion: according to the MCS membership: weekly SARIMA: "SARIMA(5,0,4)(4,1,3)[168] (smooth)"  drops out the least of all models considered 
#            only for 04h to 05h the model is dropped out of the set for both rolling and expanding window
#            at later hours several models stay in the set.
#            for horizons further in the future the differences between the models are less clear

# 3. DM test for every PAIR of models at every horizon (two-sided, Bartlett variance,
#    h corrected for the daily origins, as in dm.pvalues)
#    NOTE: 21 pairs x 24 horizons = 504 tests per window: at 5% about 25 "significant"
#    results by chance, so the MCS remains the main joint test; the DM matrix is descriptive
pairwise.dm<-function(forecast_list){
  models<-colnames(forecast_list[[1]][["errors"]])
  pairs<-combn(models,2,simplify=FALSE)
  do.call(rbind,lapply(pairs,function(p){
    do.call(rbind,lapply(seq_along(forecast_list),function(h){
      e<-forecast_list[[h]][["errors"]]
      ok<-complete.cases(e[,p])
      pv<-tryCatch(
        dm.test(e[[p[1]]][ok],e[[p[2]]][ok],h=max(1,ceiling(h/origin_step)),
                power=2,alternative="two.sided",varestimator="bartlett")$p.value,
        error=function(err) NA)
      data.frame(horizon=h,model1=p[1],model2=p[2],p_value=round(pv,4),
                 lower_MSE=ifelse(mean(e[[p[1]]][ok]^2)<mean(e[[p[2]]][ok]^2),p[1],p[2]))
    }))
  }))
}

# win matrix: cell [row, column] = number of horizons (out of 24) at which the row
# model is SIGNIFICANTLY more accurate than the column model (p < 0.05)
win.matrix<-function(dm_long,alpha=0.05){
  models<-unique(c(dm_long$model1,dm_long$model2))
  M<-matrix(0,length(models),length(models),dimnames=list(better=models,worse=models))
  sig<-dm_long[!is.na(dm_long$p_value) & dm_long$p_value<alpha,]
  for(k in seq_len(nrow(sig))){
    w<-sig$lower_MSE[k]
    l<-if (w==sig$model1[k]) sig$model2[k] else sig$model1[k]
    M[w,l]<-M[w,l]+1
  }
  M
}

dm_all_R<-pairwise.dm(rolling_hourly_forecasts)
dm_all_E<-pairwise.dm(expanding_hourly_forecasts)
win.matrix(dm_all_R)
win.matrix(dm_all_E)
#conclusion: according to the pairwise DM-test weekly SARIMA: "SARIMA(5,0,4)(4,1,3)[168] (smooth)" is the best 
#            considering expanding and rolling window forecast

# 4. full results table: in-sample AIC/BIC + out-of-sample statistics pooled over all
#    horizons (31 days x 24 hours = 744 forecasts per model), rolling and expanding
pooled.stats<-function(forecast_list){
  models<-colnames(forecast_list[[1]][["errors"]])
  do.call(rbind,lapply(models,function(m){
    e<-unlist(lapply(forecast_list,function(x) x[["errors"]][[m]]))
    e<-e[!is.na(e)]
    data.frame(
      model=m,
      RMSE=round(sqrt(mean(e^2)),1),
      MAE=round(mean(abs(e)),1),
      n_overshoot=sum(e<0),
      mean_overshoot=round(if (any(e<0)) -mean(e[e<0]) else 0,1),
      n_forecasts=length(e)
    )
  }))
}

# in-sample fits of the same models (training sample). NOTE: AIC/BIC only comparable
# within a model class
insample_fits<-list(
  arima_2_0_4_notrend=b,
  arima_3_0_5_t=b2,
  arima_2_0_4_t=b1,
  sarima_2_0_2_101_t=c1,
  msarima_102_003_103=c3,
  ets_MAM=d2,
  sarima_weekly=e2
)
ic<-data.frame(
  model=names(insample_fits),
  AIC=sapply(insample_fits,function(f) round(if (inherits(f,"ets")) f$aic else AIC(f),1)),
  BIC=sapply(insample_fits,function(f) round(if (inherits(f,"ets")) f$bic else BIC(f),1))
)

oos_rolling<-pooled.stats(rolling_hourly_forecasts)
oos_expanding<-pooled.stats(expanding_hourly_forecasts)
names(oos_rolling)[-1]<-paste0(names(oos_rolling)[-1],"_rolling")
names(oos_expanding)[-1]<-paste0(names(oos_expanding)[-1],"_expanding")
full_table<-merge(merge(ic,oos_rolling,by="model"),oos_expanding,by="model")
full_table$model<-hourly_model_labels[full_table$model]

# 5. the free auto.arima with period 168 (e1 = ARIMA(4,0,1) + trend, no weekly part)
#    added to the table: orders fixed (it has no weekly part to re-select anyway),
#    parameters re-estimated every day. Table only, not in the DM tests and MCS.
e1_model<-list(auto_arima_168=function(y,xr,newxr,steps){
  o<-arimaorder(e1)
  fit<-Arima(as.numeric(y),order=as.numeric(o[c("p","d","q")]),xreg=cbind(trend=xr),
             include.mean="intercept" %in% names(coef(e1)))
  as.numeric(forecast(fit,xreg=cbind(trend=newxr))$mean)
})
e1_rolling<-window.forecast(hourly_forecast_data$area6,hourly_forecast_data$trend,freq.=24,F.start=2929,H=24,
                            models=e1_model,window="rolling",origin.step=origin_step)
e1_expanding<-window.forecast(hourly_forecast_data$area6,hourly_forecast_data$trend,freq.=24,F.start=2929,H=24,
                              models=e1_model,window="expanding",origin.step=origin_step)
oos_e1_r<-pooled.stats(e1_rolling)
oos_e1_e<-pooled.stats(e1_expanding)
names(oos_e1_r)[-1]<-paste0(names(oos_e1_r)[-1],"_rolling")
names(oos_e1_e)[-1]<-paste0(names(oos_e1_e)[-1],"_expanding")
e1_row<-merge(merge(data.frame(model="auto_arima_168",AIC=round(AIC(e1),1),BIC=round(BIC(e1),1)),
                    oos_e1_r,by="model"),oos_e1_e,by="model")
e1_row$model<-"auto.arima [168]: ARIMA(4,0,1) + trend (no weekly part)"
full_table<-rbind(full_table,e1_row)

# 5b. the forced Auto-SARIMA with period 168 (e5, forecast package, weekly difference
#     forced) added to the table: orders fixed, parameters re-estimated every day.
#     Table only, not in the DM tests and MCS.
#     Its trace showed that every model with a weekly AR or MA term failed (Inf), so e5
#     only has the weekly DIFFERENCE: ARIMA(p,0,q)(0,1,0)[168]. That model is the same as
#     an ARMA(p,q) on the weekly differences w_t = y_t - y_(t-168) (plus a mean if e5 has
#     drift), and a forecast for hour t+h (h <= 24) is y_(t+h-168) + forecast of w_(t+h).
#     Estimating it that way avoids the very slow seasonal estimation with period 168.
#     (If e5 did have weekly AR/MA terms, the slow seasonal Arima() is used instead.)
e5_model<-list(auto_sarima_168_forced=function(y,xr,newxr,steps){
  o<-arimaorder(e5)
  y<-as.numeric(y)
  n<-length(y)
  has_drift<-"drift" %in% names(coef(e5))
  if (o["P"]==0 && o["Q"]==0) {
    w<-y[169:n]-y[1:(n-168)]                                   # weekly differences
    fit<-Arima(w,order=as.numeric(o[c("p","d","q")]),include.mean=has_drift)
    y[n-168+(1:steps)]+as.numeric(forecast(fit,h=steps)$mean)  # back to levels
  } else {
    fit<-Arima(ts(y,frequency=168),order=as.numeric(o[c("p","d","q")]),
               seasonal=list(order=as.numeric(o[c("P","D","Q")]),period=168),
               include.drift=has_drift)
    as.numeric(forecast(fit,h=steps)$mean)
  }
})
e5_rolling<-window.forecast(hourly_forecast_data$area6,hourly_forecast_data$trend,freq.=24,F.start=2929,H=24,
                            models=e5_model,window="rolling",origin.step=origin_step)
e5_expanding<-window.forecast(hourly_forecast_data$area6,hourly_forecast_data$trend,freq.=24,F.start=2929,H=24,
                              models=e5_model,window="expanding",origin.step=origin_step)
oos_e5_r<-pooled.stats(e5_rolling)
oos_e5_e<-pooled.stats(e5_expanding)
names(oos_e5_r)[-1]<-paste0(names(oos_e5_r)[-1],"_rolling")
names(oos_e5_e)[-1]<-paste0(names(oos_e5_e)[-1],"_expanding")
e5_row<-merge(merge(data.frame(model="auto_sarima_168_forced",AIC=round(AIC(e5),1),BIC=round(BIC(e5),1)),
                    oos_e5_r,by="model"),oos_e5_e,by="model")
e5_row$model<-e5_name
full_table<-rbind(full_table,e5_row)

full_table<-full_table[order(full_table$RMSE_expanding),]
print(full_table,row.names=FALSE)

# same in the layout of the slide, one table per window
slide.layout<-function(tab,window){
  data.frame(
    Test=tab$model,
    `AIC/BIC`=paste(tab$AIC,"/",tab$BIC),
    `RMSE/MAE`=paste(tab[[paste0("RMSE_",window)]],"/",tab[[paste0("MAE_",window)]]),
    `# overshoot`=tab[[paste0("n_overshoot_",window)]],
    `mean overshoot`=tab[[paste0("mean_overshoot_",window)]],
    check.names=FALSE
  )
}
print(slide.layout(full_table,"rolling"),row.names=FALSE)
print(slide.layout(full_table,"expanding"),row.names=FALSE)

# 6. export: one Excel file with all tables, MCS printouts as text, RMSE plots as images
if (!requireNamespace("writexl",quietly=TRUE)) install.packages("writexl")
dir.create("results",showWarnings=FALSE)
as.sheet<-function(M) data.frame(model=rownames(M),M,check.names=FALSE)
writexl::write_xlsx(
  list(
    insample           = slide_table,
    oos_rolling        = slide.layout(full_table,"rolling"),
    oos_expanding      = slide.layout(full_table,"expanding"),
    oos_all_statistics = full_table,
    avg_rmse           = as.sheet(avg_rmse),
    mcs_counts         = as.sheet(mcs_counts),
    rmse_rolling       = metric.table(rolling_hourly_criteria,"RMSE"),
    rmse_expanding     = metric.table(expanding_hourly_criteria,"RMSE"),
    dm_vs_msarima_R    = round(rolling_hourly_dm,4),
    dm_vs_msarima_E    = round(expanding_hourly_dm,4),
    dm_all_pairs_R     = dm_all_R,
    dm_all_pairs_E     = dm_all_E,
    dm_win_matrix_R    = as.sheet(win.matrix(dm_all_R)),
    dm_win_matrix_E    = as.sheet(win.matrix(dm_all_E))
  ),
  "results/citibike_results.xlsx"
)
writeLines(capture.output(for (h in 1:24) {
  cat("\n==== Rolling window MCS, horizon",h,"====\n"); print(rolling_hourly_mcs_results[[h]])
}),"results/mcs_rolling.txt")
writeLines(capture.output(for (h in 1:24) {
  cat("\n==== Expanding window MCS, horizon",h,"====\n"); print(expanding_hourly_mcs_results[[h]])
}),"results/mcs_expanding.txt")
ggsave("results/rmse_rolling.png",
       plot.metric(rolling_hourly_criteria,"RMSE","Rolling-window RMSE by forecast horizon"),width=9,height=5)
ggsave("results/rmse_expanding.png",
       plot.metric(expanding_hourly_criteria,"RMSE","Expanding-window RMSE by forecast horizon"),width=9,height=5)

# forecast comparison daily -----------------------------------------------
# two-sided (H0: both models equally accurate)

# moving window DM test
                  dm.test(
                    rolling_daily_forecast[["errors"]][["ets_ANA"]],
                    rolling_daily_forecast[["errors"]][["arima_1_0_0"]],
                    h = 1,
                    power = 2,
                    alternative = "two.sided"
                  )
                  # conclusion: no significant difference between models

# expanding window DM test
                dm.test(
                  expanding_daily_forecast[["errors"]][["ets_ANA"]],
                  expanding_daily_forecast[["errors"]][["arima_1_0_0"]],
                  h = 1,
                  power = 2,
                  alternative = "two.sided"
                )
                # conclusion: no significant difference between models

# moving window MCS
rolling_daily_mcs_result <- MCSprocedure(
  Loss = loss.matrix(rolling_daily_forecast[["errors"]]),
  alpha = 0.05,
  B = 5000,
  statistic = "Tmax"
)

rolling_daily_mcs_result@Info
print(rolling_daily_mcs_result)

# expanding window MCS
expanding_daily_mcs_result <- MCSprocedure(
  Loss = loss.matrix(expanding_daily_forecast[["errors"]]),
  alpha = 0.05,
  B = 5000,
  statistic = "Tmax"
)

expanding_daily_mcs_result@Info
print(expanding_daily_mcs_result)

# Time series part of cross-sectional Hierarchy -----------------------------------------
# two base models for all 9 series, written as functions:
#   - "sarima_weekly": weekly SARIMA (smooth, lags 1 and 168 + drift)
#   - "sarima_daily" : daily SARIMA (period 24 + linear trend)
# For every series the orders are selected ONCE on the training sample (April-July), with
# the same searches as for area 6 (auto.msarima / auto.arima), and then kept fixed; every
# day all parameters are re-estimated. The training fits also give the WLS residuals.
# Same set-up as for area 6: all of August, origins every hierarchy_origin_step hours,
# rolling and expanding window, h = 1,...,24.
# The steps per model and window are the same as before: base forecasts -> coherence check
# -> reconciliation (bottom-up, OLS, WLS) -> coherence check -> accuracy.

# checkpoint: save everything so far (area 6 results), restore with load("after_area6.RData")
save.image("after_area6.RData")

# which base models to run (run one at a time by keeping only one name here)
hierarchy_base_models <- c("sarima_weekly", "sarima_daily")

# distance between forecast origins in the hierarchy: the same as for area 6
# (origin_step = 24: one origin per day at midnight)
hierarchy_origin_step <- origin_step
n_hierarchy_origins <- length(seq(2929, 3672, by = hierarchy_origin_step))

#setup
hourly_forecast_data <- subset1[1:(2928 + 744), ]   # all of August
hourly_forecast_data$zone1<- hourly_forecast_data$area1+hourly_forecast_data$area2+hourly_forecast_data$area3+hourly_forecast_data$area4
hourly_forecast_data$zone2<- hourly_forecast_data$area6+hourly_forecast_data$area5
hourly_forecast_data$NYC<- hourly_forecast_data$zone2+hourly_forecast_data$zone1

bs <- paste0("area", 1:6)
us <- c("zone1", "zone2", "NYC")
all_series <- c(bs, us)

# timing check: one fit of each base model (fixed area 6 orders) on the area 6 training
# data, times the number of fits in the hierarchy (9 series x origins x 2 windows)
for (m in hierarchy_base_models) {
  area6_model <- if (m == "sarima_daily") "sarima_2_0_2_101_t" else m   # area 6 name of the model
  t_one <- system.time(
    hourly_models[[area6_model]](ts(hourly_forecast_data$area6[1:2928], frequency = 24),
                       hourly_forecast_data$trend[1:2928], 2929:2952, 24)
  )["elapsed"]
  cat(m, ": one fit", round(t_one), "s -> hierarchy about",
      round(t_one * 9 * n_hierarchy_origins * 2 / 3600, 1), "hours (+ the order selection per series)\n")
}

# order selection per series on the training sample (once) + training residuals for WLS
select_hierarchy_model <- function(model_name, s) {
  cat("order selection", model_name, "series", s, "\n")
  y  <- hourly_forecast_data[[s]][1:2928]
  tr <- hourly_forecast_data$trend[1:2928]
  if (model_name == "sarima_weekly") {
    fit <- smooth::auto.msarima(ts(y, frequency = 168),
                                orders = list(ar = c(5, 5), i = c(0, 1), ma = c(5, 5)),
                                lags = c(1, 168), initial = "backcasting")
    pad2 <- function(x) c(as.numeric(x), rep(0, 2 - length(x)))[1:2]
    ord <- smooth::orders(fit)
    list(orders = list(ar = pad2(ord$ar), i = pad2(ord$i), ma = pad2(ord$ma)),
         residuals = y - as.numeric(fitted(fit)))
  } else {
    fit <- auto.arima(ts(y, frequency = 24), xreg = cbind(trend = tr))
    list(template = fit, residuals = as.numeric(residuals(fit)))
  }
}

# forecast function of one series with its FIXED orders (parameters re-estimated each day)
hierarchy_forecast_function <- function(model_name, selection) {
  force(selection)
  if (model_name == "sarima_weekly") {
    function(y, xr, newxr, steps) {
      fit <- smooth::msarima(ts(as.numeric(y), frequency = 168),
                             orders = selection$orders, lags = c(1, 168),
                             constant = TRUE, initial = "backcasting")
      as.numeric(forecast(fit, h = steps)$mean)
    }
  } else {
    function(y, xr, newxr, steps) {
      template <- selection$template
      o <- arimaorder(template)
      seas <- if (all(c("P", "D", "Q") %in% names(o))) as.numeric(o[c("P", "D", "Q")]) else c(0, 0, 0)
      fit <- Arima(ts(as.numeric(y), frequency = 24),
                   order = as.numeric(o[c("p", "d", "q")]),
                   seasonal = list(order = seas, period = 24),
                   xreg = cbind(trend = xr),
                   include.mean = "intercept" %in% names(coef(template)),
                   include.drift = "drift" %in% names(coef(template)))
      as.numeric(forecast(fit, xreg = cbind(trend = newxr))$mean)
    }
  }
}

# forecasts for all 9 series, one base model, one window
# (selections = the per-series order selections of that base model)
forecast_hierarchy <- function(model_name, window_type, selections) {
  out <- lapply(all_series, function(s) {
    cat("\n", model_name, "-", window_type, "- series", s, "\n")
    window.forecast(hourly_forecast_data[[s]], hourly_forecast_data$trend,
                    freq. = 24, F.start = 2929, H = 24,
                    models = list(base = hierarchy_forecast_function(model_name, selections[[s]])),
                    window = window_type,
                    origin.step = hierarchy_origin_step)
  })
  names(out) <- all_series
  out
}

# coherence of the base forecasts (as before: maximum deviation per horizon)
base_coherence <- function(fc) {
  coherence_results <- do.call(rbind, lapply(1:24, function(h) {
    f <- data.frame(lapply(fc, function(x) x[[h]]$forecasts[, 1]))
    check.coherence(f)
  }))
  coherence_results$horizon <- 1:24
  coherence_results
}


# Forecast Reconciliation -------------------------------------------------
#setup 
hierarchy_series <- c(bs, us)
agg_mat1 <- matrix(
  c(
    1, 1, 1, 1, 0, 0,
    0, 0, 0, 0, 1, 1,
    1, 1, 1, 1, 1, 1),nrow = 3,byrow = TRUE)
rownames(agg_mat1) <- us
colnames(agg_mat1) <- bs

cons_mat1 <- matrix(0,nrow = 3,ncol = length(hierarchy_series))
colnames(cons_mat1) <- hierarchy_series
cons_mat1[1, "zone1"] <- 1
cons_mat1[1, c("area1", "area2", "area3", "area4")] <- -1
cons_mat1[2, "zone2"] <- 1
cons_mat1[2, c("area5", "area6")] <- -1
cons_mat1[3, "NYC"] <- 1
cons_mat1[3, bs] <- -1

cons_mat_ols <- matrix(0,nrow = 3,ncol = 9)
colnames(cons_mat_ols) <- c(
  "zone1", "zone2", "NYC",
  "area1", "area2", "area3", "area4", "area5", "area6")
cons_mat_ols[1, "zone1"] <- 1
cons_mat_ols[1, c("area1", "area2", "area3", "area4")] <- -1
cons_mat_ols[2, "zone2"] <- 1
cons_mat_ols[2, c("area5", "area6")] <- -1
cons_mat_ols[3, "NYC"] <- 1
cons_mat_ols[3, bs] <- -1
hierarchy_series <- c("area1", "area2", "area3","area4", "area5", "area6","zone1", "zone2", "NYC")

#base for forec
get_base_forecasts <- function(h,area_forecasts,zone1_forecast,zone2_forecast,nyc_forecast) {
  
  base <- cbind(
    area1 = area_forecasts[[1]][[h]]$forecasts[, 1],
    area2 = area_forecasts[[2]][[h]]$forecasts[, 1],
    area3 = area_forecasts[[3]][[h]]$forecasts[, 1],
    area4 = area_forecasts[[4]][[h]]$forecasts[, 1],
    area5 = area_forecasts[[5]][[h]]$forecasts[, 1],
    area6 = area_forecasts[[6]][[h]]$forecasts[, 1],
    zone1 = zone1_forecast[[h]]$forecasts[, 1],
    zone2 = zone2_forecast[[h]]$forecasts[, 1],
    NYC   = nyc_forecast[[h]]$forecasts[, 1]
  )
  
  return(base)
}
get_bottom_forecasts <- function(h, area_forecasts) {
  
  base_bottom <- cbind(
    area1 = area_forecasts[[1]][[h]]$forecasts[, 1],
    area2 = area_forecasts[[2]][[h]]$forecasts[, 1],
    area3 = area_forecasts[[3]][[h]]$forecasts[, 1],
    area4 = area_forecasts[[4]][[h]]$forecasts[, 1],
    area5 = area_forecasts[[5]][[h]]$forecasts[, 1],
    area6 = area_forecasts[[6]][[h]]$forecasts[, 1]
  )
  
  return(base_bottom)
}

# coherence check
check_coherence <- function(f) {
  
  # f is a forecast matrix
  # columns must be in the order:
  # area1,...,area6,zone1,zone2,NYC
  
  violations <- f %*% t(cons_mat1)
  
  max_violation <- max(abs(violations), na.rm = TRUE)
  
  return(max_violation)
}

# WLS residuals of one base model: from the per-series training fits (order selection).
# Column order as in the OLS/WLS base matrix: zone1, zone2, NYC, area1,...,area6
training_residuals <- function(selections) {
  res <- sapply(c(us, bs), function(s) selections[[s]]$residuals)
  res[complete.cases(res), , drop = FALSE]
}

# reconciliation of one set of base forecasts (one model, one window).
# Origins where a base forecast is missing (model failed on that day) are left out for
# ALL methods, so base and reconciled forecasts are evaluated on the same origins.
reconcile_hierarchy <- function(fc, res) {
  
  area_forecasts <- fc[bs]
  
  base <- lapply(1:24, function(h) get_base_forecasts(
    h = h, area_forecasts = area_forecasts,
    zone1_forecast = fc$zone1, zone2_forecast = fc$zone2, nyc_forecast = fc$NYC))
  
  empty <- function(n) matrix(NA_real_, n, 9, dimnames = list(NULL, hierarchy_series))
  base_eval <- list(); bu <- list(); ols <- list(); wls <- list()
  
  for (h in 1:24) {
    b  <- base[[h]]
    ok <- complete.cases(b)
    base_eval[[h]] <- b; base_eval[[h]][!ok, ] <- NA
    bu[[h]] <- empty(nrow(b)); ols[[h]] <- empty(nrow(b)); wls[[h]] <- empty(nrow(b))
    if (any(ok)) {
      b_ok   <- b[ok, , drop = FALSE]
      b_ols  <- b_ok[, c("zone1", "zone2", "NYC", bs), drop = FALSE]
      bottom <- b_ok[, bs, drop = FALSE]
      bu[[h]][ok, ]  <- as.matrix(csbu(base = bottom, agg_mat = agg_mat1))[, hierarchy_series, drop = FALSE]
      ols[[h]][ok, ] <- as.matrix(csrec(base = b_ols, agg_mat = agg_mat1, cons_mat = cons_mat_ols,
                                        comb = "ols"))[, hierarchy_series, drop = FALSE]
      wls[[h]][ok, ] <- as.matrix(csrec(base = b_ols, agg_mat = agg_mat1, cons_mat = cons_mat_ols,
                                        comb = "wls", res = res))[, hierarchy_series, drop = FALSE]
    }
  }
  names(base_eval) <- names(bu) <- names(ols) <- names(wls) <- paste0("horizon_", 1:24)
  
  coherence <- rbind(Base = sapply(base_eval, check_coherence),
                     Bottom_Up = sapply(bu, check_coherence),
                     OLS = sapply(ols, check_coherence),
                     WLS = sapply(wls, check_coherence))
  colnames(coherence) <- paste0("horizon_", 1:24)
  
  list(base = base_eval, bu = bu, ols = ols, wls = wls, coherence = coherence)
}

calculate_accuracy <- function(forecasts, actuals) {
  actuals <- actuals[, colnames(forecasts), drop = FALSE]
  
  residuals <- actuals - forecasts
  
  accuracy <- t(apply(residuals, 2, forecast.criteria))
  
  return(as.data.frame(accuracy))
}

# accuracy of base and reconciled forecasts (one model, one window)
hierarchy_accuracy <- function(rec, fc) {
  
  methods <- list(Base = rec$base, `Bottom-Up` = rec$bu, OLS = rec$ols, WLS = rec$wls)
  
  accuracy_comparison <- do.call(rbind, lapply(1:24, function(h) {
    target <- fc$area1[[h]]$target
    actuals_h <- hourly_forecast_data[target, hierarchy_series, drop = FALSE]
    do.call(rbind, lapply(names(methods), function(m) {
      a <- calculate_accuracy(forecasts = methods[[m]][[h]], actuals = actuals_h)
      data.frame(horizon = h, series = rownames(a), method = m, a, row.names = NULL)
    }))
  }))
  
  # improvement over the base forecast of the same series and horizon (as before)
  accuracy_comparison$RMSE_improvement <- NA_real_
  accuracy_comparison$MAE_improvement <- NA_real_
  
  for (h in unique(accuracy_comparison$horizon)) {
    for (s in unique(accuracy_comparison$series)) {
      base_rmse <- accuracy_comparison$RMSE[accuracy_comparison$horizon == h &
                                              accuracy_comparison$series == s &
                                              accuracy_comparison$method == "Base"]
      base_mae <- accuracy_comparison$MAE[accuracy_comparison$horizon == h &
                                            accuracy_comparison$series == s &
                                            accuracy_comparison$method == "Base"]
      idx <- accuracy_comparison$horizon == h &
        accuracy_comparison$series == s &
        accuracy_comparison$method != "Base"
      accuracy_comparison$RMSE_improvement[idx] <- (base_rmse - accuracy_comparison$RMSE[idx]) / base_rmse * 100
      accuracy_comparison$MAE_improvement[idx]  <- (base_mae - accuracy_comparison$MAE[idx]) / base_mae * 100
    }
  }
  
  summary_by_series <- do.call(rbind, lapply(unique(accuracy_comparison$series), function(s) {
    x <- accuracy_comparison[accuracy_comparison$series == s & accuracy_comparison$method != "Base", ]
    data.frame(
      series = s,
      method = c("Bottom-Up", "OLS", "WLS"),
      mean_RMSE_improvement = c(mean(x$RMSE_improvement[x$method == "Bottom-Up"], na.rm = TRUE),
                                mean(x$RMSE_improvement[x$method == "OLS"], na.rm = TRUE),
                                mean(x$RMSE_improvement[x$method == "WLS"], na.rm = TRUE)),
      mean_MAE_improvement = c(mean(x$MAE_improvement[x$method == "Bottom-Up"], na.rm = TRUE),
                               mean(x$MAE_improvement[x$method == "OLS"], na.rm = TRUE),
                               mean(x$MAE_improvement[x$method == "WLS"], na.rm = TRUE))
    )
  }))
  
  list(accuracy_comparison = accuracy_comparison, summary_by_series = summary_by_series)
}


# run everything: every base model x window -------------------------------------
# results are saved after every model/window (hierarchy_results.rds), so a crash
# halfway does not cost the finished parts
hierarchy_results <- list()

hierarchy_selections <- list()

for (m in hierarchy_base_models) {
  
  # orders per series, selected once on April-July
  sel_m <- lapply(all_series, function(s) select_hierarchy_model(m, s))
  names(sel_m) <- all_series
  hierarchy_selections[[m]] <- sel_m
  res_m <- training_residuals(sel_m)
  
  for (w in c("rolling", "expanding")) {
    
    key <- paste(m, w, sep = "_")
    cat("\n================", key, "================\n")
    
    fc  <- forecast_hierarchy(m, w, sel_m)
    coh <- base_coherence(fc)
    rec <- reconcile_hierarchy(fc, res_m)
    acc <- hierarchy_accuracy(rec, fc)
    
    hierarchy_results[[key]] <- list(forecasts = fc, base_coherence = coh,
                                     reconciled = rec, accuracy = acc)
    saveRDS(hierarchy_results, "hierarchy_results.rds")
    
    #evaluation
    print(coh)                       # base forecasts: not coherent
    print(round(rec$coherence, 4))   # reconciled forecasts: deviations are gone
    print(acc$summary_by_series)     # improvement of reconciliation over the base forecasts
  }
}


# final comparison: base model x reconciliation method, per series -------------
# average RMSE/MAE over all horizons per series, method, base model and window
final_comparison <- do.call(rbind, lapply(names(hierarchy_results), function(key) {
  x <- hierarchy_results[[key]]$accuracy$accuracy_comparison
  agg <- aggregate(cbind(RMSE, MAE) ~ series + method, data = x, FUN = mean)
  agg$base_model <- sub("_(rolling|expanding)$", "", key)
  agg$window <- sub(".*_", "", key)
  agg[, c("window", "series", "base_model", "method", "RMSE", "MAE")]
}))
final_comparison <- final_comparison[order(final_comparison$window,
                                           match(final_comparison$series, hierarchy_series),
                                           final_comparison$RMSE), ]
final_comparison$RMSE <- round(final_comparison$RMSE, 1)
final_comparison$MAE <- round(final_comparison$MAE, 1)
print(final_comparison, row.names = FALSE)

# best combination (base model + method) per series and window
best_per_series <- do.call(rbind, lapply(
  split(final_comparison, list(final_comparison$series, final_comparison$window), drop = TRUE),
  function(d) d[which.min(d$RMSE), ]))
best_per_series <- best_per_series[order(best_per_series$window,
                                         match(best_per_series$series, hierarchy_series)), ]
print(best_per_series, row.names = FALSE)

# export the hierarchy results to Excel
short_key <- function(key) sub("expanding", "exp", sub("rolling", "roll", sub("sarima_", "", key)))
as.sheet <- function(M) data.frame(name = rownames(M), M, check.names = FALSE)
hierarchy_sheets <- list(final_comparison = final_comparison, best_per_series = best_per_series)
for (key in names(hierarchy_results)) {
  k <- short_key(key)
  hierarchy_sheets[[paste0(k, "_summary")]]      <- hierarchy_results[[key]]$accuracy$summary_by_series
  hierarchy_sheets[[paste0(k, "_base_coh")]]     <- hierarchy_results[[key]]$base_coherence
  hierarchy_sheets[[paste0(k, "_rec_coh")]]      <- as.sheet(round(hierarchy_results[[key]]$reconciled$coherence, 4))
}
dir.create("results", showWarnings = FALSE)
writexl::write_xlsx(hierarchy_sheets, "results/citibike_hierarchy_results.xlsx")

# conclusion forecast reconciliation --------------------------------------

#Are the base forecasts coherent? 
#No. Each series is forecast by its own model with its own orders and parameters, 
#and nothing forces the forecasts to add up; a forecast of a sum is not the sum of the forecasts. 
#The gaps are largest around the rush hours, when demand and forecast errors are largest.

#Why the gap differs per zone?
#Because the actual zone equals the sum of its actual areas, the gap equals the sum of the area forecast errors minus the zone's own forecast error. 
#It is large when the zone model and the area models make different mistakes, and it grows with the size of the series. 
#Relative to the mean, the weekly SARIMA gaps are 2% to 4% for every level. 
#The daily SARIMA's zone 1 gap is 12%, because its own zone 1 model is weak. 
#The daily SARIMA has smaller zone 2 gaps than the weekly SARIMA (89 against 135),
#yet its zone 2 forecasts are 38% less accurate (RMSE 1,163 against 717). Its area and zone models all miss the weekly pattern, so on a Saturday they all forecast a weekday rush hour; that shared mistake cancels in the gap.
#The forecasts are consistent, but consistently wrong. Coherence and accuracy are separate questions, which is why a good base model and reconciliation are both needed.

#Forecast Reconciliation
#Reconciliation makes all forecasts coherent.

#Classical methods. 
#Bottom-up forecasts only the six areas and adds them up to get the zones and NYC; 
#the upper-level base forecasts are not used. Top-down works the other way;
#-> we use bottomup as the classical method.

#Regression-based methods (Athanasopoulos et al., 2024). 
#Stack the nine base forecasts in a vector ŷ and write the hierarchy as y = S b, 
#where b holds the six areas and the summing matrix S adds them into zones and NYC.
#The weight matrix W sets how much each base forecast may be adjusted. 
# ->OLS uses W = I and treats every series as equally uncertain. 
# ->WLS uses the in-sample residual variances of each series on the diagonal, so noisy series absorb more of the adjustment. 
# ->MinT would use the full residual covariance matrix; we did not include it.
#We used csbu() and csrec() from FoReco (Girolimetto et al., 2026), with the WLS variances from each base model's own training fits.



#Does reconciliation improve accuracy? 
#Slightly, depending on the method and the level: WLS gives small but consistent gains, OLS is risky when series differ a lot in size, 
#and the gains (at most a few percent) are small compared with choosing the right base model. 
#Its main value is coherence: one consistent set of numbers for areas, zones and the city.
#Bottom-up leaves every area unchanged (0.0 by construction) and only replaces the zones and NYC by the sums of their areas: that helps zone 2,
#but slightly worsens zone 1 and NYC.
#reconciliation cannot replace a good model: the best weekly-based combination beats the best daily-based combination for every series.


# demand forecast usage ---------------------------------------------------

#The main use is overnight rebalancing: 
#a forecast made at midnight for every hour of the next day, per area, tells Citi Bike where bikes will be needed at the 8 am and 5 to 6 pm peaks. 
#Coherent forecasts let the same numbers serve planning at the area, zone and city level. 
#Rebalancing: move bikes before the morning rush to the areas where pick-ups will exceed the bikes in the docks, and plan truck routes per zone. 
#Staff and fleet scheduling: schedule rebalancing crews and mechanics by zone and hour, with most capacity before the peaks and repair work in the quiet hours between 2 and 5 am. 
#Station and dock capacity: identify areas that regularly run out of bikes at peak hours, where extra docks or bikes would add trips. . 
#Pricing and rider incentives: reward riders for returning bikes to areas forecast to run short, or adjust prices in the busiest hours.
#Fleet size and budgets: the NYC-level forecast and its trend (about +29% from April to July) feed decisions on fleet size and seasonal budgets, consistent with the area plans.


# conclusion --------------------------------------------------------------
#For Citi Bike we recommend the weekly SARIMA with a trend (drift) as the base model for every area, zone and the city, reconciled with WLS, re-estimated every night. 
#The most important lesson is that the weekly pattern drives accuracy, and that the standard automatic tools quietly fail to model it. 

#Lessons learned. 
#1. Check what the automatic search returns. auto.arima with a 168-hour season dropped the weekly part without a warning; only its order output (p, d, q without P, D, Q) and the Inf values in its trace revealed it. ets() does the same for seasons longer than 24. 
#2. Use a tool built for long seasons. The smooth package estimated the weekly SARIMA in seconds, and it lowered the forecast error by 27% to 34% compared with the models without a weekly component.
#3. Keep the trend. The unit root tests point to trend stationarity, and leaving out the trend more than doubled the forecast error (RMSE 3,116 against 1,431). Differencing to remove it (as auto.arima does without a trend) conflicts with the tests. 
#4. Test on the training sample only. Unit root tests, model selection and seasonal adjustment should not see the test month. 
#5. Choose metrics from the costs. RMSE plus overshoot counts reflect costly rush-hour misses; MAPE breaks down with zero-trip hours. 
#6. Compare many models with the MCS. Hundreds of DM tests produce chance results; the MCS gives one joint answer. 
#7. Reconcile with weights. WLS is a safe default; OLS distorts small series when series differ a lot in size. 

#Tips and tricks. 
#Compare AIC only within one package, one data set and one differencing order; use out-of-sample RMSE and MAE across classes. 
#Look at the residual ACF at lags 24 and 168: leftover spikes there show which season a model misses. 
#Fix the model orders on the training sample and re-estimate the parameters at each origin; it is fast, and it evaluates the models you actually specified.
#Save intermediate results (checkpoints) during long runs. Next steps. Add weather and holiday regressors, which would explain the outlier days. Treat hours with empty stations as censored demand. Forecast returns as well as pick-ups, and produce prediction intervals.
#Try MinT reconciliation, and evaluate over a longer test period. Test the weekly ETS out of sample (DM and MCS); it fitted well in-sample but was not part of the forecast exercise.
