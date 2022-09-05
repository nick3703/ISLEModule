library(readr)
library(hawkesbow)
record_table_mod<-read_rds("record_table_mod.rds")

days_between = as.numeric(diff(record_table_mod$Date_ymd))
daysfromstart <- cumsum(days_between)
daysfromstart <- c(0,daysfromstart)  ### get data in terms of days from first record (time 0)

daysfromstart_mod <- daysfromstart/max(daysfromstart)*100

set.seed(21)
optMarathon<-mle(daysfromstart_mod,"Exponential",100)  # end date picked number greater than longest times
optMarathon$par


compensate<-compensator(daysfromstart_mod, t=seq(0,100), fun = optMarathon$par[1], repr = optMarathon$par[2], 
                        family = "exp", rate = optMarathon$par[3])
plot(c(0:100),compensate-optMarathon$par[1]*seq(0,100),main="Full time period", type = "l")
points(daysfromstart_mod,rep(0,length(daysfromstart_mod)), pch=3)

optMarathon$model$loglik(daysfromstart_mod, 100)  #  log-likelihood

summary(optMarathon)
optMarathon$events
optMarathon$end

optMarathon$model$param
optMarathon$model$mean()  # expected value
optMarathon$model$dmean()  # Jacobian matrix of expected value
optMarathon$model$ddmean()  # Hessian matrix of expected value
optMarathon$model$loglik(daysfromstart_mod, optMarathon$end)  #  log-likelihood
optMarathon$model$dloglik(daysfromstart_mod, optMarathon$end)  # Jacobian matrix of log-lik
optMarathon$model$ddloglik(daysfromstart_mod, optMarathon$end) # hessian matrix of log-lik

resid<-residuals(daysfromstart_mod, fun = optMarathon$par[1], repr = optMarathon$par[2], 
                 family = "exp", rate = optMarathon$par[3])
resid
plot(resid)
abline(0, 1, col="red", lty="dashed")



compensate<-compensator(daysfromstart_mod, t=seq(0,100,by=10), fun = optMarathon$par[1], repr = optMarathon$par[2], 
                        family = "exp", rate = optMarathon$par[3])
plot(c(0:100),compensate,main="Full time period", type = "l")


points(daysfromstart_mod,rep(0,length(daysfromstart_mod)), pch=3)



simRecs <- hawkes(100, fun = optMarathon$par[1], repr = optMarathon$par[2], 
                  family = "exp", rate = optMarathon$par[3]) 

plot(simRecs, intensity = FALSE)
plot(simRecs, intensity = TRUE)
