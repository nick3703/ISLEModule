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


days_df <- data.frame(days=daysfromstart_mod)

binned_vals <- days_df %>% mutate(new_bin = cut(daysfromstart_mod,seq(0,100,5))) %>% 
  group_by(new_bin)%>%summarize(count=n())

binned_vals <- c(0,0,0,binned_vals$count)
binned_vals[4] <- 8
binned_vals <- binned_vals[-length(binned_vals)]
mean(binned_vals)
var(binned_vals)

act<-var(binned_vals)/mean(binned_vals)

M <- 5000
resp<-c()

for(m in 1:M){

  simRecs <- hawkes(100, fun = optMarathon$par[1], repr = optMarathon$par[2], 
                    family = "exp", rate = optMarathon$par[3]) 
  
  
  simDat <- data.frame(simDays = simRecs$p)
  sim_binned_vals <- simDat %>% mutate(new_bin = cut(simDays,seq(0,100,5))) %>% 
    group_by(new_bin)%>%summarize(count=n())
  
  missing_bins <- 20-nrow(sim_binned_vals)
  
  sim_bin_vals <- c(rep(0,missing_bins),sim_binned_vals$count)
  
  resp[m]<-var(sim_bin_vals)/mean(sim_bin_vals)
  
}


simulation_df <- data.frame(Var_Mean = resp)

simulation_df %>% ggplot(aes(x=Var_Mean))+geom_histogram()+
  theme_bw()+geom_vline(aes(xintercept = act,color="Variance_Mean_ratio"),size=1.5)+
  ggtitle("Simulated Variance to Mean Ratio using Hawkes Process")+
  xlab("Empirical Distribution of Variance to Mean Ratio")+
  scale_color_manual(name = "Value from Data", values = c(Variance_Mean_ratio = "red"))





lam<-mean(binned_vals)

M <- 5000
resp_pois<-c()

for(m in 1:M){
  
  simPois<-rpois(length(binned_vals),lam)
  
  resp_pois[m]<-var(simPois)/mean(simPois)
  
}



simulation_df_pois <- data.frame(Var_Mean = resp_pois)

simulation_df_pois %>% ggplot(aes(x=Var_Mean))+geom_histogram()+
  theme_bw()+geom_vline(aes(xintercept = act,color="Variance_Mean_ratio"),size=1.5)+
  ggtitle("Simulated Variance to Mean Ratio using Poisson Process")+
  xlab("Empirical Distribution of Variance to Mean Ratio")+
  scale_color_manual(name = "Value from Data", values = c(Variance_Mean_ratio = "red"))




1.26*4.2*exp(-4.2*1)


base_x<-seq(0,100,.01)
base_y<-rep(.43,length(base_x))


base_y_exc <- .126*4.2*exp(-4.2*base_x)+base_y


excite <- base_y
for(j in 1:length(daysfromstart_mod)){
  for(i in 1:length(base_x)){
    if(daysfromstart_mod[j]>base_x[i]){
      excite[i] <- excite[i]+0
    }else{
      excite[i] <- excite[i]+.126*4.2*exp(-4.2*(base_x[i]-daysfromstart_mod[j]))
    }
  }
}

plot(base_x,excite,type="l",ylim=c(0,2),col="red")
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
