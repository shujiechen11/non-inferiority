rm(list = ls())
##this r script is used to store the function used in the simulation
# delta compute function to calculate delta based on the historical information
deltacal <- function(lambdap,lambdaac,num_event,signal)
{
#lamobdap:incidence rate of placebo in the historical trial
#lambdaac:incidence rate of control treatment in the historical trial
#num_event: number of events in the historical trial
#signal: to distinct the deltas 0:delta_0,1:delta_1, otherwise delta
  hazard_ratio = lambdaac/lambdap
  na = num_event/(1+hazard_ratio)*hazard_ratio
  np = num_event/(1+hazard_ratio)
  #p = rpois(1,100*lambdap)
  #c = rpois(1,100*lambdaac)
  log_hr_ul0.95_fromal = log(hazard_ratio)+1.96*sqrt(1/na+1/np)
  #log_hr_ul0.95_fromal = log(c/p)+1.96*sqrt(4/num_event)
  #quantile way
  #c_p<-c/p
  #c_p[is.infinite(c_p)]=100#outliers handling
  #c_p[is.na(c_p)]=1#outliers handling
  #cic_p<-quantile(c_p,0.975)
  #hr_ul0.95_fromcal <- quantile(c_p,0.975)
  hr_ul0.95_fromcal = exp(log_hr_ul0.95_fromal)
  if(signal == 0)
  {
    return(hr_ul0.95_fromcal^{-1}*0.7^{1/2})
  }
  if(signal == 1)
  {
    return(hr_ul0.95_fromcal^{-1}*0.7)
  }
  return(hr_ul0.95_fromcal^{-1/2})
}
################################################################################################3
#################################total person years calculation function under the null hypothesis
personyears <- function(lambdac,lambdat,alpha,beta,delta,delta1)
{
  #lambdac: incidence rate of control treatment in the non-inferiority trial
  #lambdat: incidence rate of experimental treatment in the non-inferiority trial
  #alpha: type one error
  #beta: power
  #delta: non-inferiority margin in the h0 hypothesis
  #delta1: non-inferiority margin in the alternative hypothesis
  z_alpha = qnorm(1-alpha/2)
  z_beta = qnorm(beta)
  N = 2*(z_beta*sqrt(delta*(delta1+1)/(delta1*(delta+1)))+z_alpha)^2/((log(delta)-log(delta1))^2)*(1/lambdac+1/(lambdac*delta))
  return(ceiling(N))
}
#function to compute the type one error or power in the non-inferiority trial 
typeone_power_test <- function(lambdac, lambdat,personyears,delta,alpha)
{ 
  N = personyears/2
  z_alpha = qnorm(1-alpha/2)
  t = rpois(10000, N*lambdat)
  c = rpois(10000, N*lambdac)
  tupper = log(t/N)-log(c/N)-log(delta)
  tupper[is.na(tupper)]=-log(delta)#deal with the situlation log0-log0 which gives the na in tupper
  #### the way to deal with na makes the power unrealistic
  tlower=sqrt(1/(c*delta)+1/(c))
  tn = tupper/tlower
  #deal with the situation that tupper is Inf and tlower is Inf
  tn[is.na(tn) & tupper>0] = 1
  #deal with the situation that tupper is -Inf and tlower is Inf
  tn[is.na(tn) & tupper<0] = -1
  outcome = tn< -z_alpha
  result = sum(outcome)/10000
  return(result)
}
simulation <- function(lambdat,lambdac,alpha,beta,lambdap,lambdaac,num_event,signal,delta1)
{ 
  #calculate the ni margin
  delta = deltacal(lambdap,lambdaac,num_event,signal)
  # firstly, calculate the person years to generate poisson distribution
  personyears = personyears(lambdac,lambdat,alpha,beta,delta,delta1)
  # secondly, test the margin
  result = typeone_power_test(lambdac, lambdat,personyears,delta,alpha)
 # print(paste("type one error test or power of tn is", result,
  #            ' while the ni margin is',delta, ' and the incidence rate of experimental treatment is', lambdat, 
    #          ' and the incidence rate of control treatment is',lambdac))
  return(result)
}

#############################################################################################################################3
#######################################total person years calculation function without the null hypothesis
personyears_2 <- function(lambdac,lambdat,alpha,beta,delta,delta1)
{
  #lambdac: incidence rate of control treatment in the non-inferiority trial
  #lambdat: incidence rate of experimental treatment in the non-inferiority trial
  #alpha: type one error
  #beta: power
  #delta: non-inferiority margin in the h0 hypothesis
  #delta1: non-inferiority margin in the alternative hypothesis
  z_alpha = qnorm(1-alpha/2)
  z_beta = qnorm(beta)
  N = 2*(z_beta+z_alpha)^2/((log(delta)-log(delta1))^2)*(1/lambdac+1/lambdat)
  return(ceiling(N))
}
typeone_power_test_2 <- function(lambdac, lambdat,personyears,delta,alpha)
{ 
  N = personyears/2
  z_alpha = qnorm(1-alpha/2)
  t = rpois(10000, N*lambdat)
  c = rpois(10000, N*lambdac)
  tupper = log(t/N)-log(c/N)-log(delta)
  tupper[is.na(tupper)]=-log(delta)#deal with the situlation log0-log0 which gives the na in tupper
  tlower=sqrt(1/(t)+1/(c))
  tn = tupper/tlower
  #deal with the situation that tupper is Inf and tlower is Inf
  tn[is.na(tn) & tupper>0] = 1
  #deal with the situation that tupper is -Inf and tlower is Inf
  tn[is.na(tn) & tupper<0] = -1
  outcome = tn< -z_alpha
  result = sum(outcome)/10000
  return(result)
}
simulation_2 <- function(lambdat,lambdac,alpha,beta,lambdap,lambdaac,num_event,signal,delta1)
{ 
  #calculate the ni margin
  delta = deltacal(lambdap,lambdaac,num_event,signal)
  # firstly, calculate the person years to generate poisson distribution
  personyears = personyears_2(lambdac,lambdat,alpha,beta,delta,delta1)
  # secondly, test the margin
  result = typeone_power_test_2(lambdac, lambdat,personyears,delta,alpha)
  # print(paste("type one error test or power of tn is", result,
  #            ' while the ni margin is',delta, ' and the incidence rate of experimental treatment is', lambdat, 
  #          ' and the incidence rate of control treatment is',lambdac))
  return(result)
}



###########################################way one to calculate the total person years#########################################################
signal = -1
signal_0 = 0
signal_1 = 1
delta1 = 1
alpha = 0.05
beta = 0.9
lambdaacs <- c(0.01,0.01,0.02,0.02,0.04,0.04,0.06,0.06,0.08,0.08)
lambdap <- 0.2
num_events <- c(175,350,175,350,175,350,175,350,175,350)
lambdacs <- c(0.01,0.01,0.02,0.02,0.04,0.04,0.06,0.06,0.08,0.08)
ratio_set <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
y <- c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
#par(mfcol=c(1,1))
for (i in 1:10)
{
  lambdac <- lambdacs[i]
  lambdaac <- lambdaacs[i]
  num_event <- num_events[i]
  simulations <- c()
  simulations_0 <- c()
  simulations_1 <- c()
  for (j in ratio_set)
  {
    simulations <- append(simulations,simulation(lambdac*j,lambdac,alpha,beta,lambdap,lambdaac,num_event,signal,delta1))
    simulations_0 <- append(simulations_0,simulation(lambdac*j,lambdac,alpha,beta,lambdap,lambdaac,num_event,signal_0,delta1))
    simulations_1 <- append(simulations_1,simulation(lambdac*j,lambdac,alpha,beta,lambdap,lambdaac,num_event,signal_1,delta1))
  }
  plot(simulations ,type = "o",col = "red", xlab = "EXP/AC", ylab = "Possibility to reject the null hypothesis", xaxt = 'n',
       ylim = c(0,0.95),xlim = c(1,20) ,main = "NI TRIAL (first way)")
  lines(simulations_0, type = "o", col = "blue")
  lines(simulations_1, type = "o", col = "purple")
  abline(h=0.05,col='black')
  abline(h=0.9,col='black')
  abline(v=deltacal(lambdap,lambdaac,num_event,signal),col='red')
  abline(v=deltacal(lambdap,lambdaac,num_event,signal_0),col='blue')
  abline(v=deltacal(lambdap,lambdaac,num_event,signal_1),col='purple')
  abline(v=(lambdaac/lambdap)^{-1/2},col='red',lty=3)
  abline(v=(lambdaac/lambdap)^{-1}*(0.7)^{1/2},col='blue',lty=3)
  abline(v=(lambdaac/lambdap)^{-1}*(0.7),col='purple',lty=3)
  axis(1,ratio_set) #x scale
  axis(2,y)#y scale
  legend("topright",legend=c("delta","delta0","delta1"),col=c("red","blue","purple"),lty=1,lwd=2) 
}



###########################################way two to calculate the total person years#########################################################
signal = -1
signal_0 = 0
signal_1 = 1
delta1 = 1
alpha = 0.05
beta = 0.9
lambdaacs <- c(0.01,0.01,0.02,0.02,0.04,0.04,0.06,0.06,0.08,0.08)
lambdap <- 0.2
num_events <- c(175,350,175,350,175,350,175,350,175,350)
lambdacs <- c(0.01,0.01,0.02,0.02,0.04,0.04,0.06,0.06,0.08,0.08)
ratio_set <- c(1,2,3,4,5,6,7,8,9,10,11,12)
y <- c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
#par(mfcol=c(1,1))
for (i in 1:10)
{
  lambdac <- lambdacs[i]
  lambdaac <- lambdaacs[i]
  num_event <- num_events[i]
  simulations <- c()
  simulations_0 <- c()
  simulations_1 <- c()
  for (j in ratio_set)
  {
    simulations <- append(simulations,simulation_2(lambdac*j,lambdac,alpha,beta,lambdap,lambdaac,num_event,signal,delta1))
    simulations_0 <- append(simulations_0,simulation_2(lambdac*j,lambdac,alpha,beta,lambdap,lambdaac,num_event,signal_0,delta1))
    simulations_1 <- append(simulations_1,simulation_2(lambdac*j,lambdac,alpha,beta,lambdap,lambdaac,num_event,signal_1,delta1))
  }
  plot(simulations ,type = "o",col = "red", xlab = "EXP/AC", ylab = "Possibility to reject the null hypothesis", xaxt = 'n',
       ylim = c(0,0.95),xlim = c(1,20) ,main = "NI TRIAL (second way)")
  lines(simulations_0, type = "o", col = "blue")
  lines(simulations_1, type = "o", col = "purple")
  abline(h=0.05,col='black')
  abline(h=0.9,col='black')
  abline(v=deltacal(lambdap,lambdaac,num_event,signal),col='red')
  abline(v=deltacal(lambdap,lambdaac,num_event,signal_0),col='blue')
  abline(v=deltacal(lambdap,lambdaac,num_event,signal_1),col='purple')
  abline(v=(lambdaac/lambdap)^{-1/2},col='red',lty=3)
  abline(v=(lambdaac/lambdap)^{-1}*(0.7)^{1/2},col='blue',lty=3)
  abline(v=(lambdaac/lambdap)^{-1}*(0.7),col='purple',lty=3)
  axis(1,ratio_set) #x scale
  axis(2,y)#y scale
  legend("topright",legend=c("delta","delta0","delta1"),col=c("red","blue","purple"),lty=1,lwd=2) 
}

#############################plot to summarize the power (first way)###########################################
simulations <- c()
simulations_0 <- c()
simulations_1 <- c()
lambdaacs <- c(0.01,0.02,0.04,0.06,0.08)
num_events1 <- c(175,175,175,175,175)
for (i  in 1:5)
{
  lambdac <- lambdacs[i]
  lambdaac <- lambdaacs[i]
  num_event <- num_events1[i]
  simulations <- append(simulations,simulation(lambdac,lambdac,alpha,beta,lambdap,lambdaac,num_event,signal,delta1))
  simulations_0 <- append(simulations_0,simulation(lambdac,lambdac,alpha,beta,lambdap,lambdaac,num_event,signal_0,delta1))
  simulations_1 <- append(simulations_1,simulation(lambdac,lambdac,alpha,beta,lambdap,lambdaac,num_event,signal_1,delta1))
}
plot(simulations ,type = "o",col = "red", xlab = "AC", ylab = "Possibility to reject the null hypothesis", xaxt = 'n',
     ylim = c(0.75,0.95), main = "Power(first way)")
lines(simulations_0, type = "o", col = "blue")
lines(simulations_1, type = "o", col = "purple")
axis(1,lambdacs) #x scale
legend("topright",legend=c("delta","delta0","delta1"),col=c("red","blue","purple"),lty=1,lwd=2) 

#############################plot to summarize the power (second way)###########################################
simulations <- c()
simulations_0 <- c()
simulations_1 <- c()
lambdaacs <- c(0.01,0.02,0.04,0.06,0.08)
num_events1 <- c(175,175,175,175,175)
for (i  in 1:5)
{
  lambdac <- lambdacs[i]
  lambdaac <- lambdaacs[i]
  num_event <- num_events1[i]
  simulations <- append(simulations,simulation_2(lambdac,lambdac,alpha,beta,lambdap,lambdaac,num_event,signal,delta1))
  simulations_0 <- append(simulations_0,simulation_2(lambdac,lambdac,alpha,beta,lambdap,lambdaac,num_event,signal_0,delta1))
  simulations_1 <- append(simulations_1,simulation_2(lambdac,lambdac,alpha,beta,lambdap,lambdaac,num_event,signal_1,delta1))
}
plot(simulations ,type = "o",col = "red", xlab = "AC", ylab = "Possibility to reject the null hypothesis", xaxt = 'n',
     ylim = c(0.75,0.95), main = "Power(second way)")
lines(simulations_0, type = "o", col = "blue")
lines(simulations_1, type = "o", col = "purple")
axis(1,lambdacs) #x scale
legend("topright",legend=c("delta","delta0","delta1"),col=c("red","blue","purple"),lty=1,lwd=2) 

