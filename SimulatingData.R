#Simple MCMC simulation of SIR dinyamics
simSIR.Markov <- function(N, beta, gamma) {
  
  # initial number of infectives and susceptibles;
  I <- 1
  S <- N-1;
  R<-0
  # recording time;
  t <- 0;
  times <- c(t);
  
  # a vector which records the type of event (1=infection, 2=removal)
  type <- c(1);
  
  #
  infec<-c(1)
  
  while (I > 0) {
    
    # time to next event;
    t <- t + rexp(1, (beta)*I*S + gamma*I);
    
    times <- append(times, ceiling(t));
    
    
    if (runif(1) <= beta*S/(beta*S + gamma)) {
      # infection
      I <- I+1;
      S <- S-1;
      type <- append(type, 1);
    }
    else {
      #removal
      I <- I-1
      R<-R+1
      type <- append(type, 2);
    }
  }
  
  # record the final size , i.e. the number of initially susceptlbles who contracted the disease sometime during the epidemic.
  
  # record the times of events (infections/removals) as well as the type
  # add the following lines in each of the simulation functions
  
  final.size <- sum(type==1) - 1
  duration <- max(times)
  
  # modify the existing `res` object to store the final size and the duration
  
  res <- list("t"=times, "final.size"=final.size, "duration" = duration,"type"=type)
  res
  
}

# Example
trial=simSIR.Markov(1000,0.0005,0.2,0.5,20)
trial
#Looking at the daily infected and removed
infec<-c(1)
remov<-c(0)
for(i in 1:trial$duration){
  if (i %in% trial$t){
    remov<-append(remov,remov[length(remov)]+sum(trial$type[trial$t==i]==2));
    infec<-append(infec,infec[length(infec)]+sum(trial$type[trial$t==i]==1)-sum(trial$type[trial$t==i]==2))
  } 
  else{
    remov<-append(remov,remov[length(remov)]);
    infec<-append(infec,infec[length(infec)])
  }  
}
plot(c(0:trial$duration),infec)

#Write the table, simulated additionally 999 times
write.table(t(infec), "infec_1000_25.csv", row.names=FALSE,col.names = FALSE,append=F,sep=',')
write.table(t(remov), "remov_1000_25.csv", row.names=F, col.names = F, append=F,sep=',')
i<-0
while (i<999){
  trial=simSIR.Markov(1000,0.0005,0.2,0.5,20)
  infec<-c(1)
  remov<-c(0)
  for(j in 1:trial$duration){
    if (j %in% trial$t){
      remov<-append(remov,remov[length(remov)]+sum(trial$type[trial$t==j]==2));
      infec<-append(infec,infec[length(infec)]+sum(trial$type[trial$t==j]==1)-sum(trial$type[trial$t==j]==2))
    } 
    else{
      remov<-append(remov,remov[length(remov)]);
      infec<-append(infec,infec[length(infec)])
    }  
  }
  if (max(infec)>50 && (sum(infec[1:10])+sum(remov[1:10]))>10){
    write.table(t(infec), "infec_1000_25.csv", row.names=FALSE,col.names = FALSE,append=T,sep=',');
    write.table(t(remov), "remov_1000_25.csv", row.names=F, col.names = F, append=T,sep=',');
    i<-i+1
  } 
}

library(plyr)
library(MASS)
library(smfsb)
library(ramcmc) 
library(SDSMCMC)
library(tictoc)

#Directly copied example 3 from https://github.com/cbskust/SDS.Epidemic and slightly modfied
# Third example data: simulation data. The data is simulated by Sellke construction.
# This example has two MCMC simulations using synthetic data by Sellke construction. 
# Running the other methods on the paper except SDS approach

burn = 1000; #  burning period
thin =1; # tinning period 
nrepeat = 2000; # number of posterior sample;  
sim.num= burn + nrepeat * thin; # total number of simulation 

#initial parameter setting
k1 = 0.5; k2 = 0.2 ; k3 = 0.001; n=1000; T.max = 200; 
beta=k1; gamma=k2; rho=k3;
pop.data = Sellke(n=n, rho=k3, beta=k1, gamma=k2, Tmax = T.max)

#converting Sellke epidemic data to SIR trajectory
emp.sir <- Sellke.to.trajectory(pop.data, Tmax = T.max)

#saving IR-data
discrI<-c(1)
for (i in 1:length(emp.sir$I)){
  discrI<-append(discrI,tail(emp.sir$I[ceiling(emp.sir$time)==i],n=1))
} 
write.table(t(discrI), "infec_Sellke_1000_25.csv", row.names=FALSE,col.names = FALSE,append=T,sep=',')


discrR<-c(0)
for (i in 1:length(emp.sir$R)){
  discrR<-append(discrR,tail(emp.sir$R[ceiling(emp.sir$time)==i],n=1))
} 
write.table(t(discrR), "remov_Sellke_1000_25.csv", row.names=FALSE,col.names = FALSE,append=T,sep=',')

#Saving SDS data
sds  <- SDS.MCMC(data = pop.data, Tmax=T.max, fitn = T, nrepeat = sim.num, 
                 prior.a=c(0.001,0.001,0.001), prior.b=c(0.001,0.001,0.001), ic = c(k1, k2, k3, n))
write.table(t(c(mean(sds[,1]),mean(sds[,2]),mean(sds[,4]))),"sds_Sellke_25_1000.csv",row.names = F,
            col.names = F,append = F,sep=",")
#Combine the third column from each setting into one, named "Sellke_N.csv" to compare.

#If more simulations are needed
i<-0
while (i<900){
  pop.data = Sellke(n=n, rho=k3, beta=k1, gamma=k2, Tmax = T.max)
  emp.sir <- Sellke.to.trajectory(pop.data, Tmax = T.max)
  infec<-c(1)
  remov<-c(0)
  for (j in 1:length(emp.sir$I)){
    infec<-append(infec,tail(emp.sir$I[ceiling(emp.sir$time)==j],n=1))
  } 
  for (j in 1:length(emp.sir$R)){
    remov<-append(remov,tail(emp.sir$R[ceiling(emp.sir$time)==j],n=1))
  } 
  if (max(infec)>50 && (sum(infec[1:10])+sum(remov[1:10]))>10&&length(infec)<=150){
    write.table(t(infec), "infec_Sellke_25_1000.csv", row.names=FALSE,col.names = FALSE,append=T,sep=',');
    write.table(t(remov), "remov_Sellke_25_1000.csv", row.names=F, col.names = F, append=T,sep=',');
    #sds  <- SDS.MCMC(data = pop.data, Tmax=T.max, fitn = T, nrepeat = sim.num, 
    #                 prior.a=c(0.001,0.001,0.001), prior.b=c(0.001,0.001,0.001), ic = c(k1, k2, k3, n));
    #write.table(t(c(mean(sds[,1]),mean(sds[,2]),mean(sds[,4]))),"sds_Sellke_25_1000.csv",row.names = F,
    #            col.names = F,append = T,sep=",");
    i<-i+1
  } 
}



