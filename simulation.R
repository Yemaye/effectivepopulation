#Original comments are kept/slightly modified

################################################################################
# R-code to simulate from a Markov and Non-Markov Stochastic epidemic models
################################################################################


# This function assumes 1 initial infective and N-1 initially susceptibles
# Per-person infection rate is beta
# This function is for lockdown. Lockdown takes effect after 'lock' days. beta is reduced by 'frac' ratio
simSIR.Markov <- function(N, beta, gamma,frac,lock) {
  
  # initial number of infectives, revovered and susceptibles;
  I <- 1
  S <- N-1;
  R <-0
  # recording time (as a variable and an array);
  t <- 0;
  times <- c(t);
  
  # a vector which records the type of event (1=infection, 2=removal)
  type <- c(1);
  
  #
  infec<-c(1)
  
  while (I > 0) {
    if (t>lock){
      adjbeta<-beta*frac
    } else{
      adjbeta<-beta
    }
    # time to next event;
    t <- t + rexp(1, (adjbeta)*I*S + gamma*I);
    
    times <- append(times, ceiling(t));
    
    if (runif(1) <= adjbeta*S/(adjbeta*S + gamma)) {
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

# This function assumes no lockdown
source("C:\\Users\\madiy\\Desktop\\R_Files\\simulation.R")
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
write.table(t(infec), "infec_1000_25_cut_20.csv", row.names=FALSE,col.names = FALSE,append=F,sep=',')
write.table(t(remov), "remov_1000_25_cut_20.csv", row.names=F, col.names = F, append=F,sep=',')
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
  write.table(t(infec), "infec_1000_25_cut_20.csv", row.names=FALSE,col.names = FALSE,append=T,sep=',');
  write.table(t(remov), "remov_1000_25_cut_20.csv", row.names=F, col.names = F, append=T,sep=',');
  i<-i+1
} 
}


library("ggplot2")
#This part is for box-plot analysis for methods.
#N
simul<-read.table("sim_fitting_25_1000_lock_40.csv",header=F,sep=',')

data<-data.frame(simul[,1],simul[,2],simul[,3],simul[,4],simul[,5],simul[,6],simul[,7])
data<-data.frame(simul[,1],simul[,2],simul[,4])
names(data)[names(data) == "simul...1."] <- "Fit to I-data"
names(data)[names(data) == "simul...2."] <- "Fit to R-data"
names(data)[names(data) == "simul...3."] <- "Beta via R0"
names(data)[names(data) == "simul...4."] <- "Gamma=0.2"
names(data)[names(data) == "simul...5."] <- "Combination"
names(data)[names(data) == "simul...6."] <- "Final size"
names(data)[names(data) == "simul...7."] <- "Imax"

par(cex.axis=1)
boxplot(data)
abline(h=1000,col="red")

#if you want to look separately
boxplot(simul[,1],main="Fit to I-data")
boxplot(simul[,2],main="Fit to R-data")
boxplot(simul[,3],main="Beta via R0")
boxplot(simul[,4], main="Gamma=0.2")
boxplot(simul[,5],main="Beta via R0, gamma=0.2")
boxplot(simul[,6],main="Final size")
boxplot(simul[,7],main="Imax")
summary(simul[,2])

#sse
sse<-read.table("sse_fitting_25_1000.csv",header=F,sep=',')
data_sse<-data.frame(sse[,1],sse[,2],sse[,3],sse[,4],sse[,5],sse[,6],sse[,7])
names(data_sse)[names(data_sse) == "sse...1."] <- "Fit to I-data"
names(data_sse)[names(data_sse) == "sse...2."] <- "Fit to R-data"
names(data_sse)[names(data_sse) == "sse...3."] <- "Beta via R0"
names(data_sse)[names(data_sse) == "sse...4."] <- "Gamma=0.2"
names(data_sse)[names(data_sse) == "sse...5."] <- "Combination"
names(data_sse)[names(data_sse) == "sse...6."] <- "Final size"
names(data_sse)[names(data_sse) == "sse...7."] <- "Imax"
boxplot(data_sse,log="y")


#seprately
boxplot(sse[,1],main="Fit to I-data",log="y")
boxplot(sse[,2],main="Fit to R-data",log="y")
boxplot(sse[,3],main="Beta via R0",log="y")
boxplot(sse[,4], main="Gamma=0.2",log="y")
boxplot(sse[,5],main="combination",log="y")
boxplot(sse[,6],main="Final size",log="y")
boxplot(sse[,7],main="Imax",log="y")

#beta or gamma
betas<-read.table("beta_fitting_25_1000.csv",header=F,sep=',')
data_beta<-data.frame(betas[,1],betas[,2],betas[,3],betas[,4],betas[,5])
names(data_beta)[names(data_beta) == "betas...1."] <- "Fit to I-data"
names(data_beta)[names(data_beta) == "betas...2."] <- "Fit to R-data"
names(data_beta)[names(data_beta) == "betas...3."] <- "Beta via R0"
names(data_beta)[names(data_beta) == "betas...4."] <- "Gamma=0.2"
names(data_beta)[names(data_beta) == "betas...5."] <- "Combination"
boxplot(data_beta)
abline(h=0.0005,col="red")
