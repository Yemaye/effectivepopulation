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


