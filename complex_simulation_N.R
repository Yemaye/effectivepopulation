#Multi-patch simulation, extension of simple one, it follows the original lines with some modifications
simSIRComplex3.Markov <- function(N, betamatrix, gammas) {
  n<-length(N)
  # initial number of infectives and susceptibles;
  S=N
  S[1]=N[1]-1
  I=rep(0,n)
  I[1]=1
  R=rep(0,n)
  # recording time;
  periods<-rep(0,n)
  t<-0;
  times<-c(t);
  # a vector which records the type of event (1=infection, 2=removal)
  type <- c(1);
  
  #
  
  while (any(I>0)) {
    # time to next event;
    for(i in 1:n){
      if(S[i]*sum(betamatrix[1:n,i]*I)+gammas[i]*I[i]>0){
        periods[i] <-rexp(1,S[i]*sum(betamatrix[1:n,i]*I)+gammas[i]*I[i])
      } else {
        periods[i]<-99999
      }
    }
    k=which.min(periods);
    if (runif(1) <= (S[k]*sum(betamatrix[1:n,k]*I)/(S[k]*sum(betamatrix[1:n,k]*I)+gammas[k]*I[k])) && S[k]>0) {
      # infection
      I[k] <- I[k]+1;
      S[k] <- S[k]-1;
      type <- append(type, k);
      t<-t+min(periods);
      times<-append(times,ceiling(t));
    }
    else if (I[k]>0){
      #removal
      I[k] <- I[k]-1
      R[k]<-R[k]+1
      type <- append(type, k+n);
      t<-t+min(periods);
      times<-append(times,ceiling(t));
    } 
  }
  
  # record the final size , i.e. the number of initially susceptlbles who contracted the disease sometime during the epidemic.
  
  # record the times of events (infections/removals) as well as the type
  # add the following lines in each of the simulation functions
  
  duration <- max(times)
  
  # modify the existing `res` objetct to store the final size and the duration
  
  res <- list("t"=times, "duration" = duration,"type"=type)
  res
}

#examples of parameters, use inputs appropriately
betas=matrix(c(0.0005,0.00001,0.0000001,
               0.00002,0.0005,0.00001,
               0.0000001,0.00001,0.0005),ncol=3)
betas=matrix(c(0.0006,0.000067,0.00005,0.00004,0.000033,
               0.0001,0.0004,0.00005,0.00004,0.000033,
               0.0001,0.000067,0.0003,0.00004,0.000033,
               0.0001,0.000067,0.00005,0.00024,0.000033,
               0.0001,0.000067,0.00005,0.00004,0.0002),ncol=5)
betas=matrix(runif(25,0.00005,0.0001),nrow=5)
Ns=c(1000,1000,1000)
Ns=c(500,750,1000,1250,1500)
gammas=runif(3,0.2,0.2)
gammas=runif(5,0.2,0.2)
#example of simulation
trial<-simSIRComplex3.Markov(Ns,betas,gammas)
trial
#looking at the dynamics
infec<-matrix(rep(0,length(Ns)*(trial$duration+1)),nrow=length(Ns))
infec[1,1]=1
remov<-matrix(rep(0,length(Ns)*(trial$duration+1)),nrow=length(Ns))
for(i in 1:trial$duration){
  if (i %in% trial$t){
    for(j in 1:length(Ns)){
      remov[j,i+1]<-remov[j,i]+sum(trial$type[trial$t==i]==(j+length(Ns)));
      infec[j,i+1]<-infec[j,i]+sum(trial$type[trial$t==i]==(j))-sum(trial$type[trial$t==i]==(j+length(Ns)));
    }
  } 
  else{
    for(j in 1:length(Ns)){
    remov[j,i+1]<-remov[j,i];
    infec[j,i+1]<-infec[j,i];
    }
  }  
}

#plotting
plot(infec[1,1:trial$duration+1])
points(infec[2,1:trial$duration+1])
points(infec[3,1:trial$duration+1])
points(infec[4,1:trial$duration+1])
points(infec[5,1:trial$duration+1])
inf=colSums(infec)
plot(inf)
rem=colSums(remov)
plot(rem)

#recording
write.table(t(infec[1,]), "complex/com_infec_3_2_1.csv", row.names=FALSE,col.names = FALSE,append=F,sep=',')
write.table(t(remov[1,]), "complex/com_remov_3_2_1.csv", row.names=F, col.names = F, append=F,sep=',')
write.table(t(infec[2,]), "complex/com_infec_3_2_2.csv", row.names=FALSE,col.names = FALSE,append=F,sep=',')
write.table(t(remov[2,]), "complex/com_remov_3_2_2.csv", row.names=F, col.names = F, append=F,sep=',')
write.table(t(infec[3,]), "complex/com_infec_3_2_3.csv", row.names=FALSE,col.names = FALSE,append=F,sep=',')
write.table(t(remov[3,]), "complex/com_remov_3_2_3.csv", row.names=F, col.names = F, append=F,sep=',')
write.table(t(inf), "complex/com_infec_3_2.csv", row.names=F, col.names = F, append=F,sep=',')
write.table(t(rem), "complex/com_remov_3_2.csv", row.names=F, col.names = F, append=F,sep=',')

#simulating 999 additionally
k<-0
while (k<999){
  trial=simSIRComplex3.Markov(Ns,betas,gammas)
  infec<-matrix(rep(0,length(Ns)*(trial$duration+1)),nrow=length(Ns))
  infec[1,1]=1
  remov<-matrix(rep(0,length(Ns)*(trial$duration+1)),nrow=length(Ns))
  for(i in 1:trial$duration){
    if (i %in% trial$t){
      for(j in 1:length(Ns)){
        remov[j,i+1]<-remov[j,i]+sum(trial$type[trial$t==i]==(j+length(Ns)));
        infec[j,i+1]<-infec[j,i]+sum(trial$type[trial$t==i]==(j))-sum(trial$type[trial$t==i]==(j+length(Ns)));
      }
    } 
    else{
      for(j in 1:length(Ns)){
        remov[j,i+1]<-remov[j,i];
        infec[j,i+1]<-infec[j,i];
      }
    }  
  }
  if (max(infec)>50 && (sum(infec[1:length(Ns),1:10])+sum(remov[1:length(Ns),1:10]))>10){
    write.table(t(infec[1,]), "complex/com_infec_3_2_1.csv", row.names=FALSE,col.names = FALSE,append=T,sep=',')
    write.table(t(remov[1,]), "complex/com_remov_3_2_1.csv", row.names=F, col.names = F, append=T,sep=',')
    write.table(t(infec[2,]), "complex/com_infec_3_2_2.csv", row.names=FALSE,col.names = FALSE,append=T,sep=',')
    write.table(t(remov[2,]), "complex/com_remov_3_2_2.csv", row.names=F, col.names = F, append=T,sep=',')
    write.table(t(infec[3,]), "complex/com_infec_3_2_3.csv", row.names=FALSE,col.names = FALSE,append=T,sep=',')
    write.table(t(remov[3,]), "complex/com_remov_3_2_3.csv", row.names=F, col.names = F, append=T,sep=',')
    inf=colSums(infec)
    rem=colSums(remov)
    write.table(t(inf), "complex/com_infec_3_2.csv", row.names=F, col.names = F, append=T,sep=',')
    write.table(t(rem), "complex/com_remov_3_2.csv", row.names=F, col.names = F, append=T,sep=',')
    k<-k+1
  } 
}



