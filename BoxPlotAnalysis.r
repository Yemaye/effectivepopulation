library("ggplot2")
#This part is for box-plot analysis for methods.
#N
simul<-read.table("complex_results\\complex_sim_n_3_1.csv",header=F,sep=',')

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
abline(h=3000,col="red")

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

sse<-read.table("complex_results\\complex_sim_sse_3_1.csv",header=F,sep=',')

data_sse<-data.frame(sse[,1],sse[,2],sse[,3],sse[,4],sse[,5],sse[,6],sse[,7])
names(data_sse)[names(data_sse) == "sse...1."] <- "Fit to I-data"
names(data_sse)[names(data_sse) == "sse...2."] <- "Fit to R-data"
names(data_sse)[names(data_sse) == "sse...3."] <- "Beta via R0"
names(data_sse)[names(data_sse) == "sse...4."] <- "Gamma=0.2"
names(data_sse)[names(data_sse) == "sse...5."] <- "Combination"
names(data_sse)[names(data_sse) == "sse...6."] <- "Final size"
names(data_sse)[names(data_sse) == "sse...7."] <- "Imax"
boxplot(data_sse,log="y")


#separately
boxplot(sse[,1],main="Fit to I-data",log="y")
boxplot(sse[,2],main="Fit to R-data",log="y")
boxplot(sse[,3],main="Beta via R0",log="y")
boxplot(sse[,4], main="Gamma=0.2",log="y")
boxplot(sse[,5],main="combination",log="y")
boxplot(sse[,6],main="Final size",log="y")
boxplot(sse[,7],main="Imax",log="y")

#beta or gamma
betas<-read.table("complex_results\\complex_sim_gamma_3_1.csv",header=F,sep=',')
data_beta<-data.frame(betas[,2],betas[,5],betas[,8])
names(data_beta)[names(data_beta) == "betas...2."] <- "Fit to I-data"
names(data_beta)[names(data_beta) == "betas...5."] <- "Fit to R-data"

names(data_beta)[names(data_beta) == "betas...8."] <- "Gamma=0.2"

boxplot(data_beta)
abline(h=0.0005,col="red")

#N* for selected method and setting
simul<-read.table("sim_fitting_25_500.csv",header=F,sep=',')
res<-data.frame(simul[,4])
names(res)[names(res)=="simul...4."]<-"Population"
p4<-ggplot(res,aes(y=Population,x=""))+
  geom_boxplot(outlier.size=2,lwd=1.2)+
  stat_summary(fun=mean,geom='point',color='red',size=3)+
  xlab("Smaller population") + ylab("Effective population size")+
  geom_hline(yintercept=500,color='blue',size=1)+
  theme(aspect.ratio = 1.9)
p4

grid.arrange(p1,p2,p3,p4,p5,p6,nrow=1)

#more sophisticated plotting for N*
simul<-read.table("complex_results\\complex_sim_n_9_1.csv",header=F,sep=',')
res<-data.frame(meth=factor(rep(c("Fitting to I","Fitting to R","Beta via R0","Fixed gamma","Combination",
                                  "Final size", "Imax"),each=1000)),
                pop=c(simul[,1],simul[,2],simul[,3],simul[,4],simul[,5],simul[,6],simul[,7]))
p<-ggplot(res,aes(y=pop,x=meth,fill=meth))+
  geom_boxplot(outlier.size=2,lwd=1)+
  stat_summary(fun=mean,geom='point',color='red',size=4)+
  ylab("Effective population size")+xlab("Method")+
  geom_hline(yintercept=9000,color='blue',size=1.5)+
  xlim("Fitting to I","Fitting to R","Beta via R0","Fixed gamma","Combination",
         "Final size", "Imax")+
  theme_light()+
  scale_fill_brewer(palette="Set2")+
  coord_cartesian(ylim=c(0, 18000))+
theme(text=element_text(size=20),legend.position="none")
p

#more sophisticated plotting for N*, NOTE not all methods present
simul<-read.table("complex_results\\complex_sim_n_3_16.csv",header=F,sep=',')
res<-data.frame(meth=factor(rep(c("Fitting to I","Beta via R0","Fixed gamma","Combination",
                                  "Final size", "Imax"),each=1000)),
                pop=c(simul[,1],simul[,3],simul[,4],simul[,5],simul[,6],simul[,7]))
p<-ggplot(res,aes(y=pop,x=meth,fill=meth))+
  geom_boxplot(outlier.size=2,lwd=1)+
  stat_summary(fun=mean,geom='point',color='red',size=4)+
  ylab("Effective population size")+xlab("Method")+
  geom_hline(yintercept=3000,color='blue',size=1.5)+
  xlim("Fitting to I","Beta via R0","Fixed gamma","Combination",
       "Final size", "Imax")+
  theme_light()+
  scale_fill_brewer(palette="Set2")+
  coord_cartesian(ylim=c(0, 6000))+
  theme(text=element_text(size=20),legend.position="none")
p

#more sophisticated for beta
simul<-read.table("complex_results\\complex_sim_beta_3_9.csv",header=F,sep=',')
res<-data.frame(meth=factor(rep(c("Fitting to I","Fitting to R","Fixed gamma"),each=1000)),
                pop=c(simul[,2],simul[,5],simul[,8]))
p<-ggplot(res,aes(y=pop,x=meth,fill=meth))+
  geom_boxplot(outlier.size=2,lwd=1)+
  stat_summary(fun=mean,geom='point',color='red',size=4)+
  ylab(expression(beta))+xlab("Method")+
  xlim("Fitting to I","Fitting to R","Fixed gamma")+
  theme_light()+
  scale_fill_brewer(palette="Set2")+
  theme(text=element_text(size=20),legend.position="none")
p

#more sophisticated for gamma
simul<-read.table("complex_results\\complex_sim_gamma_3_9.csv",header=F,sep=',')
res<-data.frame(meth=factor(rep(c("Fitting to I","Fitting to R","Beta via R0"),each=1000)),
                pop=c(simul[,2],simul[,5],simul[,8]))
p<-ggplot(res,aes(y=pop,x=meth,fill=meth))+
  geom_boxplot(outlier.size=2,lwd=1)+
  stat_summary(fun=mean,geom='point',color='red',size=4)+
  ylab(expression(gamma))+xlab("Method")+
  xlim("Fitting to I","Fitting to R","Beta via R0")+
  theme_light()+
  scale_fill_brewer(palette="Set2")+
  geom_hline(yintercept=0.2,color='blue',size=1.5)+
  theme(text=element_text(size=20),legend.position="none")
p

#more sophisticated for average SSE
simul<-read.table("complex_results\\complex_sim_sse_9_1.csv",header=F,sep=',')
res<-data.frame(meth=factor(rep(c("Fitting to I","Fitting to R","Beta via R0","Fixed gamma","Combination",
                                  "Final size", "Imax"),each=1000)),
                pop=c(simul[,1],simul[,2],simul[,3],simul[,4],simul[,5],simul[,6],simul[,7]))
p<-ggplot(res,aes(y=pop,x=meth,fill=meth))+
  geom_boxplot(outlier.size=2,lwd=1)+
  stat_summary(fun=mean,geom='point',color='red',size=4)+
  ylab("SSE over outbreak length")+xlab("Method")+
  xlim("Fitting to I","Fitting to R","Beta via R0","Fixed gamma","Combination",
       "Final size", "Imax")+
  theme_light()+
  scale_fill_brewer(palette="Set2")+
  scale_y_continuous(trans="log10")+
  theme(text=element_text(size=20),legend.position="none")
p

#more sophisticated for SSE if average is not computed
simul<-read.table("sse_fitting_15_1000.csv",header=F,sep=',')
lan<-read.table("length.csv",header=F,sep=',')
mlan<-lan[,2]
res<-data.frame(meth=factor(rep(c("Fitting to I","Fitting to R","Fixed gamma","Combination",
                                  "Final size","Imax"),each=1000)),
                pop=c(simul[,1]/mlan,simul[,4]/mlan,simul[,3]/mlan,
                      simul[,5]/mlan,simul[,6]/mlan,simul[,7]/mlan))
p<-ggplot(res,aes(y=pop,x=meth,fill=meth))+
  geom_boxplot(outlier.size=2,lwd=1)+
  stat_summary(fun=mean,geom='point',color='red',size=4)+
  ylab("SSE over outbreak length")+xlab("Method")+
  xlim("Fitting to I","Fitting to R","Fixed gamma","Combination",
       "Final size","Imax")+
  theme_light()+
  scale_y_continuous(trans="log10")+
  scale_fill_brewer(palette="Set2")+
  theme(text=element_text(size=20),legend.position="none")
p

#plotting selected cluster
mainarray<-largeo
t<-75
simul<-read.table("City_Confirmed_0115_0816_infected.csv",header=T,sep=',',fill=T)
df<-data.frame(city=factor(rep(paste(simul[mainarray,1],simul[mainarray,2],sep=" "),each=(t+1))),
               infec=c(t(as.vector(simul[mainarray,3:(t+3)]))),
               days=c(0:t))

p<-ggplot(df,aes(y=infec,x=days,group=city,color=city))+geom_line(size=1.25)+
  ylab("Number of currently infected")+xlab("Days since January 15th, 2021")+
  scale_x_continuous(limits=c(0,75),breaks=c(0,25,50,75))+
  theme_light()+
  theme(text=element_text(size=20))+
  scale_colour_discrete("City Province")
p

#beta*N for fixed gamma
df=data.frame(Location=locpop[,3],pop=locpop[,2],Method="Fixed gamma",InForce=betas[,8]*meann[,4])
df$Location<-factor(df$Location,levels=df$Location[order(df$pop)])

p<-ggplot(df,aes(y=InForce,x=Method,fill=Method))+
  geom_boxplot(outlier.size=2,lwd=1)+
  stat_summary(fun=mean,geom='point',color='red',size=4)+
  ylab(expression(beta*N^'*'))+xlab("Method")+
  theme_light()+
  scale_fill_brewer(palette="Set2")+
  theme(text=element_text(size=20),legend.position="none")
p
