library(ggplot2)
library(ggstance)
library(tidyverse)
library(corrplot)
library(gridExtra)
library(devtools)
library(lemon)
library(reshape)
library(ggpubr)

#Reading the data
locpop<-read.table("Nresults\\china_pop.csv",header=F,sep=',') #Contains city information on a census pop. size (CPS) and a combined format
pop<-read.table("Nresults\\china_pop2.csv",header=F, sep=",") #Contains city information on CPS, distance to Wuhan, I_max and R_infinity
meann<-read.table("Nresults\\China_Ns.csv",header=F,sep=',') #N effective obtained from Fixed gamma method
lowern<-read.table("Nresults\\China_Ns_lower.csv",header=F,sep=',') #lower CI
uppern<-read.table("Nresults\\China_Ns_upper.csv",header=F,sep=',') #Uppes CI
betas<-read.table("Nresults\\China_Betas.csv", header=F,sep=',') #beta obtained from Fixed gamma method
trusol<-read.table("true_sol.csv",header=F,sep=",") #True solution for given settings

#Indices of cities clustered according to Imax
largeo<-c(35,36,38,39,40,42) 
bigo<-c(9,32,33,34,37,41,52)
mediumo<-c(2,7,8,10:12,14,19,25:27,29,31,43,46)
smallo<-c(1,3:6,13,15:18,20:24,28,30,44,45,47:51,53)

#Figure 1, plots simulated data alongside with a true solution for the given setting
infecdata<-read.table("infec_1000_35.csv",header=F,sep=",") #reading the data
nc=50 #Limit your data here
df<-data.frame(t(infecdata[,1:nc])) #Modify so it could be usef for ggplot
df<-melt(df)
dg<-bind_rows(rep(list(data.frame(t=c(0:(nc-1)))),1000))
df<-data.frame(df,dg)

trusol<-read.table("true_sol.csv",header=F,sep=",") #Reading true solution provided by .jl code, 
dt<-data.frame(trusol[6,1:nc]) #change the first index depending on the index
dt<-melt(dt) #time data
dt<-data.frame(dt,data.frame(t=c(0:(nc-1))))

#plotting, change label and some parameters appropriately
p<-ggplot()+
  geom_line(df,mapping=aes(x=t,y=value,color=variable),size=1)+
  geom_line(dt,mapping=aes(x=t,y=value),color='black',size=3)+
  theme_light()+
  xlab("Days since introduction")+ylab("Active infections")+
  xlim(0,50)+
  labs(caption=expression("f) N = 1000, "~beta~"= 0.0007, "~gamma~"= 0.2"))+
  #scale_color_hue(l=60,c=80)+
  theme(legend.position = "none",axis.text=element_text(size=14),
        axis.title=element_text(size=15),plot.caption = element_text(size=16,hjust=0.5,vjust=-1))
p

#Figure 3, produces boxplots for Neffective for a given setting
#Change the file according to the setting. The file should contain Neff computed by the Fixed gamma method
fixedg<-read.table("sim_fitting_35_1000.csv",header=F,sep=',')
res<-data.frame(fixedg[,4])
names(res)[names(res)=="fixedg...4."]<-"Population"
#Removing outliers
q1<-unname(quantile(res$Population,1/4))
q2<-unname(quantile(res$Population,3/4))
Iq<-IQR(res$Population)
rowsneed<-res$Population>=(q1-1.5*Iq)&res$Population<=(q2+1.5*Iq)
res1<-subset(res,rowsneed)
#Plotting, change yintercept and labels appropriately
p<-ggplot(res1,aes(y=Population,x=""))+
  geom_boxplot(outlier.size=2,lwd=1.2)+
  stat_summary(fun=mean,geom='point',color='red',size=4)+
  xlab(expression("")) + ylab("Effective population size")+
  theme_light()+
  geom_hline(yintercept=1000,color='blue',size=2)+
  labs(caption=expression("f) N = 1000, "~beta~"= 0.0007, "~gamma~"= 0.2"))+
  #coord_cartesian(ylim=c((q1-1.5*Iq),(q2+1.5*Iq)))+
  theme(aspect.ratio = 1.9,legend.position = "none")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),plot.caption = element_text(size=16,hjust=0.5))
p

#Figures 4 and 5, produces Neffective (lower CI, optimum, upper CI) for a given city and allows to compare with others
#note that it would be difficult to represent all the clusters simultaneously, so we have to group
df1=data.frame(Location=locpop[largeo,3],Mean_N=meann[largeo,4],Method="Fixed gamma",
               lowerCI=lowern[largeo,4],upperCI=uppern[largeo,4], pop=locpop[largeo,2],sizeo="largeo")
df2=data.frame(Location=locpop[bigo,3],Mean_N=meann[bigo,4],Method="Fixed gamma",
               pop=locpop[bigo,2],sizeo="bigo",lowerCI=lowern[bigo,4],upperCI=uppern[bigo,4])
df3=data.frame(Location=locpop[mediumo,3],Mean_N=meann[mediumo,4],Method="Fixed gamma",
               pop=locpop[mediumo,2],sizeo="mediumo",lowerCI=lowern[mediumo,4],upperCI=uppern[mediumo,4])
df4=data.frame(Location=locpop[smallo,3],Mean_N=meann[smallo,4],Method="Fixed gamma",
               pop=locpop[smallo,2],sizeo="smallo",lowerCI=lowern[smallo,4],upperCI=uppern[smallo,4])
df<-rbind(df1,df2,df3) #Choose this for 3 clusters with Imax>=100
df<-df4 #Choose this for the cluster with Imax<100

#Change labels and breaks appropriately for the selected cluster
p<-ggplot(df, aes(y=Location, x=Mean_N, group=Method, color=sizeo)) +
  xlab(expression(N^'*')) + ylab("Location (N, millions)") +
  geom_errorbarh(mapping=aes(xmin=lowerCI, xmax=upperCI), size = 1, alpha=0.99, position=position_dodgev(height=0.8), height=0) +
  geom_point(position=position_dodgev(height=0.8),size=2)  +
  guides(colour = guide_legend(nrow = 2))+ 
  theme_light()+ 
  #scale_x_log10()+
  guides(color=guide_legend(("Maximum active cases")))+
  scale_color_manual(breaks=c("largeo","bigo","mediumo"),labels = c("1000+","250-999","100-249"), values = c("red","orange","yellow"))+
  #scale_color_manual(labels = c("50-99"), values = c("green"))+
  theme(legend.position="bottom",text=element_text(size=14))
p

#Figure 6, boxplot for the infection force beta*N from all cities. Adds clustering
df1=data.frame(Location=locpop[largeo,3],Mean_N=meann[largeo,4],Method="Fixed gamma",
            pop=locpop[largeo,2],sizeo="largeo",InForce=betas[largeo,8]*meann[largeo,4],bet=betas[largeo,8])
df2=data.frame(Location=locpop[bigo,3],Mean_N=meann[bigo,4],Method="Fixed gamma",
               pop=locpop[bigo,2],sizeo="bigo",InForce=betas[bigo,8]*meann[bigo,4],bet=betas[bigo,8])
df3=data.frame(Location=locpop[mediumo,3],Mean_N=meann[mediumo,4],Method="Fixed gamma",
               pop=locpop[mediumo,2],sizeo="mediumo",InForce=betas[mediumo,8]*meann[mediumo,4],bet=betas[mediumo,8])
df4=data.frame(Location=locpop[smallo,3],Mean_N=meann[smallo,4],Method="Fixed gamma",
               pop=locpop[smallo,2],sizeo="smallo",InForce=betas[smallo,8]*meann[smallo,4],bet=betas[smallo,8])
df<-rbind(df1,df2,df3,df4)
df$Location<-factor(df$Location,levels=df$Location[order(df$pop)])
#Plotting, can change to see beta alone
p<-ggplot(df,aes(y=bet*Mean_N,x=Method,fill=Method))+
  geom_boxplot(outlier.size=2,lwd=1,fill="grey")+
  geom_point(position=position_jitterdodge(),size=4,aes(color=sizeo))+
  stat_summary(fun=mean,geom='point',color='blue',size=6)+
  ylab(expression(beta))+xlab("")+
  theme_light()+
  scale_color_manual(breaks=c("largeo","bigo","mediumo","smallo"),
                     labels = c("1000+","250-999","100-249","50-99"), 
                     values = c("red","orange","yellow","green"))+
  guides(color=guide_legend(("Maximum active cases"),nrow=2,byrow=TRUE),fill="none")+
  theme(text=element_text(size=15),legend.position="bottom",axis.text.x=element_blank())
p

#Figure 7 and S1, correlation plot for N*/N and distance to Wuhan (and N^* and final size)
#need slightly a different database, but f and g are simular and could be combined
dg1=data.frame(Mean_N=meann[largeo,4],pop_n=pop[largeo,2],dist=pop[largeo,3],finsize=pop[largeo,5],outbreak_size="Large outbreak")
dg2=data.frame(Mean_N=meann[bigo,4],pop_n=pop[bigo,2],dist=pop[bigo,3],finsize=pop[bigo,5],outbreak_size="Big outbreak")
dg3=data.frame(Mean_N=meann[mediumo,4],pop_n=pop[mediumo,2],dist=pop[mediumo,3],finsize=pop[mediumo,5],outbreak_size="Medium outbreak")
dg4=data.frame(Mean_N=meann[smallo,4],pop_n=pop[smallo,2],dist=pop[smallo,3],finsize=pop[smallo,5],outbreak_size="Small outbreak")
#Plotting those correlations
dg<-rbind(dg1,dg2,dg3,dg4)
dg$outbreak_size<-factor(dg$outbreak_size,levels=c("Large outbreak","Big outbreak","Medium outbreak","Small outbreak"))
q<-ggplot(dg,aes(y=Mean_N/pop_n*10^-6,x=dist,group=outbreak_size,color=outbreak_size))+
#q<-ggplot(dg,aes(y=Mean_N,x=finsize,group=outbreak_size,color=outbreak_size))+ #Change to this if you want to see other correlation of N and fin size
  geom_point(size=2)+
  xlab("Distance from Wuhan, km")+ylab(expression(N^'*'/N))+
  #xlab("Final size")+ylab(expression(N^'*'))+ #Change to this for finsize/N^* correlation
  theme_light()+
  #stat_summary(fun=mean,geom='point',color='blue',size=3)+ #Adds mean
  scale_color_manual(breaks=c("Large outbreak","Big outbreak","Medium outbreak","Small outbreak"),
                     labels = c("1000+","250-999","100-249","50-99"), values = c("red","orange","yellow","green"))+
  guides(color=guide_legend(("Maximum active cases")))+
  theme(legend.position = "bottom",text=element_text(size=15))
q

#Figure S2, creates a boxplots to compare N^* produced by fixed gamma (LSQ) and survival analysis (SDS) methods
sellke<-read.table("Sellke_N.csv",header=F,sep=',') #reading N's for a given setting
fixedg<-read.table("Fixedg_N.csv",header=F,sep=',') #Our file contained N's for all settings
#res<-data.frame(sellke[,1],fixedg[,1])
i=2 #Fix the setting
#Change the last part so that it matches with the setting and the number of data's
res<-data.frame(meth=factor(rep(c("SDS","LSQ"),each=20)),
                pop=abs(c(sellke[,i],fixedg[,i])-1000)/10)
#Ploting, change labels and parameters accordingly
p<-ggplot(res,aes(y=pop,x=meth,fill=meth))+
  geom_boxplot(outlier.size=2,lwd=1.2)+
  stat_summary(fun=mean,geom='point',color='red',size=3)+
  xlab(expression("e) N = 1000, "~beta~"= 0.0003, "~gamma~"= 0.2"))+
  #xlab(expression(atop("d) N = 1000, "~beta~"= 0.0003,",paste(~gamma~"= 0.2"))))+
  ylab("Relative error (%)")+
  coord_cartesian(ylim=c(0,50))+
  theme(aspect.ratio = 1.9,legend.position = "none",axis.text=element_text(size=14),
        axis.title.x=element_text(size=15,hjust=0.5,vjust=-1), axis.title.y=element_text(size=15))
p

