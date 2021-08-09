library(ggplot2)
library(ggstance)
library(tidyverse)
library(corrplot)
#indeces for clustering
largeo<-c(35,36,38,39,40,42)
bigo<-c(9,32,33,34,37,41,52)
mediumo<-c(2,7,8,10:12,14,19,25:27,29,31,43,46)
mediumo2<-c(2,7,8,10,11,19,25,27,29,43,46)
smallo<-c(1,3:6,13,15:18,20:24,28,30,44,45,47:51,53)
smallo2<-smallo[c(1:6,8:9,11:12,14,15,17:21,24:25)]
#reading data
locpop<-read.table("Nresults\\china_pop.csv",header=F,sep=',')
pop<-read.table("Nresults\\china_pop2.csv",header=F, sep=",")
meann<-read.table("Nresults\\China_Ns.csv",header=F,sep=',')
lowern<-read.table("Nresults\\China_Ns_lower.csv",header=F,sep=',')
uppern<-read.table("Nresults\\China_Ns_upper.csv",header=F,sep=',')
betas<-read.table("Nresults\\China_Betas.csv", header=F,sep=',')

mainarray=bigo
#creating data frames
df1=data.frame(Location=locpop[mainarray,3],Mean_N=meann[mainarray,1],Method="Fitting to I",
              lowerCI=lowern[mainarray,1],upperCI=uppern[mainarray,1],pop=locpop[mainarray,2])
df2=data.frame(Location=locpop[mainarray,3],Mean_N=meann[mainarray,2],Method="Fitting to R",
              lowerCI=lowern[mainarray,2],upperCI=uppern[mainarray,2],pop=locpop[mainarray,2])

df2=data.frame(Location=locpop[mediumo2,3],Mean_N=meann[mediumo2,2],Method="Fitting to R",
               lowerCI=lowern[mediumo2,2],upperCI=uppern[mediumo2,2],pop=locpop[mediumo2,2])
df2=data.frame(Location=locpop[smallo2,3],Mean_N=meann[smallo2,2],Method="Fitting to R",
               lowerCI=lowern[smallo2,2],upperCI=uppern[smallo2,2],pop=locpop[smallo2,2])

df3=data.frame(Location=locpop[mainarray,3],Mean_N=meann[mainarray,3],Method="Beta via R0",
               lowerCI=lowern[mainarray,3],upperCI=uppern[mainarray,3],pop=locpop[mainarray,2])
df4=data.frame(Location=locpop[mainarray,3],Mean_N=meann[mainarray,4],Method="Fixed gamma",
               lowerCI=lowern[mainarray,4],upperCI=uppern[mainarray,4],pop=locpop[mainarray,2])
df5=data.frame(Location=locpop[mainarray,3],Mean_N=meann[mainarray,5],Method="Combination",
               lowerCI=lowern[mainarray,5],upperCI=uppern[mainarray,5],pop=locpop[mainarray,2])
df6=data.frame(Location=locpop[mainarray,3],Mean_N=meann[mainarray,6],Method="Imax",
               lowerCI=lowern[mainarray,6],upperCI=uppern[mainarray,6],pop=locpop[mainarray,2])
df7=data.frame(Location=locpop[mainarray,3],Mean_N=meann[mainarray,7],Method="Final size",
               lowerCI=lowern[mainarray,7],upperCI=uppern[mainarray,7],pop=locpop[mainarray,2])
#combining and factoring
df<-rbind(df1,df2,df3,df4,df5)
df<-df4
df1$Location<-factor(df1$Location,levels=df1$Location[order(df1$pop)])

#plotting, change labels and breaks appropriately for the selected cluster
p<-ggplot(df, aes(y=Location, x=Mean_N, group=Method, color=Method)) +
  xlab(expression(N^'*')) + ylab("Location (N)") +
  geom_errorbarh(mapping=aes(xmin=lowerCI, xmax=upperCI), size = 1, alpha=0.99, position=position_dodgev(height=0.8), height=0) +
  geom_point(position=position_dodgev(height=0.8),size=2)  +
  guides(colour = guide_legend(nrow = 2))+ scale_x_continuous( limits=c(0,500), 
  labels=c("0","100","200","300","400","500"), breaks=c(0,100,200,300,400,500))+
  theme_light()+ 
  theme(legend.position="bottom",text=element_text(size=20))
p

#correlations between ratios and distances

dg1=data.frame(Mean_N=meann[largeo,4],pop_n=pop[largeo,2],dist=pop[largeo,3],outbreak_size="Large outbreak")
dg2=data.frame(Mean_N=meann[bigo,4],pop_n=pop[bigo,2],dist=pop[bigo,3],outbreak_size="Big outbreak")
dg3=data.frame(Mean_N=meann[mediumo,4],pop_n=pop[mediumo,2],dist=pop[mediumo,3],outbreak_size="Medium outbreak")
dg4=data.frame(Mean_N=meann[smallo,4],pop_n=pop[smallo,2],dist=pop[smallo,3],outbreak_size="Small outbreak")
#plotting those correlations
dg<-rbind(dg1,dg2,dg3,dg4)
dg$outbreak_size<-factor(dg$outbreak_size,levels=c("Large outbreak","Big outbreak","Medium outbreak","Small outbreak"))
q<-ggplot(dg,aes(y=Mean_N/pop_n*10^-6,x=dist,group=outbreak_size,color=outbreak_size))+
  geom_point()+
  xlab("Distance from Wuhan, km")+ylab(expression(N^'*'/N))+
  theme(legend.position = "bottom",text=element_text(size=15))
q


#correlations between N*s from different methods
meann<-read.table("Nresults\\China_Ns.csv",header=F,sep=',')
names(meann)[1]<-"F to I"
names(meann)[2]<-"F to R"
names(meann)[3]<-"B via R0"
names(meann)[4]<-"Fixed g"
names(meann)[5]<-"Comb"
names(meann)[6]<-"Imax"
names(meann)[7]<-"Final"

N<-cor(meann)
corrplot.mixed(N)

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
