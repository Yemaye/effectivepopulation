library(tidyverse)
library(ggmap)
library(forcats)
library(maps)
library(maptools)

#map of whole China
chinamap<-map_data("world", "China")
#China cities data
chinacities<-world.cities[world.cities$country.etc=="China",] 
#Getting only necessary info
chinac<-select(chinacities,name,long,lat)
#Manual searching necessary cities. Each greek letter stands for the size of outbreak
alphas<-c("Xiangfan","Ezhou","Xiaogan","Jingzhou","Huanggang","Suizhou")
betas<-c("Wenzhou","Huangshi","Shiyan","Yichang","Jingmen","Xianning","Chongqing")
gammas<-c("Harbin","Hangzhou","Ningbo","Taizhou","Hefei","Bengbu","Fuyang","Nanchang","Shangrao","Zhengzhou",
            "Nanyang","Xinyang","Zhumadian","Changsha","Yueyang")
deltas<-c("Tianjin","Nanjing","Xuzhou","Huaiyin","Anqing","Luan","Bozhou","Putian","Jiujiang",
          "Xinyu","Ganzhou","Yichun","Shangqiu","Zhoukou","Zhuzhou","Shaoyang","Changde","Loudi",
          "Zhuhai","Huizhou","Dongguan","Xian")
#converting to dataframe
df_alpha <- data.frame(chinac[chinac$name%in%alphas,])
df_beta <- data.frame(chinac[chinac$name%in%betas,])
df_gamma <- data.frame(chinac[chinac$name%in%gammas,])
df_delta <- data.frame(chinac[chinac$name%in%deltas,])
#we still missing or having duplicates
df_gamma<-df_gamma[c(1:3,5:16),]
df_delta<-df_delta[c(1:19,21:23),]
df_delta<-rbind(df_delta,data.frame(name = "Fuzhou", long = 116.21,lat=27.56))
df_delta<-rbind(df_delta,data.frame(name = "Fuzhou", long=119.17, lat=27.04))
df_delta<-rbind(df_delta,data.frame(name="Suzhou", long=120.35, lat=31.17))
#plotting
ggplot() +
geom_polygon(data = chinamap, aes(x = long, y = lat, group = group))+
geom_point(data = df_alpha, aes(x = long, y = lat), color = 'red',size=4)+
geom_point(data = df_beta, aes(x = long, y = lat), color = 'orange',size=3)+
geom_point(data = df_gamma, aes(x = long, y = lat), color = 'yellow',size=2)+
geom_point(data = df_delta, aes(x = long, y = lat), color = 'white',size=1.5)