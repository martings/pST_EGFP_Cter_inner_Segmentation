###library
library(ggplot2)

#######Cargo cvs verticales
setwd("/media/maxpower/Datos/Martin/Laboratorio/Microscopio/2019-06-17/citomCherry/TM/procesadas/buenas/procesadas")
temp = list.files(pattern="*tifhorizontal.csv")
for (i in 1:length(temp)) assign(temp[i], read.csv(temp[i] ,header = TRUE))


p<-ggplot(data=`4MAX_TM3.tifhorizontal.csv`, aes(x=X, y=EGFP,colour=cell_msk))+geom_point()+geom_smooth(method="auto", se=T)
p

TM<-rbind(`10MAX_TM3.tifhorizontal.csv`, `12MAX_TM3.tifhorizontal.csv`, `15MAX_TM3.tifhorizontal.csv`, `16MAX_TM3.tifhorizontal.csv`,`1MAX_TM3.tifhorizontal.csv`,`4MAX_TM3.tifhorizontal.csv`, `6MAX_TM3.tifhorizontal.csv`, `7MAX_TM3.tifhorizontal.csv`)

p<-ggplot(data= TM, aes(x=X, y=EGFP,colour=factor(cell_msk)))+geom_point()+geom_smooth(method="auto", se=T)+
  facet_wrap(~ cell_msk)
p


p<-ggplot(data= TM, aes(x=X, y=EGFP,colour=factor(cell_msk)))+geom_point()+geom_line()
p

p<-ggplot(data= TM, aes(x=factor(cell_msk), y= log2(EGFP)))+geom_boxplot()

p
TM_sinCero<-subset(TM, cell_msk!=0)

p<-ggplot(data= TM_sinCero, aes(x=X, y=EGFP))+geom_point()+geom_smooth()
p


p<-ggplot(data= TM_sinCero, aes(x=X, y=EGFP,colour=factor(cell_msk)))+geom_point()+geom_smooth()
p

######Cargo cvs horizontales
setwd("/media/maxpower/Datos/Martin/Laboratorio/Microscopio/2019-06-17/citomCherry/TM/procesadas/buenas/procesadas")
temp = list.files(pattern="*.csv")
for (i in 1:length(temp)) assign(temp[i], read.csv(temp[i] ,header = TRUE))

TM_ver<-rbind(`10MAX_TM3.tifvertical.csv`, `12MAX_TM3.tifvertical.csv`, `15MAX_TM3.tifvertical.csv`, `16MAX_TM3.tifvertical.csv`,`1MAX_TM3.tifvertical.csv`,`4MAX_TM3.tifvertical.csv`, `6MAX_TM3.tifvertical.csv`, `7MAX_TM3.tifvertical.csv`)

p<-ggplot(data= TM_ver, aes(x=X, y=EGFP,colour=factor(cell_msk)))+geom_point()+geom_smooth(method="auto", se=T)+
  facet_wrap(~ cell_msk)
p

TM_ver_sinCero<-subset(TM_ver, cell_msk!=0)

p<-ggplot(data= TM_ver_sinCero, aes(x=X, y=log2(EGFP),colour=factor(cell_msk)))+geom_point()+geom_smooth(method="auto", se=T)+
  facet_wrap(~ cell_msk)
p

p<-ggplot(data= TM_ver_sinCero, aes(x=X, y=log2(EGFP)))+geom_point()+geom_smooth(method="auto", se=T)
p

TM_sinCero[,4]<-"Horizontal"
TM_ver_sinCero[,4]<-"Vertical"
TM_todo<-rbind(TM_sinCero,TM_ver_sinCero)

p<-ggplot(data= TM_todo, aes(x=X, y=log2(EGFP)))+geom_point()+geom_smooth(method="auto", se=T)
p

p<-ggplot(data= TM_todo, aes(x=X, y=log2(EGFP)))+geom_point()+geom_smooth(method="auto", se=T)+
  facet_grid(~V4)
p

p<-ggplot(data= TM_todo, aes(x=X, y=log2(EGFP)))+geom_point()+geom_smooth(method="auto", se=T)+
  facet_grid(~cell_msk)
p
