###library
library(ggplot2)
library(data.table)
library(plyr)
#######Cargo cvs verticales
setwd("/media/maxpower/Datos/Martin/Laboratorio/Microscopio/2019-06-17/citomCherry/short/stack/buenas/procesadas_short")
temp = list.files(pattern="*tifhorizontal.csv")
for (i in 1:length(temp)) assign(temp[i], read.csv(temp[i] ,header = TRUE))


setwd("/media/maxpower/Datos/Martin/Laboratorio/Microscopio/2019-06-17/citomCherry/TM/procesadas/buenas/procesadas/")


temp<-as.data.frame(temp)
temp$temp <- paste0("`", temp$temp, "`")
temp<-list(temp)
short<-do.call(rbind, list(temp))
short<-ldply(list(temp), rbind)


TM<-rbind(`10MAX_TM3.tifhorizontal.csv`, `12MAX_TM3.tifhorizontal.csv`, `15MAX_TM3.tifhorizontal.csv`, `16MAX_TM3.tifhorizontal.csv`,`1MAX_TM3.tifhorizontal.csv`,`4MAX_TM3.tifhorizontal.csv`, `6MAX_TM3.tifhorizontal.csv`, `7MAX_TM3.tifhorizontal.csv`)

short<-rbind(`10MAX_short2.tifhorizontal.csv`,`11MAX_short2.tifhorizontal.csv`,`12MAX_short2.tifhorizontal.csv`,`14MAX_short2.tifhorizontal.csv`,
      `15MAX_short2.tifhorizontal.csv`,`17MAX_short2.tifhorizontal.csv`,`19MAX_short2.tifhorizontal.csv`,`21MAX_short2.tifhorizontal.csv`,
      `22MAX_short2.tifhorizontal.csv`,`24MAX_short2.tifhorizontal.csv`,`25MAX_short2.tifhorizontal.csv`,`26MAX_short2.tifhorizontal.csv`,
      `28MAX_short2.tifhorizontal.csv`,`29MAX_short2.tifhorizontal.csv`,`30MAX_short2.tifhorizontal.csv`,`31MAX_short2.tifhorizontal.csv`,
      `32MAX_short2.tifhorizontal.csv`,`33MAX_short2.tifhorizontal.csv`,`34MAX_short2.tifhorizontal.csv`,`4MAX_short2.tifhorizontal.csv`,
      `6MAX_short2.tifhorizontal.csv`,`7MAX_short2.tifhorizontal.csv`,`8MAX_short2.tifhorizontal.csv`)

p<-ggplot(data= short, aes(x=X, y=EGFP,colour=factor(cell_msk)))+geom_point()+geom_smooth(method="auto", se=T)+
  facet_wrap(~ cell_msk)
p


p<-ggplot(data= short, aes(x=X, y=EGFP,colour=factor(cell_msk)))+geom_point()+geom_smooth()
p
short_sinCero<-subset(short, cell_msk!=0)

p<-ggplot(data= short_sinCero, aes(x=factor(cell_msk), y= log2(EGFP)))+geom_boxplot()
p
p<-ggplot(data= short_sinCero, aes(x=X, y=EGFP))+geom_point()+geom_smooth()
p
p<-ggplot(data= short_sinCero, aes(x=X, y=(EGFP),colour=factor(cell_msk)))+geom_point()+geom_smooth(se=F)
p

#funciona
p<-ggplot(data= short_sinCero, aes(x=X, y=log10(EGFP)))+geom_smooth(se=F)+xlim(0,100)+scale_y_log10()#+geom_point()
p

p<-ggplot(data= short_sinCero, aes(x=EGFP, ,colour=factor(cell_msk)))+geom_histogram()
p
summary(short)

######Cargo cvs horizontales
setwd("/media/maxpower/Datos/Martin/Laboratorio/Microscopio/2019-06-17/citomCherry/short/stack/buenas/procesadas_short")
temp = list.files(pattern="*tifvertical.csv")
for (i in 1:length(temp)) assign(temp[i], read.csv(temp[i] ,header = TRUE))

short_ver<-rbind(`10MAX_short2.tifvertical.csv`,`11MAX_short2.tifvertical.csv`,`12MAX_short2.tifvertical.csv`,`14MAX_short2.tifvertical.csv`,
             `15MAX_short2.tifvertical.csv`,`17MAX_short2.tifvertical.csv`,`19MAX_short2.tifvertical.csv`,`21MAX_short2.tifvertical.csv`,
             `22MAX_short2.tifvertical.csv`,`24MAX_short2.tifvertical.csv`,`25MAX_short2.tifvertical.csv`,`26MAX_short2.tifvertical.csv`,
             `28MAX_short2.tifvertical.csv`,`29MAX_short2.tifvertical.csv`,`30MAX_short2.tifvertical.csv`,`31MAX_short2.tifvertical.csv`,
             `32MAX_short2.tifvertical.csv`,`33MAX_short2.tifvertical.csv`,`34MAX_short2.tifvertical.csv`,`4MAX_short2.tifvertical.csv`,
             `6MAX_short2.tifvertical.csv`,`7MAX_short2.tifvertical.csv`,`8MAX_short2.tifvertical.csv`)

setwd("/media/maxpower/Datos/Martin/Laboratorio/Microscopio/2019-06-17/citomCherry/TM/procesadas/buenas/procesadas/")
temp = list.files(pattern="*tifvertical.csv")
for (i in 1:length(temp)) assign(temp[i], read.csv(temp[i] ,header = TRUE))


TM_ver<-rbind(`10MAX_TM3.tifvertical.csv`, `12MAX_TM3.tifvertical.csv`, `15MAX_TM3.tifvertical.csv`, `16MAX_TM3.tifvertical.csv`,`1MAX_TM3.tifvertical.csv`,`4MAX_TM3.tifvertical.csv`, `6MAX_TM3.tifvertical.csv`, `7MAX_TM3.tifvertical.csv`)


p<-ggplot(data= TM_ver, aes(x=X, y=EGFP,colour=factor(cell_msk)))+geom_point()+geom_smooth(method="auto", se=T)+
  facet_wrap(~ cell_msk)
p

short_ver_sinCero<-subset(short_ver, cell_msk!=0)

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
