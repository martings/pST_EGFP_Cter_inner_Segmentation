library(ggplot2)

TM_sinCero<-subset(TM, cell_msk!=0)
short_sinCero<-subset(short, cell_msk!=0)

short_ver_sinCero<-subset(short_ver, cell_msk!=0)
TM_ver_sinCero<-subset(TM_ver, cell_msk!=0)





TM_sinCero[,4]<-"TM"
short_sinCero[,4]<-"Short"
todo<-rbind(TM_sinCero,short_sinCero)

p<-ggplot(data= todo, aes(x=X, y=log10(EGFP),colour=V4))+geom_point()+geom_smooth()+xlim(10,100)
p

p<-ggplot(data= todo, aes(x=X, y=log10(EGFP),colour=V4))+geom_smooth()+xlim(10,100)
p

p<-ggplot(data= todo, aes(x=X, y=EGFP,colour=V4))+geom_smooth()+xlim(10,100)
p

p<-ggplot(data= todo, aes(x=V4,y=EGFP))+geom_boxplot()
p

p<-ggplot(data= todo, aes(x=EGFP, colour=V4, fill= V4))+geom_histogram(bins = 50, alpha=0.5, position="identity")
p

p<-ggplot(data= todo, aes(x=EGFP, colour=V4, fill= V4))+labs(title= "Horizontal")+geom_histogram(bins = 50, alpha=0.5, position="identity")
p


t.test(EGFP~V4, data = todo)
wilcox.test(EGFP~V4, data = todo)


TM_ver_sinCero[,4]<-"TM"
short_ver_sinCero[,4]<-"Short"
todo_ver<-rbind(TM_ver_sinCero,short_ver_sinCero)

p<-ggplot(data= todo_ver, aes(x=EGFP, colour=V4, fill= V4))+labs(title= "Vertical")+geom_histogram(bins = 50, alpha=0.5, position="identity")
p

todo[,5]<-"horizontal"
todo_ver[,5]<-"vertical"

todo_todo<-rbind(todo,todo_ver)

p<-ggplot(data= todo_todo, aes(x=EGFP, colour=V4, fill= V4))+labs(title= "Todo")+geom_histogram(bins = 50, alpha=0.5, position="identity")
p

p<-ggplot(data= todo_todo, aes(x=EGFP, colour=V4, fill= V4))+labs(title= "Todo")+geom_density(bins = 50, alpha=0.5, position="identity")
p


p<-ggplot(data= todo_todo, aes(x=EGFP, y=X, colour=V4, fill= V4))+labs(title= "Todo")+geom_density2d(bins = 50, alpha=0.5, position="identity")
p


p<-ggplot(data= todo_todo, aes(x=EGFP, colour=V4, fill= V4))+labs(title= "Todo")+geom_freqpoly(bins = 50, alpha=0.5, position="identity")
p


p<-ggplot(data= todo_todo, aes(x=EGFP, colour=V4, fill= V5))+labs(title= "Todo")+geom_histogram(bins = 50, alpha=0.5, position="identity")
p


t.test(EGFP~V4, data = todo_todo)
wilcox.test(EGFP~V4, data = todo_todo)

#Por partes

p<-ggplot(data= todo_ver, aes(x=EGFP, colour=V4, fill= V4))+labs(title= "Vertical")+geom_density(bins = 50, alpha=0.5, position="identity")
p

p<-ggplot(data= todo, aes(x=EGFP, colour=V4, fill= V4))+labs(title= "Horizontal")+geom_density(bins = 50, alpha=0.5, position="identity")
p

