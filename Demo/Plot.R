#以最大值最小值为顶的箱图
x <- matrix(c(36,97,33,89,45,99,51,93,47,88),2,5)
boxplot(x,medlty="blank",names=c("A","B","C","D","E"),col="pink", boxwex=0.35)
abline(h=71,col="navy", lwd=2, lty=5)

#表示数据在最大最小之间位置的线型图
X <- matrix(c(36,88,97,33,86,89,45,77,99,51,90,93,47,65,88),3,5)
plot(c(X[1,],X[3,],X[2,]),c(Y,Y,Y),pch = c(rep(19,10),rep(4,5)),cex = 1.5,col = c(rep("seagreen",10),rep("magenta",5)),lwd = 2,xlab = "成绩",ylab= "科目",yaxt = "n")
axis(2,at=c(1:5),labels=c('A','B','C','D','E'))
arrows(c(X[2,],X[2,]),c(Y,Y),c(X[1,],X[3,]), c(Y,Y),col = "springgreen",lwd = 2,length = 0.15,angle = 20)

#Polar barChart:
library(ggplot2)
GiniData<-read.csv('IncomeInequality.csv',header = T)
Gini<-ggplot(GiniData,aes(x=paste(GiniIndex,Country),y=GiniIndex,fill=GiniIndex%/%10))
Gini<-ggplot(GiniData,aes(x=paste(GiniIndex,Country),y=GiniIndex,fill=(sign(GiniIndex-41.60)+sign(Country=="Coted'lvoire")*2)*sign(Country!="China")))
Gini<- Gini + geom_bar(stat="identity",position="dodge")+coord_polar()
Gini<- Gini + scale_fill_continuous(high="darkred",low="darkgreen")
Gini<- Gini + theme(
panel.background=element_rect(fill="white",colour = "white", size=0),
axis.text=element_blank(),
axis.title=element_blank(),
legend.title=element_blank())
x <-c(1:dim(GiniData)[1])
Gini +geom_text(
aes(
x=x,
label=paste(GiniData$GiniIndex,GiniData$Country),
angle=270-x/134*360,
hjust=1),
y=GiniData$GiniIndex+3,
size=3,
vjust=0)
