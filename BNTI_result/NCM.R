rm(list = ls())
setwd("E:/博士/合作/周欣/bnti")
getwd()
tw<-t(read.csv("wh.csv",row.names = 1))
write.csv(tw,"tw.csv")
tw<-read.csv("tw.csv")
names(tw)[1]<-c("SampleID")
group<-read.csv("Group.csv")
wh<-merge(group,tw,by="SampleID")


sp<-wh[which(wh$Seasonal=="Winter"),]
sp<-sp[-1:-11]
spp<-sp
library(Hmisc)
library(minpack.lm)
library(stats4)
source("sncm.fit.R")

N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
##Fit model parameter m (or Nm) using Non-linear least squares (NLS) ??? (NLS)拟合模型参数 m (or Nm)
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr
N*coef(m.fit)
N
ncol(spp)
m<-coef(m.fit)
m

#再输???
#write.csv(p, file = "j:/iue/aehg11/neutr_model/p.csv")
#write.csv(freq, file = "j:/iue/aehg11/neutr_model/freq.csv")
#write.csv(freq.pred, file = "j:/iue/aehg11/neutr_model/freq.pred.csv")

#得到3个csv文件，把这三个文件里面的内容拷贝出来，如附件所示，再输入以下命令拟合???
#argnls=read.table("neutrl_file.txt",row.names=1,header=T)
bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-"#e31a1c"
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-"#1f78b4"
library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='orangered3',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='#00B76D',lwd=2,lty=2),default='native')  ###左边置信区间
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='#00B76D',lwd=2,lty=2),default='native')  ###右边置信区间
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) #添加横坐标log10(p)
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) #添加纵坐标Frequency
#grid.text(x=unit(0,'npc')-unit(-1,'lines'), y=unit(0,'npc')-unit(-15,'lines'),label='Mean Relative Abundance (log)', gp=gpar(fontface=2)) #添加横坐标log10(p)
#grid.text(round(coef(m.fit)*N),x=unit(0,'npc')-unit(-5,'lines'), y=unit(0,'npc')-unit(-15,'lines'),gp=gpar(fontface=2)) #添加横坐标log10(p)
#grid.text(label = "Nm=",x=unit(0,'npc')-unit(-3,'lines'), y=unit(0,'npc')-unit(-15,'lines'),gp=gpar(fontface=2)) #添加横坐标log10(p)
#grid.text(round(Rsqr,2),x=unit(0,'npc')-unit(-5,'lines'), y=unit(0,'npc')-unit(-16,'lines'),gp=gpar(fontface=2)) #添加横坐标log10(p)
#grid.text(label = "Rsqr=",x=unit(0,'npc')-unit(-3,'lines'), y=unit(0,'npc')-unit(-16,'lines'),gp=gpar(fontface=2)) #添加横坐标log10(p)
draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,4),"\n","m=",round(m,4)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)

