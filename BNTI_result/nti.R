rm(list = ls())
setwd("E:/博士/合作/周欣/bnti")
getwd()
library(ggplot2)
library(reshape2)
library(vegan)
library(ape)
library(cowplot)
Beta_NTI<-function(phylo,comun,beta.reps=99){
  require(picante)
  
  comun=t(comun)
  match.phylo.comun = match.phylo.data(phylo, t(comun))
  beta.mntd.weighted = as.matrix(comdistnt(t(match.phylo.comun$data),cophenetic(match.phylo.comun$phy),abundance.weighted=T))
  
  rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(match.phylo.comun$data),ncol(match.phylo.comun$data),beta.reps))
  for (rep in 1:beta.reps) {
    rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(match.phylo.comun$data),taxaShuffle(cophenetic(match.phylo.comun$phy)),abundance.weighted=T,exclude.conspecifics = F))
    print(c(date(),rep))
  }
  
  weighted.bNTI = matrix(c(NA),nrow=ncol(match.phylo.comun$data),ncol=ncol(match.phylo.comun$data))
  for(columns in 1:(ncol(match.phylo.comun$data)-1)) {
    for(rows in (columns+1):ncol(match.phylo.comun$data)) {
      rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
      weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals)
      rm("rand.vals")
    }
  }
  
  rownames(weighted.bNTI) = colnames(match.phylo.comun$data);
  colnames(weighted.bNTI) = colnames(match.phylo.comun$data);
  return(as.dist(weighted.bNTI))
}

raup_crick= function(comun, reps=9){
  require(ecodist) 
  ## count number of sites and total species richness across all plots (gamma)
  n_sites<-nrow(comun)
  gamma<-ncol(comun)
  ##build a site by site matrix for the results, with the names of the sites in the row and col names:
  results<-matrix(data=NA, nrow=n_sites, ncol=n_sites, dimnames=list(row.names(comun), row.names(comun)))
  ##make the comun matrix into a new, pres/abs. matrix:
  ceiling(comun/max(comun))->comun.inc
  ##create an occurrence vector- used to give more weight to widely distributed species in the null model:
  occur<-apply(comun.inc, MARGIN=2, FUN=sum)
  ##create an abundance vector- used to give more weight to abundant species in the second step of the null model:
  abundance<-apply(comun, MARGIN=2, FUN=sum)
  ##make_null:
  ##looping over each pairwise community combination:
  for(null.one in 1:(nrow(comun)-1)){
    for(null.two in (null.one+1):nrow(comun)){
      null_bray_curtis<-NULL
      for(i in 1:reps){
        ##two empty null communities of size gamma:
        com1<-rep(0,gamma)
        com2<-rep(0,gamma)
        ##add observed number of species to com1, weighting by species occurrence frequencies:
        com1[sample(1:gamma, sum(comun.inc[null.one,]), replace=FALSE, prob=occur)]<-1
        com1.samp.sp = sample(which(com1>0),(sum(comun[null.one,])-sum(com1)),replace=TRUE,prob=abundance[which(com1>0)]);
        com1.samp.sp = cbind(com1.samp.sp,1); # head(com1.samp.sp);
        com1.sp.counts = as.data.frame(tapply(com1.samp.sp[,2],com1.samp.sp[,1],FUN=sum)); colnames(com1.sp.counts) = 'counts'; # head(com1.sp.counts);
        com1.sp.counts$sp = as.numeric(rownames(com1.sp.counts)); # head(com1.sp.counts);
        com1[com1.sp.counts$sp] = com1[com1.sp.counts$sp] + com1.sp.counts$counts; # com1;
        #sum(com1) - sum(spXsite[null.one,]); ## this should be zero if everything work properly
        rm('com1.samp.sp','com1.sp.counts');      
        ##same for com2:
        com2[sample(1:gamma, sum(comun.inc[null.two,]), replace=FALSE, prob=occur)]<-1
        com2.samp.sp = sample(which(com2>0),(sum(comun[null.two,])-sum(com2)),replace=TRUE,prob=abundance[which(com2>0)]);
        com2.samp.sp = cbind(com2.samp.sp,1); # head(com2.samp.sp);
        com2.sp.counts = as.data.frame(tapply(com2.samp.sp[,2],com2.samp.sp[,1],FUN=sum)); colnames(com2.sp.counts) = 'counts'; # head(com2.sp.counts);
        com2.sp.counts$sp = as.numeric(rownames(com2.sp.counts)); # head(com2.sp.counts);
        com2[com2.sp.counts$sp] = com2[com2.sp.counts$sp] + com2.sp.counts$counts; # com2;
        # sum(com2) - sum(spXsite[null.two,]); ## this should be zero if everything work properly
        rm('com2.samp.sp','com2.sp.counts');
        null.comun = rbind(com1,com2); # null.comun;
        ##calculate null bray curtis
        null_bray_curtis[i] = distance(null.comun,method='bray-curtis');
      }; # end reps loop
      ## empirically observed bray curtis
      obs.bray = distance(comun[c(null.one,null.two),],method='bray-curtis');
      ##how many null observations is the observed value tied with?
      num_exact_matching_in_null = sum(null_bray_curtis==obs.bray);
      ##how many null values are smaller than the observed *dissimilarity*?
      num_less_than_in_null = sum(null_bray_curtis<obs.bray);
      rc = ((num_less_than_in_null +(num_exact_matching_in_null)/2)/reps)
      ##modification of raup crick standardizes the metric to range from -1 to 1 instead of 0 to 1
      rc = (rc-.5)*2
      results[null.two,null.one] = round(rc,digits=2); ##store the metric in the results matrix
      print(c(null.one,null.two,date()));
    }; ## end null.two loop
  }; ## end null.one loop
  
  results<-as.dist(results)
  return(results)
}
#########################
library(picante)
phylo<-read.tree("otu_tree.tre")
comun<-read.csv("wh.csv",row.names = 1)
dim(comun)
comun[1:5,1:5]
match.phylo.comun = match.phylo.data(phylo, comun)
str(match.phylo.comun)
beta.mntd.weighted = as.matrix(comdistnt(t(match.phylo.comun$data),cophenetic(match.phylo.comun$phy),abundance.weighted=T))
dim(beta.mntd.weighted)
write.csv(beta.mntd.weighted,'betaMNTD_ar.csv', quote=F)
Bnti<-Beta_NTI(phylo,comun)
Bnti=as.matrix(Bnti)
write.csv(Bnti,'Bnti.csv',quote=F)




library(vegan)
library(ecodist)
comun<-read.csv("wh.csv",row.names = 1)
#RC_bray
comun=t(comun)
comun<-as.data.frame(comun)
ab<- raup_crick(comun) 
ab<- as.matrix(ab)
rownames(ab)<-colnames(ab)
colnames(ab)<-colnames(ab)
write.csv(ab, 'rc_bray.csv')

nti.ar<-read.csv("Bnti.csv",row.names=1)
rc.ar<-read.csv("rc_bray.csv",row.names=1)
nti.ar<-as.matrix(as.dist(nti.ar))
rc.ar<-as.matrix(as.dist(rc.ar))

pro1<-function(dat1=matrix(rpois(80,1),nrow=9),dat2=matrix(rpois(80,1),nrow=9))
{df<-as.vector(as.dist(dat1))
df1<-as.vector(as.dist(dat2))
VS<-length(df[df>1.9999999])/length(df)
HS<-length(df[df<(-1.9999999)])/length(df)
df1[abs(df)>1.9999999]<-0
DL<-length(df1[df1>0.9499999])/length(df)
HD<-length(df1[df1<(-0.9499999)])/length(df)
UN<-1-VS-HS-DL-HD
a<-c(VS,HS,DL,HD,UN)
names(a)<-c(" Variable selection"," Homogeneous selection","Diserpal limitation","Homogenizing dispersal","Undominated")
return(a)}

group <- read.delim('group.txt', stringsAsFactors = FALSE)
group<-group[1:2]


dis<-read.csv("Bnti.csv",row.names=1)
env1 <- subset(group, Treat == 'SI')$Names
dis_env1 <- dis[env1,env1]

env2 <- subset(group, Treat == 'SM')$Names
dis_env2 <- dis[env2,env2]

env3 <- subset(group, Treat == 'NI')$Names
dis_env3 <- dis[env3,env3]

env4<- subset(group, Treat == 'NM')$Names
dis_env4 <- dis[env4,env4]


write.csv(dis_env1, 'SI.csv')
write.csv(dis_env2, 'SM.csv')
write.csv(dis_env3, 'NI.csv')
write.csv(dis_env4, 'NM.csv')

dis<-read.csv("rc_bray.csv",row.names=1)
env1 <- subset(group, Treat == 'SI')$Names
dis_env1 <- dis[env1,env1]

env2 <- subset(group, Treat == 'SM')$Names
dis_env2 <- dis[env2,env2]

env3 <- subset(group, Treat == 'NI')$Names
dis_env3 <- dis[env3,env3]

env4<- subset(group, Treat == 'NM')$Names
dis_env4 <- dis[env4,env4]


write.csv(dis_env1, 'SI1.csv')
write.csv(dis_env2, 'SM1.csv')
write.csv(dis_env3, 'NI1.csv')
write.csv(dis_env4, 'NM1.csv')

dis_env1 <- as.vector(as.dist(dis_env1))
dis_env2 <- as.vector(as.dist(dis_env2))
dis_env3 <- as.vector(as.dist(dis_env3))
dis_env4 <- as.vector(as.dist(dis_env4))


dis_env1 <- as.matrix(dis_env1)
write.csv(dis_env1, '1.csv')

dis_env2 <- as.matrix(dis_env2)
write.csv(dis_env2, '2.csv')

dis_env3 <- as.matrix(dis_env3)
write.csv(dis_env3, '3.csv')

dis_env4 <- as.matrix(dis_env4)
write.csv(dis_env4, '4.csv')





b<-read.csv("SI.csv",row.names=1)
b1<-read.csv("SI1.csv",row.names=1)
b<-as.matrix(as.dist(b))
b1<-as.matrix(as.dist(b1))

c<-read.csv("SM.csv",row.names=1)
c1<-read.csv("SM.csv",row.names=1)
c<-as.matrix(as.dist(c))
c1<-as.matrix(as.dist(c1))

d<-read.csv("NI.csv",row.names=1)
d1<-read.csv("NI1.csv",row.names=1)
d<-as.matrix(as.dist(d))
d1<-as.matrix(as.dist(d1))

e<-read.csv("NM.csv",row.names=1)
e1<-read.csv("NM1.csv",row.names=1)
e<-as.matrix(as.dist(e))
e1<-as.matrix(as.dist(e1))

A<-read.csv("Bnti.csv",row.names=1)
A1<-read.csv("rc_bray.csv",row.names=1)
A<-as.matrix(as.dist(A))
A1<-as.matrix(as.dist(A1))


dat1<-b[rownames(b),rownames(b)]
dat2<-b1[rownames(b1),rownames(b1)]
null_SI<-pro1(dat1,dat2)
dat1<-c[rownames(c),rownames(c)]
dat2<-c1[rownames(c1),rownames(c1)]
null_SM<-pro1(dat1,dat2)
dat1<-d[rownames(d),rownames(d)]
dat2<-d1[rownames(d1),rownames(d1)]
null_NI<-pro1(dat1,dat2)
dat1<-e[rownames(e),rownames(e)]
dat2<-e1[rownames(e1),rownames(e1)]
null_NM<-pro1(dat1,dat2)
dat1<-A[rownames(A),rownames(A)]
dat2<-A1[rownames(A1),rownames(A1)]
null_Whole<-pro1(dat1,dat2)


library(RColorBrewer)
library(gplots)
library(ggplot2)
library(reshape2)


SES_ar<-rbind(null_Whole,null_SI,null_SM,null_NI,null_NM)
rownames(SES_ar)<-c("Whole","SI","SM","NI","NM")
write.csv(SES_ar,"ratio.csv",row.names = T)
col<-c("#FF5900", "#00B76D","#0D0DBB","#1f78b4","#e31a1c")
pdf("seed3.pdf")
p<-barplot2(t(SES_ar*100),col=col,las=1,border=NA,space=0.4,
            beside=F,xlab="Group",ylab="Community assembly processes (%)")
dev.off()


df<-melt(SES_ar)
df$Var1<-factor(df$Var1,levels=c("Whole","SI","SM","NI","NM"))
pdf("seed4.pdf",width=8,height=6)
p<-ggplot(df,aes(Var2,value,fill=Var2))
p+geom_bar(stat="identity",position=position_dodge())+
  facet_wrap(~Var1,nrow=2)+
  scale_fill_manual(values=col)+
  theme_bw()
dev.off()
######################
bb <- read.csv("nti.csv",header = TRUE,sep = ",")
library(reshape2)

bb1 <- melt(bb)
##????????
bb$Group <- factor(bb$Group,levels = c("SI","SM","NI","NM"))
bb1$Group <- factor(bb1$Group,levels = c("SI","SM","NI","NM"))

library(multcomp)
bb.sample <- colnames(bb)[2:3]
##?????????????Բ???
test.b <- c()
for (i in bb.sample) {
  fit1 <- aov(as.formula(sprintf("%s ~ Group",i)),
              data = bb)
  tuk1<-glht(fit1,linfct=mcp(Group="Tukey"))
  res1 <- cld(tuk1,alpah=0.05)
  test.b <- cbind(test.b,res1$mcletters$Letters)
}

colnames(test.b) <- colnames(bb)[2:3]
test.b <- melt(test.b)
colnames(test.b) <- c("Group","variable","value")


library(tidyverse)
test.b1 <- bb %>% gather(variable,value,-Group) %>% group_by(variable,Group) %>% 
  summarise(Max = max(value))
test.b11 <- dcast(test.b1,Group~variable)
for (i in 2:ncol(test.b11)) {
  test.b11[,i] <- test.b11[,i] + max(test.b11[,i])*0.1
}
test.b11 <- melt(test.b11)
test.b1 <- merge(test.b1,test.b11,by = c("variable","Group"))
test.b2 <- merge(test.b,test.b1,by = c("variable","Group"))



library(ggplot2)
library(ggpubr)
cbbPalette <- c("#e31a1c","#1f78b4","#FF5900", "#00B76D","#FF5900","#FF5900","#FF5900","#FF5900","#FF5900","#FF5900","#FF5900","#FF5900","#FF5900","#FF5900","#FF5900","#FF5900","#FF5900")
p<-ggplot(bb1,aes(Group,value)) + 
  geom_boxplot(aes(fill = Group)) +
  #geom_point()+
  geom_text(data = test.b2,aes(x = Group,y = value.y,label = value.x),
            size = 5,color = "black",fontface = "bold") +
  scale_fill_manual(values = cbbPalette) +
  ylab("") +
  facet_wrap(.~variable,ncol = 2,scales = "free_y") +
  theme_bw()+
  geom_hline(yintercept = c(-2,2),linetype="dotted")+
  theme(axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_blank(),
        axis.title.y=element_text(colour='black', size=18,face = "bold",vjust = 1.5),
        axis.text.y=element_text(colour='black',size=10),
        axis.text.x=element_text(colour = "black",size = 15,
                                 angle = 0,hjust = 1,vjust = 0.5),
        strip.text = element_text(colour = "black",size = 10,face = "bold"),
        legend.position = "none")
p

ggsave('??NTI1.pdf', p, width = 8, height = 6)

