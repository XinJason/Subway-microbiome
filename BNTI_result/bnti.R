rm(list = ls())
setwd("E:/博士/合作/周欣/bnti")
getwd()
library(ggplot2)
library(reshape2)
library(vegan)
library(ape)
library(picante)
library(NST)
tree<-read.tree("otu_tree.tre")

comun<-read.csv("wh.csv",row.names = 1)
otu <- data.frame(t(comun))

group <- read.csv('Gro.csv', row.names = 1)

set.seed(123)
pnst <- pNST(comm = otu, tree = tree, group = group, phylo.shuffle = TRUE, taxo.null.model = NULL, 
             pd.wd = tempdir(), abundance.weighted = TRUE, rand = 100, nworker = 2, SES = TRUE, RC = T)

#本篇主要是计算 betaMNTD 和 betaNTI，暂且不管其它的（比方说 pNST），因此在结果项里面，主要看这个就行了
#names(pnst)
betaMNTD <- pnst$index.pair
head(betaMNTD)

#输出两两样本的 betaMNTD 和 betaNTI
write.csv(betaMNTD, 'betaMNTD.csv', quote = FALSE, row.names= FALSE)

group<-read.csv("Group.csv",row.names = 1)

Winter <- rownames(subset(group, Seasonal=='Winter'))
bnti_Winter <- subset(betaMNTD, name1 %in% Winter & name2 %in% Winter)
bnti_Winter$Seasonal<-("Winter")

Spring <- rownames(subset(group, Seasonal=='Spring'))
bnti_Spring <- subset(betaMNTD, name1 %in% Spring & name2 %in% Spring)
bnti_Spring$Seasonal<-("Spring")

Summer <- rownames(subset(group, Seasonal=='Summer'))
bnti_Summer <- subset(betaMNTD, name1 %in% Summer & name2 %in% Summer)
bnti_Summer$Seasonal<-("Summer")

Autumn <- rownames(subset(group, Seasonal=='Autumn'))
bnti_Autumn <- subset(betaMNTD, name1 %in% Autumn & name2 %in% Autumn)
bnti_Autumn$Seasonal<-("Autumn")

bnti<-rbind(bnti_Spring,bnti_Summer,bnti_Autumn,bnti_Winter)
bnti<-bnti[11:13]
names(bnti)[1:2]<-c("BNTI","RC")

#write.csv(bnti,"BNTI_RC.csv")


b <- within(bnti, A <- ifelse(BNTI > 2, "VS",   
                           ifelse(BNTI < -2, "HS",   
                                  ifelse(abs(BNTI) < 2 & RC > 0.95, "DL",   
                                         ifelse(abs(BNTI) < 2 & RC < -0.95, "HD",   
                                                ifelse(abs(BNTI) < 2 & abs(RC) < 0.95, "Un", NA))))))  

write.csv(b,"BNTI_RC_pro.csv")

sp<-b[which(b$Seasonal=="Spring"),]
VS<-length(sp[which(sp$A=="VS"),]$A)/length(sp$A)
HS<-length(sp[which(sp$A=="HS"),]$A)/length(sp$A)
HD<-length(sp[which(sp$A=="HD"),]$A)/length(sp$A)
DL<-length(sp[which(sp$A=="DL"),]$A)/length(sp$A)
UN<-1-VS-HS-DL-HD
sp<-as.data.frame(cbind(VS,HS,HD,DL,UN))
sp$Seasonal<-c("Spring")
w<-sp

sp<-b[which(b$Seasonal=="Summer"),]
VS<-length(sp[which(sp$A=="VS"),]$A)/length(sp$A)
HS<-length(sp[which(sp$A=="HS"),]$A)/length(sp$A)
HD<-length(sp[which(sp$A=="HD"),]$A)/length(sp$A)
DL<-length(sp[which(sp$A=="DL"),]$A)/length(sp$A)
UN<-1-VS-HS-DL-HD
sp<-as.data.frame(cbind(VS,HS,HD,DL,UN))
sp$Seasonal<-c("Summer")
su<-sp


sp<-b[which(b$Seasonal=="Autumn"),]
VS<-length(sp[which(sp$A=="VS"),]$A)/length(sp$A)
HS<-length(sp[which(sp$A=="HS"),]$A)/length(sp$A)
HD<-length(sp[which(sp$A=="HD"),]$A)/length(sp$A)
DL<-length(sp[which(sp$A=="DL"),]$A)/length(sp$A)
UN<-1-VS-HS-DL-HD
sp<-as.data.frame(cbind(VS,HS,HD,DL,UN))
sp$Seasonal<-c("Autumn")
au<-sp

sp<-b[which(b$Seasonal=="Winter"),]
VS<-length(sp[which(sp$A=="VS"),]$A)/length(sp$A)
HS<-length(sp[which(sp$A=="HS"),]$A)/length(sp$A)
HD<-length(sp[which(sp$A=="HD"),]$A)/length(sp$A)
DL<-length(sp[which(sp$A=="DL"),]$A)/length(sp$A)
UN<-1-VS-HS-DL-HD
sp<-as.data.frame(cbind(VS,HS,HD,DL,UN))
sp$Seasonal<-c("Winter")
wi<-sp

w<-rbind(w,su,au,wi)

b<-melt(w)

names(b)[2:3]<-c("Pattern","value")
library(ggplot2)
library(ggthemes)

b$Seasonal= factor(b$Seasonal,levels=c("Spring","Summer","Autumn","Winter"),order=TRUE)
#设置颜色
mycolors<-c("#f49128","#194a55","#187c65","#f26115","#c29f62","#83ba9e")
p1<-ggplot(data=b,aes(x=Seasonal,y=value,fill=Pattern)) + 
  geom_bar(stat="identity",position="fill") + 
  scale_fill_manual(values=mycolors)+
  labs(x="Seasonal",y="Percentage (%)",
       fill=" ",title="")+
  theme_bw()+
  theme(axis.title.y=element_text(size=14))+
  theme(legend.text=element_text(size=10))+
  theme(axis.text.x = element_text(size = 12, color = "black"))+
  theme(axis.text.y = element_text(size = 12, color = "black"))+
  theme(axis.ticks.length=unit(0.3,"cm"))+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1,size=11))

p1

ggsave("组装过程.pdf",p1,width = 5, height = 6)














bb <- bnti[-2]
library(reshape2)
names(bb)[2]<-c("Group")
bb1 <- melt(bb)
##????????
bb$Group <- factor(bb$Group,levels = c("Spring","Summer","Autumn","Winter"))
bb1$Group <- factor(bb1$Group,levels = c("Spring","Summer","Autumn","Winter"))

library(multcomp)
bb.sample <- colnames(bb)[1]
##?????????????Բ???
test.b <- c()
for (i in bb.sample) {
  fit1 <- aov(as.formula(sprintf("%s ~ Group",i)),
              data = bb)
  tuk1<-glht(fit1,linfct=mcp(Group="Tukey"))
  res1 <- cld(tuk1,alpah=0.05)
  test.b <- cbind(test.b,res1$mcletters$Letters)
}

colnames(test.b) <- colnames(bb)[1]
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


color <- c("#CD950C", "#104E8B", "#8d9fbb", "#48597e")
library(ggplot2)
cbbPalette <- c("#f49128","#194a55","#187c65","#f26115","#c29f62","#83ba9e")
p<-ggplot(bb1,aes(Group,value)) + 
  geom_boxplot(aes(fill = Group)) +
  geom_text(data = test.b2,aes(x = Group,y = value.y,label = value.x),
            size = 5,color = "black",fontface = "bold") +
  #geom_jitter(size=0.01,aphla=0.5)+
  scale_fill_manual(values = cbbPalette) +
  ylab("BetaNTI") +
  facet_wrap(.~variable,ncol = 3,scales = "free_y") +
  theme_bw()+
  geom_hline(yintercept = c(-2,2),linetype="dotted")+
  theme_bw(base_size = 16) +
  scale_color_manual(values = color)+
  theme(axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=20,vjust = 1.5),
        axis.title.y=element_text(colour='black', size=20,vjust = 1.5),#face = "bold"),
        axis.text.y=element_text(colour='black',size=18),#face = "bold"),
        axis.text.x=element_text(colour='black',size=18,angle = 90,hjust = 1),
        strip.text = element_text(colour = "black",size = 15,face = "bold"),
        legend.position = 'null')
p
ggsave("bnti.pdf",p,width = 4, height = 8)


color <- c("#CD950C", "#104E8B", "#8d9fbb", "#48597e")

library(ggbeeswarm)
#install.packages("ggbeeswarm")
p<-ggplot(bb1,aes(Group,value,color=Group)) + 
  geom_boxplot(aes(fill = Group)) +
  geom_jitter(aes(fill = Group),size=0.001)+
  #geom_beeswarm(cex = 0.2,size=0.01, priority = "descending") +
  geom_text(data = test.b2,aes(x = Group,y = value.y,label = value.x),
            size = 5,color = "black",fontface = "bold") +
  #geom_jitter(size=0.01,aphla=0.5)+
  scale_fill_manual(values = cbbPalette) +
  ylab("BetaNTI") +
  facet_wrap(.~variable,ncol = 3,scales = "free_y") +
  theme_bw()+
  geom_hline(yintercept = c(-2,2),linetype="dotted")+
  theme_bw(base_size = 16) +
  scale_color_manual(values = color)+
  theme(axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=20,vjust = 1.5),
        axis.title.y=element_text(colour='black', size=20,vjust = 1.5),#face = "bold"),
        axis.text.y=element_text(colour='black',size=18),#face = "bold"),
        axis.text.x=element_text(colour='black',size=18,angle = 90,hjust = 1),
        strip.text = element_text(colour = "black",size = 15,face = "bold"),
        legend.position = 'null')
p
ggsave("bnti1.pdf",p,width = 4, height = 6)



  