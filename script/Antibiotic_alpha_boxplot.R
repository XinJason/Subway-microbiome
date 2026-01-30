#!/usr/bin/env Rscript
rm(list = ls())
# # 完整参数，输出文件名默认为alpha指数类型
# Rscript ./script/alpha_boxplot.r -i alpha/alpha.txt -t richness \
# -d doc/design.txt -n group \
# -o alpha/richness \
# -w 4 -e 2.5
# Rscript ./script/alpha_boxplot.r -t chao1 
options(warn = -1) # Turn off warning
# 1.2 解析命令行
# 设置清华源加速下载
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 判断命令行解析是否安装，安装并加载
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages("optparse", repos=site)
  require("optparse",character.only=T) 
}
# 解析参数-h显示帮助信息chao1,dominance, richness, simpson, shannon, ACE
if (TRUE){
  option_list <- list(
    make_option(c("-i", "--input"), type="character", default="result/alpha/vegan.txt",
                help="Input table file to read; Alpha指数文件 [default %default]"),
    make_option(c("-t", "--type"), type="character", default="richness",
                help="type of alpha index; 指数列名 [default %default]"),
    make_option(c("-d", "--design"), type="character", default="result/metadata.txt",
                help="design file; 实验设计文件 [default %default]"),
    make_option(c("-n", "--group"), type="character", default="Season",
                help="name of group type; 分组列名 [default %default]"),
    make_option(c("-w", "--width"), type="numeric", default=5,
                help="Width of figure; 图片宽 [default %default]"),
    make_option(c("-e", "--height"), type="numeric", default=4,
                help="Height of figure; 图片高 [default %default]"),
    make_option(c("-o", "--output"), type="character", default="",
                help="output directory or prefix; 输出文件前缀, 有统计表txt和矢量图pdf [default %default]")
  )
  opts <- parse_args(OptionParser(option_list=option_list))
  
  # 如果没有设置输出，根据默认文件夹和类型参数设置输出位置
  if (opts$output==""){
    opts$output=paste("result/alpha/",opts$type, sep = "")}
  
  # 显示输入输出确认是否正确
  print(paste("The input file is ", opts$input,  sep = ""))
  print(paste("Type of alpha index is ", opts$type,  sep = ""))
  print(paste("The design file is ", opts$design,  sep = ""))
  print(paste("The group name is ", opts$group,  sep = ""))
  print(paste("Output figure width ", opts$width,  sep = ""))
  print(paste("Output figure height ", opts$height,  sep = ""))
  print(paste("The output file prefix is ", opts$output, sep = ""))
}
# 2.1 安装CRAN来源常用包
# 依赖包列表：参数解析、数据变换、绘图和开发包安装、安装依赖、ggplot主题
package_list <- c("reshape2","ggplot2","devtools","bindrcpp","RColorBrewer","ggrepel","ggsci",
                  "ggthemes","agricolae","dplyr")
# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

# 2.2 安装bioconductor常用包
package_list <- c("digest","ggsci")
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    source("https://bioconductor.org/biocLite.R")
    biocLite(p)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}
# 2.3 安装Github常用包
# 参数解析、数据变换、绘图和开发包安装
package_list <- c("kassambara/ggpubr")
for(p in package_list){
  q=unlist(strsplit(p,split = "/"))[2]
  if(!suppressWarnings(suppressMessages(require(q, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install_github(p)
    suppressWarnings(suppressMessages(library(q, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}
# 3. 读取输入文件
# 读取usearch alpha文件
alpha = read.table(opts$input, header=T, row.names=1, sep="\t", comment.char="") 
# 读取实验设计
design = read.table(opts$design, header=T, row.names=1, sep="\t", comment.char="") 

#design0<-design[design$Station_Type == "Residential",]
#design = design0
# 提取样品组信息,默认为group可指定
sampFile = data.frame(group=design[,opts$group],
                      sample=row.names(design), 
                      row.names = row.names(design))
# 数据筛选，筛选两文件中共有
idx = rownames(sampFile) %in% rownames(alpha) # match design with alpha
sampFile = sampFile[idx,]
alpha = alpha[rownames(sampFile),] 
## richness index
# add design to alpha
#index = cbind(alpha[rownames(design),]$richness, sampFile) 
index = cbind(alpha[rownames(sampFile),][[opts$type]], sampFile) 
colnames(index) = c(opts$type,"group","sample") # add richness colname is value
# 4. 统计与绘图
# 统计各组间差异
#model = aov(richness ~ group, data=index)
model = aov(index[[opts$type]] ~ group, data=index)
# 计算Tukey显著性差异检验
Tukey_HSD <- TukeyHSD(model, ordered = TRUE, conf.level = 0.95)
# 提取比较结果
Tukey_HSD_table <- as.data.frame(Tukey_HSD$group) 
Tukey_HSD_table
# LSD检验，添加差异组字母
out <- LSD.test(model,"group", p.adj="none")
stat = out$groups
# 分组结果添入Index
index$stat=stat[as.character(index$group),]$groups
# 设置分组位置为各组y最大值+高的3%
max=max(index[,c(opts$type)])
min=min(index[,opts$type])
x = index[,c("group",opts$type)]

# 下方的richness如何替换为变量
# y = x%>% group_by(group)%>% summarise(Max=max(richness))
y = x %>% group_by(group) %>% summarise_(Max=paste('max(',opts$type,')',sep=""))

y=as.data.frame(y)
rownames(y)=y$group
index$y=y[as.character(index$group),]$Max + (max-min)*0.05
##自己加的，调整分组顺序命令
index$group <- factor(index$group, levels=c("Spring", "Summer", "Autumn","Winter", ordered=TRUE))
##自己定义颜色 "#226600","#3e5780","#594359","Intercity", "UrbanHub", "Residential"
#color <- c( "#e3e5d7","#e7dbc1","#e7cfd4", "#c6e3d8", "#e7f3ee","#c8ccc1","#aabdc4","#8d9fbb", "#8d958f", "#6e746a", "#426600", "#48597e")
color <- c("#C3C5C7","#B5E1E2","#1D4C7C","#F28034")
##加载颜色包，自动配色
#color <- c(brewer.pal(8,"Set1"))
p = ggplot(index, aes(x=group, y=index[[opts$type]], color=group)) +
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
  labs(x="Groups", y=paste(opts$type, "index"),fill="Diversity variation") + 
  geom_text(data=index, aes(x=group, y=y, color=group, label= stat)) +
  geom_jitter(position=position_jitter(0.17), size=1, alpha=0.5)+
  theme_bw(base_line_size = 1.05,base_rect_size = 1.05)+
  theme(axis.text=element_text(colour='black',size=9))+
# 不需要填充的时候去掉fill及对fill的补充参数geom,alpha等即可
# 主题调整
# 颜色可以自己设置，或者直接用scale_color_brewer()

 # scale_fill_manual(values=colorRampPalette(c("gray98","#6191BF","#3D3E83"))(4))+
  scale_color_manual(values = color)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())  
#p=p+ scale_color_npg()  
# 颜色可以自己设置，或者直接用scale_color_brewer()
# 5. 保存图表
# 保存pdf和png格式方便查看和编辑
ggsave(paste(opts$output, "_Date_alpha.pdf", sep=""), p, width = opts$width, height = opts$height)
# ggsave(paste(opts$output, ".png", sep=""), p, width = opts$width, height = opts$height)
# 保存一个制表符，解决存在行名时，列名无法对齐的问题
write.table("\t", file=paste(opts$output,".txt",sep=""),append = F, quote = F, eol = "", row.names = F, col.names = F)
# 保存统计结果，有waring正常
write.table(Tukey_HSD_table, file=paste(opts$output,".txt",sep=""), append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T)
# 提示工作完成
print(paste("Output in ", opts$output, ".txt/pdf finished.", sep = ""))

#########################################################
# 提取样品组信息,默认为group可指定
# 数据筛选，筛选两文件中共有
idx = rownames(design) %in% rownames(index) # match design with alpha
design1 = design[idx,]
index1 = index[idx,] 
## richness index
#合并pH数据，和随机构建数据框
df<-cbind(index1,design1)

df1<-df[df$Type == "Rhizo",]

df<-df1

design0<-design[design$Group == "GT44",]
index0<-index[rownames(design0),]
design0<-design0[rownames(df),]
index0<-as.data.frame(index0)
#df1<-df[df$group2 == "Year one" | env.data$Location == "middle",]
##自己加的，调整分组顺序命令
index0$group <- factor(index0$group, levels=c("Winter", "Spring", "Summer", "Autumn", ordered=TRUE))
##自己定义颜色
color <- c("#B5E1E2","#1D4C7C","#F28034","#C3C5C7")

#index0 <- na.omit(index0)
p0=ggplot(index0, aes(x=group, y=index0[[opts$type]], color=group),na.rm = TRUE) +
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
  labs(x="GT44", y=paste(opts$type, "index"),fill="Temperature zone") + 
  geom_text(data=index, aes(x=group, y=y, color=group, label= stat)) +
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.5)+
  theme_bw(base_line_size = 1.05,base_rect_size = 1.05)+
  theme(axis.text=element_text(colour='black',size=9), legend.position = "none")+
  # 不需要填充的时候去掉fill及对fill的补充参数geom,alpha等即可
  # scale_fill_manual(values=colorRampPalette(c("gray98","#6191BF","#3D3E83"))(4))+
  scale_color_manual(values = color)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())  
# 颜色可以自己设置，或者直接用scale_color_brewer()

ggsave(paste(opts$output, "_GT44_alpha.pdf", sep=""), p0, width = opts$width, height = opts$height)


df<-cbind(index1,design1)
design0<-design[design$Group == "GT42",]
index0<-index[rownames(design0),]
design0<-design0[rownames(index0),]
index0<-as.data.frame(index0)

#df1<-df[df$group2 == "Year one" | env.data$Location == "middle",]
##自己加的，调整分组顺序命令
index0$group <- factor(index0$group, levels=c("Winter", "Spring", "Summer", "Autumn", ordered=TRUE))
##自己定义颜色
color <- c("#B5E1E2","#1D4C7C","#F28034","#C3C5C7")

#index0 <- na.omit(index0)
p1=ggplot(index0, aes(x=group, y=index0[[opts$type]], color=group),na.rm = TRUE) +
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
  labs(x="GT42", y=paste(opts$type, "index"),fill="Temperature zone") + 
  geom_text(data=index, aes(x=group, y=y, color=group, label= stat)) +
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.5)+
  theme_bw(base_line_size = 1.05,base_rect_size = 1.05)+
  theme(axis.text=element_text(colour='black',size=9), legend.position = "none")+
  # 不需要填充的时候去掉fill及对fill的补充参数geom,alpha等即可
  # scale_fill_manual(values=colorRampPalette(c("gray98","#6191BF","#3D3E83"))(4))+
  scale_color_manual(values = color)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())  


ggsave(paste(opts$output, "_GT42_alpha.pdf", sep=""), p1, width = opts$width, height = opts$height)
######################################
df<-cbind(index1,design1)
design0<-design[design$Group == "GL05136",]
index0<-index[rownames(design0),]
design0<-design0[rownames(index0),]
index0<-as.data.frame(index0)

#df1<-df[df$group2 == "Year one" | env.data$Location == "middle",]
##自己加的，调整分组顺序命令
index0$group <- factor(index0$group, levels=c("Jun", "Jul", "Aug", "Sep", "Oct", "Nov", ordered=TRUE))
##自己定义颜色
color <- c( "#B5E1E2","#426600","#8d9fbb", "#6e746a",  "#48597e", "#1D4C7C")

#index0 <- na.omit(index0)
p2=ggplot(index0, aes(x=group, y=index0[[opts$type]], color=group),na.rm = TRUE) +
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
  labs(x="GL05136", y=paste(opts$type, "index"),fill="Temperature zone") + 
  geom_text(data=index, aes(x=group, y=y, color=group, label= stat)) +
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.5)+
  theme_bw(base_line_size = 1.05,base_rect_size = 1.05)+
  theme(axis.text=element_text(colour='black',size=9), legend.position = "none")+
  # 不需要填充的时候去掉fill及对fill的补充参数geom,alpha等即可
  # scale_fill_manual(values=colorRampPalette(c("gray98","#6191BF","#3D3E83"))(4))+
  scale_color_manual(values = color)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())  

ggsave(paste(opts$output, "_GL05136_alpha.pdf", sep=""), p2, width = opts$width, height = opts$height)

#图片组合
library(cowplot)
g=cowplot::plot_grid(p, p0,p1,p2,ncol= 2, rel_widths = c(1, 1), labels=LETTERS[1:4])
ggsave(paste(opts$output, "Combined_alpha.pdf", sep=""), g, width = 8, height = 5)

########################################################
# 读取usearch alpha文件
alpha = read.table("result/alpha/TimeLine.txt", header=T, row.names=1, sep="\t", comment.char="") 

# 读取实验设计
design = read.table("result/metadata.txt", header=T, row.names=1, sep="\t", comment.char="") 

data("alpha")
##自己加的，调整分组顺序命令"#226600","#3e5780","#594359",
alpha$Date <- factor(alpha$Date, levels=c("Jun", "Jul", "Aug", "Sep", "Oct", "Nov",  ordered=TRUE))
color <- c( "#B5E1E2","#426600","#8d9fbb", "#6e746a",  "#48597e", "#1D4C7C")
##自己加的，调整分组顺序命令"#d1c7ae", "#8F7C00", "#104E8B", "#7E89A4","Winter", "Spring", "Summer", "Autumn"
alpha$Group <- factor(alpha$Group, levels=c("GT42", "GT44", "GL05136", ordered=TRUE))

#"#e3e5d7","#e7dbc1","#e7cfd4","#b5aabd", "#dB9c7c","#e7f3ee","#c8ccc1","#c6e3d8", "#aabdc4","#8d9fbb", "#8d958f", "#6e746a", "#426600","#544559", "#48597e"
#alpha$group <- factor(alpha$group, levels=c("day1", "day4", "day7", "day14","day21","day28","day35", "day42",ordered=TRUE))
p1 <- ggplot(
  subset(alpha),#(alpha, shannon > 2.2 & richness > 600 )
  aes(x = Date, y = richness, colour = Group)
) +
  geom_point(alpha = 1) +
  #size控制线的粗细
  geom_smooth(aes(group=Group),span=0.5,size=1.5, se=FALSE,method='lm',formula=y~I(poly(x,2)))+ theme_classic()+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.5)+
  theme_bw(base_line_size = 1.05,base_rect_size = 1.05)+
  theme(axis.text=element_text(colour='black',size=9))+
  # 不需要填充的时候去掉fill及对fill的补充参数geom,alpha等即可
  scale_color_manual(values = color)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())  
p2 <-p1+ scale_color_npg(palette = c("nrc"))

# 保存pdf和png格式方便查看和编辑
ggsave("result/alpha/TimeLine_richness.pdf", p1, width = 9, height = 5)


###############################################
