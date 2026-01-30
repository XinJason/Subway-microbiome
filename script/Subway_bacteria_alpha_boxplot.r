#!/usr/bin/env Rscript
rm(list = ls())
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
    make_option(c("-t", "--type"), type="character", default="shannon",
                help="type of alpha index; 指数列名 [default %default]"),
    make_option(c("-d", "--design"), type="character", default="result/metadata.txt",
                help="design file; 实验设计文件 [default %default]"),
    make_option(c("-n", "--group"), type="character", default="Station",
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
##自己加的，调整分组顺序命令"#B5E1E2","#1D4C7C","#F28034","#C3C5C7","Winter", "Spring", "Summer", "Autumn"
index$group <- factor(index$group, levels=c("IHS", "UHS", "SHS", ordered=TRUE))
##自己定义颜色 "#226600","#3e5780","#594359","Intercity", "UrbanHub", "Residential"
#color <- c( "#e3e5d7","#e7dbc1","#e7cfd4", "#c6e3d8", "#e7f3ee","#c8ccc1","#aabdc4","#8d9fbb", "#8d958f", "#6e746a", "#426600", "#48597e")
color <- c("#C23B2F","#FF9600","#015F96")
##加载颜色包，自动配色
#color <- c(brewer.pal(8,"Set1"))
p = ggplot(index, aes(x=group, y=index[[opts$type]], color=group)) +
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
  labs(x="Groups", y=paste(opts$type, "index"),fill="Diversity variation") + 
  geom_text(data=index, aes(x=group, y=y, color=group, label= stat)) +
  geom_jitter(position=position_jitter(0.17), size=1, alpha=0.5)+
  #theme_bw(base_line_size = 1.05,base_rect_size = 1.05)+
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
ggsave(paste(opts$output, "_station_alpha.pdf", sep=""), p, width = opts$width, height = opts$height)
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
#design0<-design[design$Sample == "All",]
design0<-design
index0<-index[rownames(design0),]
design0<-design0[rownames(index0),]
index0<-as.data.frame(index0)
##自己加的，调整分组顺序命令
index0$group <- factor(index0$group, levels=c( "Jan", "Feb", "Mar","Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec", ordered=TRUE))
##自己定义颜色
color <- c("#e3e5d7","#e7dbc1","#e7cfd4", "#e7f3ee","#c8ccc1","#c6e3d8", "#aabdc4","#8d9fbb", "#8d958f", "#6e746a", "#426600", "#48597e")

#index0 <- na.omit(index0)
p0=ggplot(index0, aes(x=group, y=index0[[opts$type]], color=group),na.rm = TRUE) +
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
  labs(x="group", y=paste(opts$type, "index"),fill="Temperature zone") + 
  geom_text(data=index, aes(x=group, y=y, color=group, label= stat)) +
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.5)+
 # theme_bw(base_line_size = 1.05,base_rect_size = 1.05)+
  theme(axis.text=element_text(colour='black',size=9), legend.position = "none")+
  # 不需要填充的时候去掉fill及对fill的补充参数geom,alpha等即可
  # scale_fill_manual(values=colorRampPalette(c("gray98","#6191BF","#3D3E83"))(4))+
  scale_color_manual(values = color)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())  
# 颜色可以自己设置，或者直接用scale_color_brewer()
#p0=p0+ scale_color_npg()  
# 读取usearch alpha文件
alpha = read.table("result/alpha/TimeLine.txt", header=T, row.names=1, sep="\t", comment.char="") 
# 读取实验设计
design = read.table("result/metadata.txt", header=T, row.names=1, sep="\t", comment.char="") 

data("alpha")
##自己加的，调整分组顺序命令"#226600","#3e5780","#594359",
alpha$group <- factor(alpha$group, levels=c("Jan", "Feb", "Mar","Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov",  "Dec", ordered=TRUE))
color <- c("#e3e5d7","#e7dbc1","#e7cfd4","#b5aabd", "#dB9c7c","#e7f3ee","#c8ccc1","#c6e3d8", "#aabdc4","#8d9fbb", "#8d958f", "#6e746a", "#426600","#544559", "#48597e")
#color <- c("#C23B2F","#FF9600","#015F96")
##自己加的，调整分组顺序命令"IHS", "UHS", "SHS",
alpha$Station <- factor(alpha$Station, levels=c( ordered=TRUE))
#alpha$group <- factor(alpha$group, levels=c("day1", "day4", "day7", "day14","day21","day28","day35", "day42",ordered=TRUE))
p1 <- ggplot(
  subset(alpha),#(alpha, shannon > 2.2 & richness > 600 )
  aes(x = group, y = shannon, colour = Station)
) +
  geom_point(alpha = 1) +
  #size控制线的粗细
  geom_smooth(aes(group=Station),span=0.5,size=1.5, se=FALSE,method='lm',formula=y~I(poly(x,2)))+ theme_classic()+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.5)+
 # theme_bw(base_line_size = 1.05,base_rect_size = 1.05)+
  theme(axis.text=element_text(colour='black',size=9))+
  # 不需要填充的时候去掉fill及对fill的补充参数geom,alpha等即可
  scale_color_manual(values = color)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())  
p2 <-p1+ scale_color_npg(palette = c("nrc"))

# 保存pdf和png格式方便查看和编辑
ggsave("result/alpha/Stations_TimeLine_shannon.pdf", p1, width = 9, height = 4)



#################################################################################


# 1.2 解析命令行
# 设置清华源加速下载
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 判断命令行解析是否安装，安装并加载
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages("optparse", repos=site)
  require("optparse",character.only=T) 
}

# 解析参数-h显示帮助信息invsimpson, richness, simpson, shannon, ACE,chao1
if (TRUE){
  option_list <- list(
    make_option(c("-i", "--input"), type="character", default="result/alpha/TimeLine.txt",
                help="Input table file to read; Alpha指数文件 [default %default]"),
    make_option(c("-t", "--type"), type="character", default="richness",
                help="type of alpha index; 指数列名 [default %default]"),
    make_option(c("-d", "--design"), type="character", default="result/metadata.txt",
                help="design file; 实验设计文件 [default %default]"),
    make_option(c("-n", "--group"), type="character", default="group",
                help="name of group type; 分组列名 [default %default]"),
    make_option(c("-w", "--width"), type="numeric", default=10,
                help="Width of figure; 图片宽 [default %default]"),
    make_option(c("-e", "--height"), type="numeric", default=8,
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
package_list <- c("reshape2","ggplot2","devtools","bindrcpp",
                  "ggthemes","agricolae","dplyr","gridExtra","ggsci")
# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}
# 2.2 安装bioconductor常用包
package_list <- c("digest")
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
data("alpha")
##自己加的，调整分组顺序命令
color <- c("#C23B2F","#FF9600","#015F96")
alpha$group <- factor(alpha$group, levels=c("Jan", "Feb", "Mar","Apr", "May", "Jun", "Jul", "Aug", "Sept","Oct", "Nov", "Dec",ordered=TRUE))
p1 <- ggplot(
  subset(alpha),#(alpha, shannon > 2.2 & richness > 600 )
  aes(x = group, y = richness, colour = Station_Type)
) +
  geom_point(alpha = 1) +
  scale_color_manual(values = color)+
  #size控制线的粗细
  geom_smooth(aes(group=Station_Type),span=0.5,size=1.5, se=FALSE,method='lm',formula=y~I(poly(x,2)))+ theme_classic()+
  geom_jitter(height = 0.05) 
#p2=p1 + scale_color_npg()
p1
# 保存pdf和png格式方便查看和编辑
ggsave("result/alpha/Station_Type_TimeLine_richness.pdf", p1, width = 8, height = 5)

##################################################3
# 读取usearch alpha文件
alpha = read.table(opts$input, header=T, row.names=1, sep="\t", comment.char="") 
# 读取实验设计
design = read.table(opts$design, header=T, row.names=1, sep="\t", comment.char="") 
data("alpha")
##自己加的，调整分组顺序命令
##自己加的，调整分组顺序命令"#226600","#3e5780","#594359",
alpha$group <- factor(alpha$group, levels=c("Jan", "Feb", "Mar","Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov",  "Dec", ordered=TRUE))
color <- c("#e3e5d7","#e7dbc1","#e7cfd4","#b5aabd", "#dB9c7c","#e7f3ee","#c8ccc1","#c6e3d8", "#aabdc4","#8d9fbb", "#8d958f", "#6e746a", "#426600","#544559", "#48597e")
p1 <- ggplot(
  subset(alpha),#(alpha, shannon > 2.2 & richness > 600 )
  aes(x = group, y = shannon, colour = Station)
) +
  geom_point(alpha = 1) +
  scale_color_manual(values = color)+
  #size控制线的粗细
  geom_smooth(aes(group=Station),span=0.5,size=1.5, se=FALSE,method='lm',formula=y~I(poly(x,2)))+ theme_classic()+
  geom_jitter(height = 0.05) 
#p2=p1 + scale_color_npg()
p1
# 保存pdf和png格式方便查看和编辑
ggsave("result/alpha/Station_TimeLine_shannon.pdf", p1, width = 8, height = 5)

