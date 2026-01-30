#!/usr/bin/env Rscript
options(warn = -1) # Turn off warning
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
    make_option(c("-t", "--type"), type="character", default="shannon",
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
#地铁细菌颜色 colour <- c("#8d9fbb","#CD950C","#104E8B","#48597e")
colour <- c("#38761d", "#16537e", "#4c1130")
my_colors <- c(
     "Intercity_hub_station" = "#38761d",
     "Residential_area_station" = "#16537e",
     "Urban_hub_station" = "#4c1130")
##自己加的，调整分组顺序命令
alpha$group <- factor(alpha$group, levels=c("Dec","Jan", "Feb", "Mar","Apr", "May", "Jun", "Jul", "Aug", "Sept","Oct", "Nov", ordered=TRUE))
p1 <- ggplot(
  subset(alpha),#(alpha, shannon > 2.2 & richness > 600 )
  aes(x = group, y = richness, colour = Station_Type)
) +
  geom_point(alpha = 1) +
#size控制线的粗细
  geom_smooth(aes(group=Station_Type),span=0.5,size=1.5, se=FALSE,method='lm',formula=y~I(poly(x,2)))+ theme_classic()+
  geom_jitter(height = 0.05)+
# 3. 应用自定义颜色（关键步骤）
scale_color_manual(values = my_colors) +  # 线条和点的边框色
scale_fill_manual(values = my_colors) +
  theme_classic() +  # 基础主题（先使用经典主题，再修改边框）
  theme(
    # 1. 取消默认的轴线（避免与边框重叠）
    axis.line = element_blank(),
    # 2. 添加完整的面板边框（四周边框）
    panel.border = element_rect(
      colour = "black",  # 边框颜色（可改为灰色#999999更柔和）
      fill = NA,         # 不填充面板（必须设为NA，否则会覆盖背景）
      size = 1           # 边框粗细（根据需要调整）
    ),
    # 3. 调整刻度线长度（避免刻度超出边框）
    axis.ticks.length = unit(-0.15, "cm"),  # 负号表示刻度线向内
    # 4. 调整坐标轴文本位置（避免被刻度线遮挡）
    axis.text.x = element_text(
      angle = 0, hjust = 1, size = 10,
      margin = margin(t = 5)  # 顶部留白，避免文本与刻度重叠
    ),
    axis.text.y = element_text(
      size = 10,
      margin = margin(r = 5)  # 右侧留白
    ),
    # 其他主题参数（保持与之前一致）
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
    legend.position = "top"
  ) +
  labs(x = "Month", y = "Bacterial Richness", colour = "Station Type", fill = "Station Type")
  # 点的填充色
#p2=p1 + scale_color_npg()
p1
# 保存pdf和png格式方便查看和编辑
ggsave("result/alpha/TimeLine_shannon.pdf", p1, width = 7, height = 5)

##绘制箱线图初步观察群落Beta多样性高低水平
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

# ggsave(paste(opts$output, ".png", sep=""), p, width = opts$width, height = opts$height)

# 保存一个制表符，解决存在行名时，列名无法对齐的问题
write.table("\t", file=paste(opts$output,".txt",sep=""),append = F, quote = F, eol = "", row.names = F, col.names = F)
# 保存统计结果，有waring正常
write.table(Tukey_HSD_table, file=paste(opts$output,".txt",sep=""), append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T)



# 提示工作完成
print(paste("Output in ", opts$output, ".txt/pdf finished.", sep = ""))

