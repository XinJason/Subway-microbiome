#!/usr/bin/env Rscript
options(warn = -1)
# 设置清华源加速下载
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 判断命令行解析是否安装，安装并加载
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages("optparse", repos=site)
  require("optparse",character.only=T) 
}
# 解析命令行
if (TRUE){
  option_list = list(
    make_option(c("-t", "--type"), type="character", default="p",
                help="Taxonomy level p c o f g; 分类学级别, 门phylum 纲class 目order 科family 属genus [default %default]"),
    make_option(c("-i", "--input"), type="character", default="",
                help="Merged taxonomy file; 分类学合并结果 [default %default]"),
    make_option(c("-d", "--design"), type="character", default="result/metadata.txt",
                help="design file; 实验设计文件 [default %default]"),
    make_option(c("-n", "--group"), type="character", default="Month",
                help="name of group type; 分组列名 [default %default]"),
    make_option(c("-b", "--number"), type="numeric", default=12,
                help="Number taxonomy for showing; 展示分类数量 [default %default]"),
    make_option(c("-w", "--width"), type="numeric", default=6,
                help="Width of figure; 图片宽 [default %default]"),
    make_option(c("-e", "--height"), type="numeric", default=4,
                help="Height of figure; 图片高 [default %default]"),
    make_option(c("-o", "--output"), type="character", default="",
                help="output directory or prefix; 输出文件前缀, 通常会有统计表txt、矢量图pdf和位图png [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))
  # 调置如果无调设置输出，根据其它参数设置默认输出
  if (opts$input==""){opts$input=paste("result/tax/sum_",opts$type,".txt", sep = "")}
  if (opts$output==""){opts$output=paste("result/tax/stackplot_",opts$type,"",sep = "")}
  # 显示输入输出确认是否正确
  print(paste("Merged taxonomy file is ", opts$input,  sep = ""))
  print(paste("Taxonomy level is ", opts$type,  sep = ""))
  print(paste("Number taxonomy for showing is ", opts$number,  sep = ""))
  print(paste("The design file is ", opts$design,  sep = ""))
  print(paste("The group name is ", opts$group,  sep = ""))
  print(paste("Output figure width ", opts$width,  sep = ""))
  print(paste("Output figure height ", opts$height,  sep = ""))
  print(paste("The output file prefix is ", opts$output, sep = ""))
}
# 2.1 安装CRAN来源常用包
# 依赖包列表：参数解析、数据变换、绘图和开发包安装、安装依赖、ggplot主题
package_list = c("reshape2","ggplot2","vegan","ggsci","ggalluvial")
# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}
# 2.2 安装bioconductor常用包
package_list = c("digest","ggrepel")
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    source("https://bioconductor.org/biocLite.R")
    biocLite(p)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}
# 2.3 安装Github常用包
package_list = c("kassambara/ggpubr")
for(p in package_list){
  q=unlist(strsplit(p,split = "/"))[2]
  if(!suppressWarnings(suppressMessages(require(q, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install_github(p)
    suppressWarnings(suppressMessages(library(q, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}
# 3. 读取输入文件
# 读取样品分类学文件
tax_sample = read.table(opts$input, header=T, row.names= 1, sep="\t", comment.char="") 

# 读取实验设计
design = read.table(opts$design, header=T, row.names= 1, sep="\t", comment.char="") 
# 提取样品组信息,默认为genotype可指定
#design0<-design[design$Region == "Hulunbeier",]
#tax_sample0<-tax_sample[,rownames(design0)]
#index0<-as.data.frame(index0)
#design<- design0
#tax_sample<- tax_sample0
sampFile = data.frame(group=design[,opts$group],
                      sample=row.names(design), 
                      row.names = row.names(design))

# 数据筛选，筛选两文件中共有
idx = rownames(sampFile) %in% colnames(tax_sample) # match design with alpha
sampFile = sampFile[idx,]
tax_sample = tax_sample[,rownames(sampFile)] 
# 检查，每组是否标准化为100%
# colSums(tax_sample)
# 按丰度降序排序
mean_sort = tax_sample[(order(-rowSums(tax_sample))), ]
mean_sort = as.data.frame(mean_sort)

# 筛选前X类，其它归为other，可设置不同组数
other = colSums(mean_sort[opts$number:dim(mean_sort)[1], ])
mean_sort = mean_sort[1:(opts$number - 1), ]
mean_sort = rbind(mean_sort,other)
rownames(mean_sort)[opts$number] = c("Low abundance")
# 再次检验计算是否出错
# colSums(mean_sort)
# 保存变量备份，并输出至文件
merge_tax=mean_sort
write.table("\t", file=paste(opts$output,"_sample.txt",sep=""),append = F, quote = F, eol = "", row.names = F, col.names = F)
write.table(merge_tax, file=paste(opts$output,"_sample.txt",sep=""), append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T)

# 添加分类学列
mean_sort$tax = rownames(mean_sort)
data_all = as.data.frame(melt(mean_sort, id.vars=c("tax")))
# data_all$tax  = factor(data_all$tax, levels=rownames(mean_sort))   # set taxonomy order
data_all = merge(data_all, sampFile, by.x="variable", by.y = "sample")
#colourCount = length(unique(p$data$Taxonomy))
#getPalette = colorRampPalette(brewer.pal(9, "Set1"))
# 更改颜色，把彩虹色改为CMYK色
#"#544559","#b5aabd", "#aabdc4", "#d1c7ae", "#e5d9da" , 
# "#8d9fbb", "#c1c2d2", "#c6d7e3", "#e0e7f5"
#"#aa3e53", "#dB9c7c", "#cf9198", "#e7cfd4", "#f6f0ee"
#"#1a3b30", "#509579", "#99bdbd", "#c6e3d8", "#e7f3ee"
#"#957064", "#c1aea1", "#c4c7b4", "#ecd6bc", "#ebe7e4"
#"#6e746a", "#8d958f" ,"#c8ccc1", "#e7dbc1", "#e3e5d7" 

colours <- c("#ebe7e4","#d1c7ae", "#e7f3ee","#c6e3d8", "#c8ccc1","#aabdc4","#e7cfd4", "#8d958f",  "#6e746a",  "#8F7C00","#FFCC99", "#48597e","#E69F00", "#cf9198", "#dB9c7c");
#class.otu.2 <- melt(t(class.tab.otu), measure.vars = colnames(class.tab.otu))
#class.seq.2 <- melt(t(class.tab.seq), measure.vars = colnames(class.tab.seq))
##自己加的，调整分组顺序命令"Jan", "Feb", "Mar","Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec"
data_all$group <- factor(data_all$group,levels=c(  "Mar","Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", ordered=TRUE))
#data_all$group <- factor(data_all$group,levels=c( "Intercity", "UrbanHub", "Residential", ordered=TRUE))

p = ggplot(data_all, aes(x=variable, y = value, fill = tax )) + 
  geom_bar(stat = "identity",position="fill", width=1)+ 
  scale_y_continuous(labels = scales::percent) + 
  #这条命令是按自己颜色进行设置
  scale_fill_manual(values=colours[1:14])+
  theme_light(base_line_size = 0)+
  # 分面，进一步按group分组，x轴范围自由否则位置异常，swith设置标签底部，并调置去除图例背景
  facet_grid( ~ group, scales = "free_x", switch = "x") +  theme(strip.background = element_blank())+
  # 关闭x轴刻度和标签
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+ xlab("Groups")+
  ylab("Percentage (%)")+ theme_classic()+theme(axis.text.x=element_text(size=8,angle=45,vjust=1, hjust=1), axis.title.y=element_text(size=20),
                                                               axis.text.y=element_text(size=16),legend.title=element_text(size=20),legend.text=element_text(size=20))
#p= p+scale_color_npg(name="group")
p
# 保存pdf和png格式方便查看和编辑
ggsave(paste(opts$output, "_sample.pdf", sep=""), p, width = opts$width, height = opts$height)
# ggsave(paste(opts$output, "_sample.png", sep=""), p, width = opts$width, height = opts$height)
print(paste(opts$output, "_sample.pdf/txt finished.", sep = ""))
# 转置样品名添加组名，并去除多余的两个样品列
mat_t = t(merge_tax)
mat_t2 = merge(sampFile, mat_t, by="row.names")
mat_t2 = mat_t2[,c(-1,-3)]

# 按组求均值，转置，再添加列名
mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean
mat_mean_final = do.call(rbind, mat_mean)[-1,]
geno = mat_mean$group
colnames(mat_mean_final) = geno
# 保存变量备份，并输出至文件
mean_sort=as.data.frame(mat_mean_final)
write.table(mean_sort, file=paste(opts$output,"_group.txt",sep=""),  append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T)

# 数据转换长表格并绘图
mean_sort$tax = rownames(mean_sort)
data_all = as.data.frame(melt(mean_sort, id.vars=c("tax")))
# 设置分类学顺序，默认字母，可选丰度或手动
# data_all$tax  = factor(data_all$tax, levels=rownames(mean_sort))   


########################################################
#从做好的物种分类信息出发
tax_sample = read.table("result/tax/stackplot_p_group.txt", comment.char="", row.names= 1, header=T, sep="\t") 
tax_sample_t = as.data.frame(t(tax_sample))

# 标准化原始reads count为百分比
tax_sample1 = t(t(tax_sample)/colSums(tax_sample,na=T)) * 100 # normalization to total 100

tax_sample=as.data.frame(tax_sample1)
# 数据转换长表格并绘图
tax_sample$tax = rownames(tax_sample)
data_all = as.data.frame(melt(tax_sample, id.vars=c("tax")))
# 设置显示顺序，否则按字母排序(11会排在1后面)
data_all$variable  = factor(data_all$variable, levels=c("Mar","Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec","Jan","Feb", ordered=TRUE))
#data_all$variable  = factor(data_all$variable, levels=c("Intercity", "UrbanHub", "Residential", ordered=TRUE))

#绘制分组柱形图
main_theme = theme(panel.background=element_blank(),
                   panel.grid=element_blank(),
                   axis.line.x=element_line(size=0.5, colour="black"),
                   axis.line.y=element_line(size=0.5, colour="black"),
                   axis.ticks=element_line(color="black"),
                   axis.text=element_text(color="black", size=12),
                   legend.position="right",
                   legend.background=element_blank(),
                   legend.key=element_blank(),
                   legend.text= element_text(size=12),
                   text=element_text(family="sans", size=12),
                   plot.title=element_text(hjust = 0.5,vjust=0.5,size=12),
                   plot.subtitle=element_text(size=12))

p = ggplot(data_all, aes(x=variable, y = value, fill = tax )) + 
  geom_bar(stat = "identity",position="fill", width=0.7)+ 
  scale_y_continuous(labels = scales::percent) + 
  xlab("Groups")+ylab("Percentage (%)")+main_theme+ theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#这条命令是按自己颜色进行设置
scale_fill_manual(values=colours[1:12])
#theme_light(base_line_size = 0)
# 保存pdf和png格式方便查看和编辑
ggsave(paste(opts$output, "group.pdf", sep=""), p, width = opts$width, height = opts$height)
# ggsave(paste(opts$output, "_group.png", sep=""), p, width = opts$width, height = opts$height)
print(paste(opts$output, "_group.pdf/txt finished.", sep = ""))
write.csv(mean_sort,"result/tax/mean_sort.csv")

# 绘制冲击图alluvium
p = ggplot(data = data_all, aes(x = variable, y = value,weight = value, alluvium = tax)) +
  geom_alluvium(aes(fill = tax), alpha = 0.75,decreasing = TRUE,width = 1/2) +
  
  ##去掉边界线
  main_theme + 
  theme(
    panel.border = element_blank(),
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  theme_classic() + 
##添加右边框线
#theme_bw(base_line_size = 1.05,base_rect_size = 1.05)+
scale_fill_manual(values=colours[1:15])+
#facet_wrap(~group, scales = 'free_x', ncol = 2) +  #分面图
ggtitle("Genus changes among groups")
ggsave(paste(opts$output, "_group_phylumpro.pdf", sep=""), p, width = opts$width, height = opts$height)

##绘制连线图
Palette <- c("#ebe7e4", "#d1c7ae", "#c4c7b4","#e5d9da" , "#ecd6bc","#c1aea1", "#c1c2d2","#c6d7e3","#8d9fbb","#48597e", "#FFCC99","#957064", "#E69F00", "#8F7C00")
p <- ggplot(data_all,aes(x = variable, y = value, alluvium = tax, stratum = tax))
p1 <- p + geom_alluvium(aes(fill = tax),alpha = .5,width = 0.5, decreasing = TRUE) + 
  geom_stratum(aes(fill = tax),width = 0.5,color ="NA", decreasing = TRUE)+
  #在各个区块中添加文字标签
  geom_text(stat = 'stratum', aes(label = after_stat(stratum)),infer.label = TRUE, size = 1.6, color="white", decreasing = TRUE)
p2 <- p1 + ylab(label = "Relative Abundance") + xlab(label = "")
p3 <- p2 + scale_fill_manual(values = Palette[1:15])
#调整绘图区主题。
p4 <- p3 + theme_bw()+ theme(panel.grid=element_blank()) + 
  theme(panel.border = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1,colour = "black",size = 12,face = "bold"),
        panel.grid.major = element_blank(),
        ) +
    
  theme(panel.background=element_rect(fill='transparent', 
                                      color='black'),plot.margin = unit(c(3,5,1,1),"mm"))
#调整坐标轴文字。
p5 <- p4 +  
  theme(axis.text.y=element_text(colour = "black",size = 18)) + 
  theme(axis.line = element_line(colour = "black"))+ 
  theme(axis.title.y = element_text(size = 12,face = "bold",
                                    margin = unit(c(0,1,0,1),"lines")))
  ##添加右边框线
  #theme_bw(base_line_size = 1.05,base_rect_size = 1.05)+

#去除横纵坐标与边框之间的空白。
#p6 <- p5 + scale_y_continuous(limits = c(0,1),expand = c(0,0))
ggsave(paste(opts$output, "_Connected_alluvium.pdf", sep=""), p5, width = opts$width, height = opts$height)

