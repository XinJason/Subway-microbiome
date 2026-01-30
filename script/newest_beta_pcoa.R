rm(list = ls())
# 2.1 安装CRAN来源常用包
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 依赖包列表：参数解析、数据变换、绘图和开发包安装、安装依赖、ggplot主题
package_list = c("reshape2","ggplot2","vegan","cowplot")
# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}
# 2.2 安装bioconductor常用包
package_list = c("digest","ggrepel","ggsci")
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    source("https://bioconductor.org/biocLite.R")
    biocLite(p)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}
# 2.3 安装Github常用包
# 参数解析、数据变换、绘图和开发包安装
package_list = c("kassambara/ggpubr")
for(p in package_list){
  q=unlist(strsplit(p,split = "/"))[2]
  if(!suppressWarnings(suppressMessages(require(q, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install_github(p)
    suppressWarnings(suppressMessages(library(q, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}
options(warn = -1)
# 1.2 解析命令行
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
    make_option(c("-t", "--type"), type="character", default="unifrac",
                help="Distance type; 距离类型, 可选bray_curtis, bray_curtis_binary, euclidean, jaccard, jaccard_binary, manhatten, unifrac, unifrac_binary [default %default]"),   
    make_option(c("-i", "--input"), type="character", default="",
                help="Input beta distance; 距离矩阵,默认beta目录下与t同名，可指定 [default %default]"),
    make_option(c("-d", "--design"), type="character", default="result/metadata.txt",
                help="design file; 实验设计文件 [default %default]"),
    make_option(c("-n", "--group"), type="character", default="Seasonal",
                help="name of group type; 分组列名 [default %default]"),
    make_option(c("-w", "--width"), type="numeric", default=7,
                help="Width of figure; 图片宽 [default %default]"),
    make_option(c("-e", "--height"), type="numeric", default=5,
                help="Height of figure; 图片高 [default %default]"),
    make_option(c("-o", "--output"), type="character", default="",
                help="output directory or prefix; 输出文件前缀, 有txt和矢量图pdf [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))
  
  # 调置如果无调设置输出，根据其它参数设置默认输出
  if (opts$input==""){opts$input=paste("result/beta/",opts$type, ".txt", sep = "")}
  if (opts$output==""){opts$output=paste("result/beta/pcoa_",opts$type, sep = "")}
  
  # 显示输入输出确认是否正确
  print(paste("The distance matrix file is ", opts$input,  sep = ""))
  print(paste("Type of distance type is ", opts$type,  sep = ""))
  print(paste("The design file is ", opts$design,  sep = ""))
  print(paste("The group name is ", opts$group,  sep = ""))
  print(paste("The output file prefix is ", opts$output, sep = ""))
}
# 3. 读取输入文件
# 读取距离矩阵文件
dis = read.table(opts$input, header=T, row.names= 1, sep="\t", comment.char="") 
# 读取实验设计
design = read.table(opts$design, header=T, row.names= 1, sep="\t", comment.char="") 
# 提取样品组信息,默认为group可指定
# sampFile = as.data.frame(design[,opts$group],row.names = row.names(design))
# colnames(sampFile)[1] = "group" # 单列数据框筛选会变为list，改为双列
sampFile = data.frame(group=design[,opts$group],
                      sample=row.names(design), 
                      row.names = row.names(design))
# 数据筛选，筛选两文件中共有
idx = rownames(sampFile) %in% colnames(dis) # match design with alpha
sampFile = sampFile[idx,]
dis = dis[rownames(sampFile),rownames(sampFile)] 
design =design[rownames(sampFile),] 
# 4. 统计与绘图
# vegan:cmdscale计算矩阵矩阵中主坐标轴坐标，取前3维
pcoa = cmdscale(dis, k=3, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
points = as.data.frame(pcoa$points) # get coordinate string, format to dataframme
##summary(points)
eig = pcoa$eig
points = cbind(points, sampFile[rownames(points),])
colnames(points) = c("x", "y", "z","group") 
#各 PCoA 轴的解释量
pcoa_exp <- pcoa$eig/sum(pcoa$eig)
##查看PCOA各轴解释率

head(pcoa_exp)
#"#c8ccc1",  "#aabdc4", "#8d958f","#8F7C00","#FFCC99", "#E69F00", "#cf9198", "#dB9c7c","#48597e"，"#8d9fbb"
site <- scores(pcoa)
#地铁细菌颜色 colour <- c("#8d9fbb","#CD950C","#104E8B","#48597e")
colour <- c("#7CCD7C", "#1E90FF", "#CD853F", "#DCDCDC")
#class.otu.2 <- melt(t(class.tab.otu), measure.vars = colnames(class.tab.otu))
#class.seq.2 <- melt(t(class.tab.seq), measure.vars = colnames(class.tab.seq))
##自己加的，调整分组顺序命令
points$group <- factor(points$group,levels=c("Spring", "Summer", "Autumn", "Winter", ordered=TRUE))
# plot PCo 1 and 2
## 在aes中添加颜色命令aes(x=x, y=y, color=design$group, shape =design$group2)
p = ggplot(points, aes(x=x, y=y, color=design$group, shape =design$StationType)) + geom_point(alpha=0.7, size=4) +
    
#stat_ellipse(level = 0.9, aes(fill = design$group),linetype = 2, type = "norm") +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       title=paste(opts$type," PCoA",sep=""))  + 

#这条命令是按自己颜色进行设置
#scale_fill_manual(values=colours[1:12])+  
#geom用于设置填充形状，alpha设置透明度。不设置则为实心填充，遮盖椭圆中的点, levels设置confidence ellipses的置信区间, 在0-1范围内。levels越小椭圆面积越小，涵盖的点越集中。
# 不需要填充的时候去掉fill及对fill的补充参数geom,alpha等即可
#theme_classic()+
#scale_color_brewer(palette="Set1")+
  #scale_fill_manual(values=colorRampPalette(c("#8d9fbb","#CD950C","#104E8B","#48597e"))(4))+
#geom_point(aes(shape = design$Group1))
theme_bw()+
theme(axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=18,face = "bold",vjust = 1.5),
        axis.title.y=element_text(colour='black', size=18,face = "bold",vjust = 1.5),
        axis.text.y=element_text(colour='black',size=15),
        axis.text.x=element_text(colour = "black",size = 15,
                                 hjust = 1,vjust = 0.5)) +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +

scale_color_manual(values = colour)
# 颜色可以自己设置，或者直接用scale_color_brewer()
#p1=p+ scale_color_npg()  
#annotate("text",0.15,0.45,label = "P=0.211***
#                                   Group:R2=0.195***",size=5)  
# 保存pdf和png格式方便查看和编辑
ggsave(paste(opts$output, ".pdf", sep=""), p, width = opts$width, height = opts$height)
#图片组合
#g=cowplot::plot_grid(p2,p2,p2,p2,ncol= 2, rel_widths = c(1, 1), labels="auto",label_size = 12)
#ggsave(paste(opts$output, "Combined_Beta.pdf", sep=""), g, width = opts$width, height = opts$height)

# Compare each group beta by vegan adonis in bray_curtis
da_adonis = function(sampleV){
  sampleA = as.matrix(sampleV$sampA)
  sampleB = as.matrix(sampleV$sampB)
  design2 = subset(sampFile, group %in% c(sampleA,sampleB))
  if (length(unique(design2$group))>1) {
    sub_dis_table = dis_table[rownames(design2),rownames(design2)]
    sub_dis_table = as.dist(sub_dis_table, diag = FALSE, upper = FALSE)
    adonis_table = adonis(sub_dis_table~group, data=design2, permutations = 10000) 
    adonis_pvalue = adonis_table$aov.tab$`Pr(>F)`[1]
    print(paste("In ",opts$type," pvalue between", sampleA, "and", sampleB, "is", adonis_pvalue, sep=" "))
    adonis_pvalue = paste(opts$type, sampleA, sampleB, adonis_pvalue, sep="\t")
    write.table(adonis_pvalue, file=paste(opts$output, ".txt", sep=""), append = TRUE, sep="\t", quote=F, row.names=F, col.names=F)
  }
}

# loop for each group pair
dis_table = as.matrix(dis)
if (TRUE) {
  compare_data = as.vector(unique(design[[opts$group]]))
  len_compare_data = length(compare_data)
  for(i in 1:(len_compare_data-1)) {
    for(j in (i+1):len_compare_data) {
      tmp_compare = as.data.frame(cbind(sampA=compare_data[i],sampB=compare_data[j]))
      print(tmp_compare)
      da_adonis(tmp_compare)
    }
  }
}else {
  compare_data = read.table("doc/compare.txt", sep="\t", check.names=F, quote='', comment.char="")
  colnames(compare_data) = c("sampA", "sampB")
  for(i in 1:dim(compare_data)[1]){da_adonis(compare_data[i,])}
}	 
print(paste("Adnois statistics result in",opts$output, ".txt is finished.", sep = ""))

# 5. 保存图表
# 提示工作完成
print(paste("Output in ", opts$output, ".txt/pdf finished.", sep = ""))
##计算群落之间的相似性
#若是已经提供好了距离矩阵，则直接使用现有的距离矩阵进行分析即可
anosim_result_dis <- anosim(dis, design$Station_Type, permutations = 999) 
###查看结果
summary(anosim_result_dis)
anosim_result_dis
names(anosim_result_dis)
#若是使用 OTU 丰度表，则需要在计算时指定所依据的距离类型，这里依然使用 Bray-Curtis 距离
#anosim_result_otu <- anosim(otu, group$site, permutations = 999, distance = 'bray')        
#dis1 <- vegdist(otu, method = 'bray')
#anosim_result_dis1 <- anosim(dis1, group$site, permutations = 999)   #同上所述
##读入文件
#现有的距离矩阵
dis <- read.delim('result/beta/bray_curtis.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
#dis <- as.dist(dis)    #将导入的样本间距离转化为 dist 类型
#或者直接使用 OTU 丰度表
otu <- read.delim('result/otutab_rare.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu <- data.frame(t(otu))
# 读取实验设计
group = read.table("result/metadata.txt", header=T, row.names= 1, sep="\t", comment.char="") 

# 交叉筛选 OTU和design共有，三者共有
idx = rownames(otu) %in% rownames(group)
otu = otu[idx, ]

dis =dis[rownames(group),]
design = group[rownames(dis),]
#ANOSIM 分析（所有分组间比较，即整体差异）
#（1）若是已经提供好了距离矩阵，则直接使用现有的距离矩阵进行分析即可
anosim_result_dis <- anosim(dis, design$group, permutations = 999)       #根据 group$site 这一列样本分组信息进行 ANOSIM 分析，随机置换检验 999 次
anosim_result_dis
summary(anosim_result_dis)
#（2）若是使用 OTU 丰度表，则需要在计算时指定所依据的距离类型，这里依然使用 Bray-Curtis 距离
anosim_result_otu <- anosim(otu, group$group, permutations = 999, distance = 'bray')         
#（3）或者首先根据丰度表计算样本距离，在将所的距离数据作为输入
dis1 <- vegdist(otu, method = 'bray')
anosim_result_dis1 <- anosim(dis1, group$group, permutations = 999)   #同上所述
#查看结果，上述 3 条命令所计算的内容一致
#或者names(anosim_result_dis1)
#可使用plot()命令简单展示统计结果。
#作图展示
#pdf(paste('anosim.all.pdf', sep = ''), width = 10, height = 5)
pdf(file = "anosim.locations.pdf",width = 10)
plot(anosim_result_dis1, col = c("#D55E00","#426600","#66CC99", "#009E73", "#0072B2","#990000","#FF0010","#2BCE48"))
dev.off()
