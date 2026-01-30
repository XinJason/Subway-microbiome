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


########################################################
#从做好的物种分类信息出发
tax_sample = read.table("Feast_result_month.txt", comment.char="", row.names= 1, header=T, sep="\t") 

# 数据转换长表格并绘图
tax_sample$tax = rownames(tax_sample)
data_all = as.data.frame(melt(tax_sample, id.vars=c("tax")))
# 设置显示顺序，否则按字母排序(11会排在1后面)
colours <- c("#48597e", "#FFCC99","#957064");

data_all$variable  = factor(data_all$variable, levels=c( "Jan", "Feb", "Mar","Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec", ordered=TRUE))
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
  scale_fill_manual(values=colours[1:12])+
  ##添加右边框线
  theme_bw(base_line_size = 1.05,base_rect_size = 1.05)
#theme_light(base_line_size = 0)
# 保存pdf和png格式方便查看和编辑
ggsave("Feast_Station_month.pdf", p, width = 8, height = 5)
# ggsave(paste(opts$output, "_group.png", sep=""), p, width = opts$width, height = opts$height)


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
  theme_bw(base_line_size = 1.05,base_rect_size = 1.05)+
  scale_fill_manual(values=colours[1:15])+
  #facet_wrap(~group, scales = 'free_x', ncol = 2) +  #分面图
  ggtitle("Genus changes among groups")
ggsave("Feast_Station_month_group_phylumpro.pdf", p, width = 8, height = 5)

##绘制连线图
Palette <- c("#48597e", "#FFCC99","#957064")
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
                                    margin = unit(c(0,1,0,1),"lines")))+
##添加右边框线
theme_bw(base_line_size = 1.05,base_rect_size = 1.05)

#去除横纵坐标与边框之间的空白。
#p6 <- p5 + scale_y_continuous(limits = c(0,1),expand = c(0,0))
ggsave("Feast_Station_month_Connected_alluvium.pdf",  p5, width = 8, height = 5)

