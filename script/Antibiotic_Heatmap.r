###heatmap###
library(pheatmap)
library(ggplot2)
library(reshape2)
otu.tab <- read.table(file = "result/stackplot_mechenism_Month_group.txt", sep = "\t", head = TRUE, row.names = 1)
##将表达量为0的值设为NA
otu.tab[otu.tab==0]=NA
##将NA值数据库中最小值*0.01
otu.tab[is.na(otu.tab)]=min(otu.tab, na.rm =T)*0.01

# 假设我们想要按照以下顺序对行进行排序  
desired_order <- c(	"Jan", "Feb",	"Mar","Apr","May","Jun",	"Jul",	"Aug","Sep",	"Oct",	"Nov",	"Dec")  
# 重新排序数据  
otu.tab_sorted <- otu.tab[, desired_order]  


pheatmap(log2(otu.tab_sorted+1),cluster_row=T,clustering_method = "average", 
         color = colorRampPalette(c("navy","white","firebrick3"))(50),
         cluster_cols = F,scale="row",border_color = "NA",filename = "result/Anti_heatmap.pdf")



otu.tab <- read.table(file = "result/Residential_Month_group.txt", sep = "\t", head = TRUE, row.names = 1)
##将表达量为0的值设为NA
otu.tab[otu.tab==0]=NA
##将NA值数据库中最小值*0.01
otu.tab[is.na(otu.tab)]=min(otu.tab, na.rm =T)*0.01

# 假设我们想要按照以下顺序对行进行排序  
desired_order <- c(	"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")  
# 重新排序数据  
otu.tab_sorted <- otu.tab[, desired_order]  


pheatmap(log2(otu.tab_sorted+1),cluster_row=T,clustering_method = "average", 
         color = colorRampPalette(c("navy","white","firebrick3"))(50),
         cluster_cols = F,scale="row",border_color = "NA",filename = "result//Residential_Anti_heatmap.pdf")
##

otu.tab <- read.table(file = "result/UrbanHub_Month_group.txt", sep = "\t", head = TRUE, row.names = 1)
##将表达量为0的值设为NA
otu.tab[otu.tab==0]=NA
##将NA值数据库中最小值*0.01
otu.tab[is.na(otu.tab)]=min(otu.tab, na.rm =T)*0.01

# 假设我们想要按照以下顺序对行进行排序  
desired_order <- c(	"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")  
# 重新排序数据  
otu.tab_sorted <- otu.tab[, desired_order]  


pheatmap(log2(otu.tab_sorted+1),cluster_row=T,clustering_method = "average", 
         color = colorRampPalette(c("navy","white","firebrick3"))(50),
         cluster_cols = F,scale="row",border_color = "NA",filename = "result/UrbanHub_Anti_heatmap.pdf")


otu.tab <- read.table(file = "result/Intercity_Month_group.txt", sep = "\t", head = TRUE, row.names = 1)
##将表达量为0的值设为NA
otu.tab[otu.tab==0]=NA
##将NA值数据库中最小值*0.01
otu.tab[is.na(otu.tab)]=min(otu.tab, na.rm =T)*0.01

# 假设我们想要按照以下顺序对行进行排序
desired_order <- c(	"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")  
# 重新排序数据  
otu.tab_sorted <- otu.tab[, desired_order]


pheatmap(log2(otu.tab_sorted+1),cluster_row=T,clustering_method = "average", 
         color = colorRampPalette(c("navy","white","firebrick3"))(50),
         cluster_cols = F,scale="row",border_color = "NA",filename = "result/Intercity_Anti_heatmap.pdf")
