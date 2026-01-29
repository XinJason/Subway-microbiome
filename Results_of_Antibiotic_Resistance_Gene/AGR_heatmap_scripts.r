
# 安装并加载必要的包
if(!require(pheatmap)) install.packages("pheatmap")
if(!require(RColorBrewer)) install.packages("RColorBrewer")
library(pheatmap)
library(RColorBrewer)

# 读取数据
# 读取数据
data <- read.table("ARG_Heatmap_matrix1.txt", header = TRUE, row.names=1, sep = "\t", stringsAsFactors = FALSE)

# 处理数据：转换为数值矩阵并进行log1p转换
data_matrix <- as.matrix(data)
data_log <- log1p(data_matrix)

# 设置配色
color_range <- colorRampPalette(c("white", "firebrick4"))(50)

# 绘制热图并保存为PDF
pdf("antibiotic_resistance_heatmap_with_borders1.pdf", width = 12, height = 13)

pheatmap(data_log,
         col = color_range,
         border_color = "gray50",  # 显示灰色边框（可改为black、gray等）
         fontsize = 8,
         fontsize_row = 10,
         fontsize_col = 10,
         scale = "none",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         legend = TRUE,
         legend_breaks = c(0, max(data_log, na.rm = TRUE)/2, max(data_log, na.rm = TRUE)),
         treeheight_row = 0,
         treeheight_col = 0,
         angle_col = 90,
         labels_row = rownames(data),
         labels_col = colnames(data),
         hjust = 1,
         vjust = 1,
         angle_row = 0
)

dev.off()
