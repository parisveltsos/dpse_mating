# Setup
library(zoo)

library(AID)
library(MASS)
library(car)
library(sjPlot)
library(ggridges)
library(ggplot2)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
library(grid) 
theme_set(theme_bw(base_size=20)) ## No gray backgroundcbred <- 2 # '#D55E00'

# Input data

datapath <- '~/git/dpse_mating/input/ridgeplots/'
outpath <- '~/git/dpse_mating/output/ridgeplots/'

mating_data <- read.table(file.path(datapath, 'mating.txt'), header=T)
str(mating_data) 

ridge_logfc <- ggplot(mating_data, aes(x=abs(logfc), y=mating, fill = mating)) + geom_density_ridges(aes(point_fill = mating), alpha = .6) + xlab("abs(log2FC)") + scale_fill_manual(values=c("gray", "gray20", "gray", "gray20", "gray", "gray20", "gray20", "gray", "gray20")) + theme(legend.position="none") + coord_cartesian(ylim = c(1.3, 10.5)) + theme(axis.text.x = element_text(size=26), text = element_text(size=34)) + 

annotate(geom="text", x=8.7, y=1.4, label=nrow(subset(mating_data, mating_data$mating=='ov_E_mated')), size=10, hjust=0) +
annotate(geom="text", x=9.2, y=2.4, label=nrow(subset(mating_data, mating_data$mating=='ov_E_virgin')), size=10, hjust=0) +
annotate(geom="text", x=8, y=6.4, label=nrow(subset(mating_data, mating_data$mating=='rt_E_mated')), size=10, hjust=0) +
annotate(geom="text", x=8, y=7.4, label=nrow(subset(mating_data, mating_data$mating=='rt_E_virgin')), size=10, hjust=0) +
annotate(geom="text", x=8, y=8.4, label=nrow(subset(mating_data, mating_data$mating=='rt_M_mated')), size=10, hjust=0) +
annotate(geom="text", x=8, y=9.4, label=nrow(subset(mating_data, mating_data$mating=='rt_M_virgin')), size=10, hjust=0)

dev.copy(pdf,file.path(outpath,'figS1.pdf'), width=16, height=12) # Fig S1
ridge_logfc
dev.off()

wilcox.test(abs(mating_data$logfc[mating_data$mating=="rt_M_virgin"]), abs(mating_data$logfc[mating_data$mating=="rt_M_mated"]))
wilcox.test(abs(mating_data$logfc[mating_data$mating=="rt_E_virgin"]), abs(mating_data$logfc[mating_data$mating=="rt_E_mated"]))
wilcox.test(abs(mating_data$logfc[mating_data$mating=="ov_E_virgin"]), abs(mating_data$logfc[mating_data$mating=="ov_E_mated"]))


mating_data2 <- read.table(file.path(datapath, 'mating2.txt'), header=T)
str(mating_data2) 

ridge_logfc <- ggplot(mating_data2, aes(x=abs(logfc), y=mating, fill = mating)) + geom_density_ridges(aes(point_fill = mating), alpha = .6) + xlab("abs(log2FC)") + scale_fill_manual(values=c("gray", "gray20", "gray", "gray20", "gray", "gray20", "gray", "gray20")) + theme(legend.position="none") + coord_cartesian(ylim = c(1.3, 10.5)) + theme(axis.text.x = element_text(size=26), text = element_text(size=34)) + 

annotate(geom="text", x=6.7, y=1.4, label=nrow(subset(mating_data2, mating_data2$mating=='ov_1Mfem_x_Emale')), size=10, hjust=0) +
annotate(geom="text", x=6.7, y=2.4, label=nrow(subset(mating_data2, mating_data2$mating=='ov_2Efem_x_Emale')), size=10, hjust=0) +
annotate(geom="text", x=6.7, y=3.4, label=nrow(subset(mating_data2, mating_data2$mating=='ov_3Efem_x_Mmale')), size=10, hjust=0) +
annotate(geom="text", x=6.7, y=4.4, label=nrow(subset(mating_data2, mating_data2$mating=='ov_4Mfem_x_Mmale')), size=10, hjust=0) +

annotate(geom="text", x=6.7, y=6.4, label=nrow(subset(mating_data2, mating_data2$mating=='rt_1Mfem_x_Emale')), size=10, hjust=0) +
annotate(geom="text", x=6.7, y=7.4, label=nrow(subset(mating_data2, mating_data2$mating=='rt_2Efem_x_Emale')), size=10, hjust=0) +
annotate(geom="text", x=6.7, y=8.4, label=nrow(subset(mating_data2, mating_data2$mating=='rt_3Efem_x_Mmale')), size=10, hjust=0) +
annotate(geom="text", x=6.7, y=9.4, label=nrow(subset(mating_data2, mating_data2$mating=='rt_4Mfem_x_Mmale')), size=10, hjust=0)


dev.copy(pdf,file.path(outpath,'figS3.pdf'), width=17.2, height=12)
ridge_logfc
dev.off()

wilcox.test(abs(mating_data2$logfc[mating_data2$mating=="rt_3Efem_x_Mmale"]), abs(mating_data2$logfc[mating_data2$mating=="rt_4Mfem_x_Mmale"]))
wilcox.test(abs(mating_data2$logfc[mating_data2$mating=="rt_1Efem_x_Emale"]), abs(mating_data2$logfc[mating_data2$mating=="rt_2Mfem_x_Emale"]))
wilcox.test(abs(mating_data2$logfc[mating_data2$mating=="ov_3Efem_x_Mmale"]), abs(mating_data2$logfc[mating_data2$mating=="ov_4Mfem_x_Mmale"]))
wilcox.test(abs(mating_data2$logfc[mating_data2$mating=="ov_2Efem_x_Emale"]), abs(mating_data2$logfc[mating_data2$mating=="ov_1Mfem_x_Emale"]))

mating_data3 <- read.table(file.path(datapath, 'mating3.txt'), header=T)
str(mating_data3) 

ridge_logfc <- ggplot(mating_data3, aes(x=abs(logfc), y=mating, fill = mating)) + geom_density_ridges(aes(point_fill = mating), alpha = .6) + xlab("abs(log2FC)") + scale_fill_manual(values=c("gray", "gray20", "gray", "gray20", "gray", "gray20", "gray", "gray20")) + theme(legend.position="none") + coord_cartesian(ylim = c(1.3, 10.5)) + theme(axis.text.x = element_text(size=26), text = element_text(size=34)) + 

annotate(geom="text", x=6.7, y=1.4, label=nrow(subset(mating_data3, mating_data3$mating=='1ov_Efem_x_Mmale')), size=10, hjust=0) +
annotate(geom="text", x=6.7, y=2.4, label=nrow(subset(mating_data3, mating_data3$mating=='2ov_Efem_x_Emale')), size=10, hjust=0) +
annotate(geom="text", x=6.7, y=3.4, label=nrow(subset(mating_data3, mating_data3$mating=='ov_Mfem_x_Emale')), size=10, hjust=0) +
annotate(geom="text", x=6.7, y=4.4, label=nrow(subset(mating_data3, mating_data3$mating=='ov_Mfem_x_Mmale')), size=10, hjust=0) +

annotate(geom="text", x=6.7, y=7.4, label=nrow(subset(mating_data3, mating_data3$mating=='rt_1Efem_x_Mmale')), size=10, hjust=0) +
annotate(geom="text", x=6.7, y=6.4, label=nrow(subset(mating_data3, mating_data3$mating=='rt_2Efem_x_Emale')), size=10, hjust=0) +
annotate(geom="text", x=6.7, y=8.4, label=nrow(subset(mating_data3, mating_data3$mating=='rt_3Mfem_x_Emale')), size=10, hjust=0) +
annotate(geom="text", x=6.7, y=9.4, label=nrow(subset(mating_data3, mating_data3$mating=='rt_4Mfem_x_Mmale')), size=10, hjust=0) 

dev.copy(pdf,file.path(outpath,'figS2.pdf'), width=17.2, height=12)
ridge_logfc
dev.off()

wilcox.test(abs(mating_data3$logfc[mating_data3$mating=="ov_Efem_x_Mmale"]), abs(mating_data3$logfc[mating_data3$mating=="ov_Efem_x_Emale"]))


# The following figures are not presented
## Virgin males
# 
# mating_data4 <- read.table(file.path(datapath, 'mating4.txt'), header=T)
# str(mating_data4) 
# 
# ridge_logfc <- ggplot(mating_data4, aes(x=abs(logfc), y=males, fill = males)) + geom_density_ridges(aes(point_fill = males), alpha = .6) + xlab("abs(log2FC)") + scale_fill_manual(values=c("gray", "gray20", "gray", "gray20")) + theme(legend.position="none") + coord_cartesian(ylim = c(1.3, 5.5)) + theme(axis.text.x = element_text(size=26), text = element_text(size=34)) + 
# 
# annotate(geom="text", x=8.7, y=1.4, label=nrow(subset(mating_data4, mating_data4$males=='ag_E_virgin')), size=10, hjust=0) +
# annotate(geom="text", x=8.7, y=2.4, label=nrow(subset(mating_data4, mating_data4$males=='ag_M_virgin')), size=10, hjust=0) +
# annotate(geom="text", x=8.3, y=3.4, label=nrow(subset(mating_data4, mating_data4$males=='test_E_virgin')), size=10, hjust=0) +
# annotate(geom="text", x=8.3, y=4.4, label=nrow(subset(mating_data4, mating_data4$males=='test_M_virgin')), size=10, hjust=0)
# 
# dev.copy(pdf,file.path(datapath,'figS4.pdf'), width=16, height=12)
# ridge_logfc
# dev.off()
# 
# wilcox.test(abs(mating_data4$logfc[mating_data4$males=="ag_E_virgin"]), abs(mating_data4$logfc[mating_data4$males=="ag_M_virgin"]))
# wilcox.test(abs(mating_data4$logfc[mating_data4$males=="test_E_virgin"]), abs(mating_data4$logfc[mating_data4$males=="test_M_virgin"]))



## Virgin females

# mating_data5 <- read.table(file.path(datapath, 'mating5.txt'), header=T)
# str(mating_data5) 
# 
# ridge_logfc <- ggplot(mating_data5, aes(x=abs(logfc), y=females, fill = females)) + geom_density_ridges(aes(point_fill = females), alpha = .6) + xlab("abs(log2FC)") + scale_fill_manual(values=c("gray", "gray20", "gray", "gray20")) + theme(legend.position="none") + coord_cartesian(ylim = c(1.3, 5.5)) + theme(axis.text.x = element_text(size=26), text = element_text(size=34)) + 
# 
# annotate(geom="text", x=8.7, y=1.4, label=nrow(subset(mating_data5, mating_data5$females=='frt_E_virgin')), size=10, hjust=0) +
# annotate(geom="text", x=8.7, y=2.4, label=nrow(subset(mating_data5, mating_data5$females=='frt_M_virgin')), size=10, hjust=0) +
# annotate(geom="text", x=8.3, y=3.4, label=nrow(subset(mating_data5, mating_data5$females=='ov_E_virgin')), size=10, hjust=0) +
# annotate(geom="text", x=8.3, y=4.4, label=nrow(subset(mating_data5, mating_data5$females=='ov_M_virgin')), size=10, hjust=0)
# 
# dev.copy(pdf,file.path(datapath,'figS5.pdf'), width=16, height=12)
# ridge_logfc
# dev.off()
# 
# wilcox.test(abs(mating_data5$logfc[mating_data5$females=="frt_E_virgin"]), abs(mating_data5$logfc[mating_data5$females=="frt_M_virgin"]))
# wilcox.test(abs(mating_data5$logfc[mating_data5$females=="ov_E_virgin"]), abs(mating_data5$logfc[mating_data5$females=="ov_M_virgin"]))
