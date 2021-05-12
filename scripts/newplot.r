datapath <- "/Users/pveltsos/git/dpse_mating/output/dedata/rtract"
kdata <- read.table(file.path(datapath, "LogCPM_0.05_rtract.txt"), header=T)
str(kdata)
kdata$cpmM <- (kdata$M_rt_1 + kdata$M_rt_2 + kdata$M_rt_3 + kdata$M_rt_4)/4
kdata$cpmE <- (kdata$E_rt_1 + kdata$E_rt_2 + kdata$E_rt_3 + kdata$E_rt_4)/4
kdata$cpmEE <- (kdata$EE_rt_1 + kdata$EE_rt_2 + kdata$EE_rt_3 + kdata$EE_rt_4)/4
kdata$cpmMM <- (kdata$MM_rt_1 + kdata$MM_rt_2 + kdata$MM_rt_3 + kdata$MM_rt_4)/4

par(mfrow=c(2,1)) 
par(mar=c(5,5,4,3))
plot(kdata$cpmM, kdata$cpmE, xlab="M counts", ylab="E counts", main="E vs M counts", cex.main=1.8, cex.lab=1.3, col='gray80')
points(kdata$cpmM[kdata$FDR.EvsM<0.05], kdata$cpmE[kdata$FDR.EvsM<0.05], pch=20, col=4)
points(kdata$cpmM[kdata$FDR.EEvsMM<0.05], kdata$cpmE[kdata$FDR.EEvsMM<0.05], pch=4, col=2)
legend("topleft", inset=0.05, legend=c("EvsM", "EEvsMM"), pch =c(20,4), col=c(4,2) ) 
abline(coef = c(0,1), col='darkgreen')

plot(kdata$cpmM + kdata$cpmE, kdata$cpmMM + kdata$cpmEE, xlab="virgin counts", ylab="mated counts", main="E+M vs EE+MM counts", cex.main=1.8, cex.lab=1.3, col='gray80')
points(kdata$cpmM[kdata$FDR.E.MMvsM.EE<0.05], kdata$cpmE[kdata$FDR.E.MMvsM.EE<0.05], pch=20, col=4)
points(kdata$cpmM[kdata$FDR.E.MvsEE.MM<0.05], kdata$cpmE[kdata$FDR.E.MvsEE.MM<0.05], pch=4, col=2)
legend("topleft", inset=0.05, legend=c("E.MMvsM.EE", "E.MvsEE.MM"), pch =c(20,4), col=c(4,2) ) 
abline(coef = c(0,1), col='darkgreen')


plot(kdata$logFC.E.MvsEE.MM, kdata$logFC.E.MMvsM.EE, xlab="virgin vs mated logFC", ylab="interaction logFC", main="virgin/mated vs interaction", cex.main=1.8, cex.lab=1.3, col='gray80')
points(kdata$logFC.E.MvsEE.MM[kdata$FDR.E.MMvsM.EE<0.05], kdata$logFC.E.MMvsM.EE[kdata$FDR.E.MMvsM.EE<0.05], pch=20, col=4)
points(kdata$logFC.E.MvsEE.MM[kdata$FDR.E.MvsEE.MM<0.05], kdata$logFC.E.MMvsM.EE[kdata$FDR.E.MvsEE.MM<0.05], pch=4, col=2)
legend("topleft", inset=0.05, legend=c("E.MMvsM.EE", "E.MvsEE.MM"), pch =c(20,4), col=c(4,2) ) 
abline(h=0, col='darkgreen')
abline(v=0, col='darkgreen')


