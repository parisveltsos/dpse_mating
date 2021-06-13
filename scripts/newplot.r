datapath <- "~/Library/Mobile Documents/com~apple~CloudDocs/MyDocuments/Work/Projects/gitOut/dpseMating/dedata/rtract"
kdata <- read.table(file.path(datapath, "LogCPM_0.05_rtract.txt"), header=T)
str(kdata)
kdata$cpmM <- (kdata$M_rt_1 + kdata$M_rt_2 + kdata$M_rt_3 + kdata$M_rt_4)/4
kdata$cpmE <- (kdata$E_rt_1 + kdata$E_rt_2 + kdata$E_rt_3 + kdata$E_rt_4)/4
kdata$cpmEE <- (kdata$EE_rt_1 + kdata$EE_rt_2 + kdata$EE_rt_3 + kdata$EE_rt_4)/4
kdata$cpmMM <- (kdata$MM_rt_1 + kdata$MM_rt_2 + kdata$MM_rt_3 + kdata$MM_rt_4)/4

# for ovary
# kdata$cpmM <- (kdata$M_ov_1 + kdata$M_ov_2 + kdata$M_ov_3 + kdata$M_ov_4)/4
# kdata$cpmE <- (kdata$E_ov_1 + kdata$E_ov_2 + kdata$E_ov_3 + kdata$E_ov_4)/4
# kdata$cpmEE <- (kdata$EE_ov_1 + kdata$EE_ov_2 + kdata$EE_ov_3 )/3
# kdata$cpmMM <- (kdata$MM_ov_1 + kdata$MM_ov_2 + kdata$MM_ov_3 + kdata$MM_ov_4)/4


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
points(kdata$logFC.E.MvsEE.MM[kdata$FDR.allEvsM<0.05], kdata$logFC.E.MMvsM.EE[kdata$FDR.allEvsM<0.05], pch=3, col=6)
legend("topleft", inset=0.05, legend=c("E.MMvsM.EE", "E.MvsEE.MM", "all E vs M"), pch =c(20,4,3), col=c(4,2, 6) ) 
abline(h=0, col='darkgreen')
abline(v=0, col='darkgreen')

plot(kdata$logFC.E.MvsEE.MM, kdata$logFC.allEvsM, xlab="virgin vs mated logFC", ylab="E vs M all", main="virgin/mated vs E.EE/M.MM", cex.main=1.8, cex.lab=1.3, col='gray80')
points(kdata$logFC.E.MvsEE.MM[kdata$FDR.E.MMvsM.EE<0.05], kdata$logFC.allEvsM[kdata$FDR.E.MMvsM.EE<0.05], pch=20, col=4)
points(kdata$logFC.E.MvsEE.MM[kdata$FDR.E.MvsEE.MM<0.05], kdata$logFC.allEvsM[kdata$FDR.E.MvsEE.MM<0.05], pch=4, col=2)
points(kdata$logFC.E.MvsEE.MM[kdata$FDR.allEvsM<0.05], kdata$logFC.allEvsM[kdata$FDR.allEvsM<0.05], pch=3, col=6)
legend("topleft", inset=0.05, legend=c("E.MMvsM.EE", "virgin/mated", "all E vs M"), pch =c(20,4,3), col=c(4,2, 6) ) 
abline(h=0, col='darkgreen')
abline(v=0, col='darkgreen')




idata <- subset(kdata, kdata$FDR.E.MMvsM.EE<0.05)
head(idata)
par(mar=c(5,5,4,3))
plot(c(1,2), c(min(idata$logFC.EvsM, idata$logFC.EEvsMM), max(idata$logFC.EvsM, idata$logFC.EEvsMM)), type='n', ylab="logFC E vs M", xlab="", main='E vs M interaction virgin vs mated')

for (i in 1:nrow(idata)) {
	if (idata$logFC.EvsM[i] > idata$logFC.EEvsMM[i]) {
		lines(c(idata$logFC.EvsM[i],idata$logFC.EEvsMM[i]), col=2)
		} else {
			lines(c(idata$logFC.EvsM[i],idata$logFC.EEvsMM[i]), col=4)
		}
	}
legend("topright", inset=0.05, legend=c("virgin > mated", "virgin < mated"), pch =c(19,19), col=c(2,4), cex=0.8 )

par(mar=c(5,5,4,3))
plot(c(1,4), c(min(idata$cpmM, idata$cpmE, idata$cpmMM, idata$cpmEE), max(idata$cpmM, idata$cpmE, idata$cpmMM, idata$cpmEE)), type='n', ylab="normalised count", xlab="", main='E vs M interaction virgin vs mated')

for (i in 1:nrow(idata)) {
	if (idata$cpmE[i] > idata$cpmM[i] & idata$cpmMM[i] > idata$cpmEE[i] & idata$cpmE[i] > idata$cpmMM[i] & (idata$cpmM[i] + idata$cpmEE[i]) < (idata$cpmE[i] + idata$cpmMM[i])) {
		lines(c(idata$cpmM[i],idata$cpmMM[i],idata$cpmE[i],idata$cpmEE[i]), col=2)
	} else if (idata$cpmE[i] > idata$cpmM[i] & idata$cpmMM[i] > idata$cpmEE[i] & idata$cpmE[i] < idata$cpmMM[i] & (idata$cpmM[i] + idata$cpmEE[i]) < (idata$cpmE[i] + idata$cpmMM[i])) {
	lines(c(idata$cpmM[i],idata$cpmMM[i],idata$cpmE[i],idata$cpmEE[i]), col=4)
	}	
		else {
		lines(c(idata$cpmM[i],idata$cpmMM[i],idata$cpmE[i],idata$cpmEE[i]), col=5)
			}
	}
	
par(mfrow=c(1,3)) 	
pat1 <- subset(idata, idata$cpmE > idata$cpmM & idata$cpmMM > idata$cpmEE & idata$cpmE > idata$cpmMM & (idata$cpmM + idata$cpmEE) < (idata$cpmE + idata$cpmMM))	
plot(c(1,4), c(min(idata$cpmM, idata$cpmE, idata$cpmMM, idata$cpmEE), max(idata$cpmM, idata$cpmE, idata$cpmMM, idata$cpmEE)), type='n', ylab="normalised count", xlab="", main=paste('(E + MM) > (M + EE) & E > MM, line number:', nrow(pat1)))
for (i in 1:nrow(pat1)) {
		lines(c(pat1$cpmM[i],pat1$cpmMM[i],pat1$cpmE[i],pat1$cpmEE[i]), col=2)
	}
	
pat2 <- subset(idata, idata$cpmE > idata$cpmM & idata$cpmMM > idata$cpmEE & idata$cpmE < idata$cpmMM & (idata$cpmM + idata$cpmEE) < (idata$cpmE + idata$cpmMM))	
plot(c(1,4), c(min(idata$cpmM, idata$cpmE, idata$cpmMM, idata$cpmEE), max(idata$cpmM, idata$cpmE, idata$cpmMM, idata$cpmEE)), type='n', ylab="normalised count", xlab="", main=paste('(E + MM) > (M + EE) & E < MM, line number:', nrow(pat2)))
for (i in 1:nrow(pat2)) {
		lines(c(pat2$cpmM[i],pat2$cpmMM[i],pat2$cpmE[i],pat2$cpmEE[i]), col=4)
	}

pat3 <- subset(idata, idata$cpmE < cpmM)
plot(c(1,4), c(min(idata$cpmM, idata$cpmE, idata$cpmMM, idata$cpmEE), max(idata$cpmM, idata$cpmE, idata$cpmMM, idata$cpmEE)), type='n', ylab="normalised count", xlab="", main=paste('E < M, line number:', nrow(pat3)))
for (i in 1:nrow(pat3)) {
		lines(c(pat3$cpmM[i],pat3$cpmMM[i],pat3$cpmE[i],pat3$cpmEE[i]), col=6)
	}



legend("topright", inset=0.05, legend=c("virgin > mated", "virgin < mated"), pch =c(19,19), col=c(2,4), cex=0.8 )
