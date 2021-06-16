datapath <- "~/Library/Mobile Documents/com~apple~CloudDocs/MyDocuments/Work/Projects/gitOut/dpseMating/dedata/ovary"
kdata <- read.table(file.path(datapath, "LogCPM_0.05_ovary.txt"), header=T)
str(kdata)

sgrey = rgb(0, 0, 0, 0.1)
cbblue <- "#2271B2"
cblblue <- '#3DB7E9'
cbpink <- "#F748A5"
cbgreen <- "#359B73"
cbbrown <- "#D55E00"
cborange <- "#E69F00"
cbyellow <- "#F0E442"
cbpurple <- "#8400CD" 
cbred <- "#E20134" # not with brown, orange, yellow


# par(mfrow=c(1,2)) 
par(mar=c(5,5,4,3))
plot(kdata$logFC.MvsMM*(-1), kdata$logFC.EvsEE*(-1), xlab="logFC M vs MM", ylab="logFC E vs EE", main="Line x Mating", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
points(kdata$logFC.MvsMM[kdata$FDR.E.MMvsM.EE<0.05]*(-1), kdata$logFC.EvsEE[kdata$FDR.E.MMvsM.EE<0.05]*(-1), pch=20, col=cbgreen)
points(kdata$logFC.MvsMM[kdata$FDR.E.MvsEE.MM<0.05]*(-1), kdata$logFC.EvsEE[kdata$FDR.E.MvsEE.MM<0.05]*(-1), pch=20, col=cbblue)
abline(v=0, lty=2, lw=1, col=1)
abline(h=0, lty=2, lw=1, col=1)

legend("topleft", inset=0.05, legend=c("DE All virgin vs mated main: congruent", "DE Interaction: non-congruent"), col =c(cbblue,cbgreen), pch=c(20,20), cex=.8 ) 

# plot(kdata$logFC.MvsMM, kdata$logFC.EvsEE, xlab="M vs MM", ylab="E vs EE", main="nothing", cex.main=1.8, cex.lab=1.3, col='gray80')
# points(kdata$logFC.MvsMM[kdata$FDR.EvsEE<0.05], kdata$logFC.EvsEE[kdata$FDR.EvsEE<0.05], pch=3, col="blue")
# points(kdata$logFC.MvsMM[kdata$FDR.MvsMM<0.05], kdata$logFC.EvsEE[kdata$FDR.MvsMM<0.05], pch=4, col="orange")
# abline(h=0, col='darkgreen')
# abline(v=0, col='darkgreen')


par(mfrow=c(1,2)) 
par(mar=c(5,5,4,3))
plot(kdata$logFC.MMvsME, kdata$logFC.EEvsEM, xlab="logFC MM vs ME", ylab="logFC EE vs EM", main="Within line effect of male", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
points(kdata$logFC.MMvsME[kdata$FDR.EE.MMvsEM.ME<0.05], kdata$logFC.EEvsEM[kdata$FDR.EE.MMvsEM.ME<0.05], pch=20, col=cbgreen)
points(kdata$logFC.MMvsME[kdata$FDR.EEvsEM<0.05], kdata$logFC.EEvsEM[kdata$FDR.EEvsEM<0.05], pch=20, col=cbblue)
points(kdata$logFC.MMvsME[kdata$FDR.MMvsME<0.05], kdata$logFC.EEvsEM[kdata$FDR.MMvsME<0.05], pch=20, col=cborange)
legend("topleft", inset=0.05, legend=c("EE.MMvsEM.ME", "EEvsEM", "MMvsME"), pch =c(20,20,20), col=c(cbgreen,cbblue,cborange), cex=.8 )
abline(v=0, lty=2, lw=1, col=1)
abline(h=0, lty=2, lw=1, col=1)

plot(kdata$logFC.MMvsME, kdata$logFC.EEvsEM, xlab="logFC MM vs ME", ylab="logFC EE vs EM", main="Within line effect of male", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
points(kdata$logFC.MMvsME[kdata$FDR.EE.MEvsEM.MM<0.05], kdata$logFC.EEvsEM[kdata$FDR.EE.MEvsEM.MM<0.05], pch=20, col=cblblue)
points(kdata$logFC.MMvsME[kdata$FDR.EEvsME<0.05], kdata$logFC.EEvsEM[kdata$FDR.EEvsME<0.05], pch=20, col=cbpink)
points(kdata$logFC.MMvsME[kdata$FDR.MMvsEM<0.05], kdata$logFC.EEvsEM[kdata$FDR.MMvsEM<0.05], pch=20, col=cbpurple)
legend("topleft", inset=0.05, legend=c("EE.MEvsEM.MM", "EEvsME", "MMvsEM"), pch =c(20,20,20), col=c(cblblue,cbpink,cbpurple), cex=.8 )
abline(v=0, lty=2, lw=1, col=1)
abline(h=0, lty=2, lw=1, col=1)


par(mfrow=c(1,2)) 
par(mar=c(5,5,4,3))
plot(kdata$logFC.EE.MMvsEM.ME, kdata$logFC.E.MvsEE.MM, xlab="logFC EE.MM vs EM.ME", ylab="logFC E.M vs EE.MM", main="Within line effect of male incl Vir/mated", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
points(kdata$logFC.EE.MMvsEM.ME[kdata$FDR.EE.MMvsEM.ME<0.05], kdata$logFC.E.MvsEE.MM[kdata$FDR.EE.MMvsEM.ME<0.05], pch=20, col=cbgreen)
points(kdata$logFC.EE.MMvsEM.ME[kdata$FDR.EEvsEM<0.05], kdata$logFC.E.MvsEE.MM[kdata$FDR.EEvsEM<0.05], pch=20, col=cbblue)
points(kdata$logFC.EE.MMvsEM.ME[kdata$FDR.MMvsME<0.05], kdata$logFC.E.MvsEE.MM[kdata$FDR.MMvsME<0.05], pch=20, col=cborange)
abline(v=0, lty=2, lw=1, col=1)
abline(h=0, lty=2, lw=1, col=1)
legend("topright", inset=0.05, legend=c("EE.MMvsEM.ME", "EEvsEM", "MMvsEM"), pch =c(20,20,20), col=c(cbgreen,cbblue,cborange), cex=.8 )

plot(kdata$logFC.EE.MMvsEM.ME, kdata$logFC.E.MvsEE.MM, xlab="logFC EE.MM vs EM.ME", ylab="logFC E.M vs EE.MM", main="Within line effect of male incl Vir/mated", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
points(kdata$logFC.EE.MMvsEM.ME[kdata$FDR.MMvsEM<0.05], kdata$logFC.E.MvsEE.MM[kdata$FDR.MMvsEM<0.05], pch=20, col=cbpurple)
points(kdata$logFC.EE.MMvsEM.ME[kdata$FDR.EEvsME<0.05], kdata$logFC.E.MvsEE.MM[kdata$FDR.EEvsME<0.05], pch=20, col=cbpink)
points(kdata$logFC.EE.MMvsEM.ME[kdata$FDR.EE.MEvsEM.MM<0.05], kdata$logFC.E.MvsEE.MM[kdata$FDR.EE.MEvsEM.MM<0.05], pch=20, col=cblblue)
abline(v=0, lty=2, lw=1, col=1)
abline(h=0, lty=2, lw=1, col=1)
legend("topright", inset=0.05, legend=c("EE.MEvsEM.MM", "EEvsME", "MMvsEM"), pch =c(20,20,20), col=c(cblblue,cbpink,cbpurple), cex=.8 )










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
points(kdata$logFC.E.MvsEE.MM[kdata$FDR.MvsEE<0.05], kdata$logFC.E.MMvsM.EE[kdata$FDR.MvsEE<0.05], pch=2, col=7)
points(kdata$logFC.E.MvsEE.MM[kdata$FDR.EvsMM<0.05], kdata$logFC.E.MMvsM.EE[kdata$FDR.EvsMM<0.05], pch=3, col=8)
legend("topleft", inset=0.05, legend=c("E.MMvsM.EE", "E.MvsEE.MM", "all E vs M"), pch =c(20,4,3), col=c(4,2, 6) ) 
abline(h=0, col='darkgreen')
abline(v=0, col='darkgreen')

plot(kdata$logFC.E.MvsEE.MM, kdata$logFC.allEvsM, xlab="virgin vs mated logFC", ylab="E vs M all", main="virgin/mated vs E.EE/M.MM", cex.main=1.8, cex.lab=1.3, col='gray80')
points(kdata$logFC.E.MvsEE.MM[kdata$FDR.E.MMvsM.EE<0.05], kdata$logFC.allEvsM[kdata$FDR.E.MMvsM.EE<0.05], pch=20, col=4)
points(kdata$logFC.E.MvsEE.MM[kdata$FDR.E.MvsEE.MM<0.05], kdata$logFC.allEvsM[kdata$FDR.E.MvsEE.MM<0.05], pch=4, col=2)
points(kdata$logFC.E.MvsEE.MM[kdata$FDR.allEvsM<0.05], kdata$logFC.allEvsM[kdata$FDR.allEvsM<0.05], pch=3, col=6)
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
