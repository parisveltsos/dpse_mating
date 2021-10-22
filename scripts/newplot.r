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

ovpath <- "~/Library/Mobile Documents/com~apple~CloudDocs/MyDocuments/Work/Projects/gitOut/dpseMating/dedata/ovary"
ovdata <- read.table(file.path(ovpath, "LogCPM_0.05_ovary.txt"), header=T)
str(ovdata)

frtpath <- "~/Library/Mobile Documents/com~apple~CloudDocs/MyDocuments/Work/Projects/gitOut/dpseMating/dedata/rtract"
frtdata <- read.table(file.path(frtpath, "LogCPM_0.05_rtract.txt"), header=T)
str(frtdata)

outpath <- '/Users/pveltsos/Library/Mobile Documents/com~apple~CloudDocs/MyDocuments/Work/Projects/gitOut/dpseMating/scatterplot'

pdf(file.path(outpath,"virginMatedFRTPlot.pdf"), width=6, height=6)
par(mar=c(5,5,4,3))
plot(frtdata$logFC.MvsMM*(-1), frtdata$logFC.EvsEE*(-1), xlab="logFC M vs MM", ylab="logFC E vs EE", main="", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
points(frtdata$logFC.MvsMM[frtdata$FDR.E.MvsEE.MM<0.05]*(-1), frtdata$logFC.EvsEE[frtdata$FDR.E.MvsEE.MM<0.05]*(-1), pch=2, col=cbblue)
# points(frtdata$logFC.MvsMM[frtdata$FDR.allEvsM<0.05 & frtdata$logFC.allEvsM>0]*(-1), frtdata$logFC.EvsEE[frtdata$FDR.allEvsM<0.05 & frtdata$logFC.allEvsM>0]*(-1), pch=4, col=cbred)
# points(frtdata$logFC.MvsMM[frtdata$FDR.allEvsM<0.05 & frtdata$logFC.allEvsM<0]*(-1), frtdata$logFC.EvsEE[frtdata$FDR.allEvsM<0.05 & frtdata$logFC.allEvsM<0]*(-1), pch=3, col=cborange)
points(frtdata$logFC.MvsMM[frtdata$FDR.E.MMvsM.EE<0.05]*(-1), frtdata$logFC.EvsEE[frtdata$FDR.E.MMvsM.EE<0.05]*(-1), pch=4, col=cbgreen)
abline(v=0, lty=2, lw=1, col=1)
abline(h=0, lty=2, lw=1, col=1)

legend("topleft", inset=0.05, legend=c("Virgin vs mated (362)", "Interaction (169)"), col =c(cbblue,cbgreen), pch=c(2,4), cex=.8 ) 
dev.off()


pdf(file.path(outpath,"virginMatedOVPlot.pdf"), width=6, height=6)
par(mar=c(5,5,4,3))

plot(ovdata$logFC.allEvsM, ovdata$logFC.E.MvsEE.MM, xlab="logFC E vs M", ylab="logFC virgin vs mated", main="", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
points(ovdata$logFC.allEvsM[ovdata$FDR.allEvsM<0.05], ovdata$logFC.E.MvsEE.MM[ovdata$FDR.allEvsM<0.05],  pch=2, col=cbblue)
points(ovdata$logFC.allEvsM[ovdata$FDR.E.MvsEE.MM<0.05], ovdata$logFC.E.MvsEE.MM[ovdata$FDR.E.MvsEE.MM<0.05],  pch=5, col=cborange)
# points(ovdata$logFC.EvsM[ovdata$FDR.E.MMvsM.EE<0.05], ovdata$logFC.E.MvsEE.MM[ovdata$FDR.E.MMvsM.EE<0.05],  pch=4, col=cbgreen)
legend("topleft", inset=0.05, legend=c("Virgin vs mated (44)", "E vs M (570)"), col =c(cbblue, cborange), pch=c(2,5), cex=.8 ) 
# plot(ovdata$logFC.MvsMM*(-1), ovdata$logFC.EvsEE*(-1), xlab="logFC M vs MM", ylab="logFC E vs EE", main="", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
# points(ovdata$logFC.MvsMM[ovdata$FDR.E.MvsEE.MM<0.05]*(-1), ovdata$logFC.EvsEE[ovdata$FDR.E.MvsEE.MM<0.05]*(-1), pch=20, col=cbblue)
# points(ovdata$logFC.MvsMM[ovdata$FDR.EvsM<0.05]*(-1), ovdata$logFC.EvsEE[ovdata$FDR.EvsM<0.05]*(-1), pch=4, col=cbred)
# points(ovdata$logFC.MvsMM[ovdata$FDR.E.MMvsM.EE<0.05]*(-1), ovdata$logFC.EvsEE[ovdata$FDR.E.MMvsM.EE<0.05]*(-1), pch=20, col=cbgreen)
# points(ovdata$logFC.MvsMM[ovdata$FDR.EvsM<0.05 & ovdata$logFC.EvsM>0]*(-1), ovdata$logFC.EvsEE[ovdata$FDR.EvsM<0.05 & ovdata$logFC.EvsM>0]*(-1), pch=4, col=cbred)
# points(ovdata$logFC.MvsMM[ovdata$FDR.EvsM<0.05 & ovdata$logFC.EvsM<0]*(-1), ovdata$logFC.EvsEE[ovdata$FDR.EvsM<0.05 & ovdata$logFC.EvsM<0]*(-1), pch=3, col=cborange)
abline(v=0, lty=2, lw=1, col=1)
abline(h=0, lty=2, lw=1, col=1)

# legend("topleft", inset=0.05, legend=c("Virgin vs mated"), col =c(cbblue), pch=c(20), cex=.8 ) 
dev.off()

par(mar=c(5,5,4,3))
plot(frtdata$logFC.allEvsM, frtdata$logFC.E.MvsEE.MM, xlab="logFC E vs M", ylab="logFC virgin vs mated", main="", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
points(frtdata$logFC.allEvsM[frtdata$FDR.allEvsM<0.05], frtdata$logFC.E.MvsEE.MM[frtdata$FDR.allEvsM<0.05],  pch=4, col=cborange)
points(frtdata$logFC.allEvsM[frtdata$FDR.E.MvsEE.MM<0.05], frtdata$logFC.E.MvsEE.MM[frtdata$FDR.E.MvsEE.MM<0.05],  pch=3, col=cbblue)
points(frtdata$logFC.allEvsM[frtdata$FDR.E.MMvsM.EE<0.05], frtdata$logFC.E.MvsEE.MM[frtdata$FDR.E.MMvsM.EE<0.05],  pch=1, col=cbgreen)
abline(v=0, lty=2, lw=1, col=1)
abline(h=0, lty=2, lw=1, col=1)
legend("bottomright", inset=0.05, legend=c("Virgin vs mated (362)", "E vs M (205)", "Interaction (169)"), col =c(cbblue,cborange,cbgreen), pch=c(3,4,1), cex=.8 ) 


pdf(file.path(outpath,"maleFemalePlot.pdf"), width=12, height=6)
par(mfrow=c(1,2)) 
par(mar=c(5,5,4,3))

plot(frtdata$logFC.EE.MEvsEM.MM, frtdata$logFC.EE.EMvsMM.ME, xlab="logFC EE.MEvsEM.MM (male effect)", ylab="logFC EE.EMvsMM.ME (female effect)", main="", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
points(frtdata$logFC.EE.MEvsEM.MM[frtdata$FDR.EE.MEvsEM.MM<0.05], frtdata$logFC.EE.EMvsMM.ME[frtdata$FDR.EE.MEvsEM.MM<0.05], pch=20, col=cbblue)
points(frtdata$logFC.EE.MEvsEM.MM[frtdata$FDR.EE.EMvsMM.ME<0.05], frtdata$logFC.EE.EMvsMM.ME[frtdata$FDR.EE.EMvsMM.ME<0.05], pch=20, col=cborange)
points(frtdata$logFC.EE.MEvsEM.MM[frtdata$FDR.EE.MMvsEM.ME<0.05], frtdata$logFC.EE.EMvsMM.ME[frtdata$FDR.EE.MMvsEM.ME<0.05], pch=20, col=cbgreen)
points(frtdata$logFC.EE.MEvsEM.MM[frtdata$FDR.MMvsEM<0.05], frtdata$logFC.EE.EMvsMM.ME[frtdata$FDR.MMvsEM<0.05], pch=4, col=cbred)
points(frtdata$logFC.EE.MEvsEM.MM[frtdata$FDR.EEvsME<0.05], frtdata$logFC.EE.EMvsMM.ME[frtdata$FDR.EEvsME<0.05], pch=3, col=cbgreen)
abline(v=0, lty=2, lw=1, col=1)
abline(h=0, lty=2, lw=1, col=1)
legend("topleft", inset=0.05, legend=c("Male effect", "Female effect", "Interaction"), col =c(cbblue,cborange,cbgreen), pch=c(20), cex=.8 ) 
legend("bottomright", inset=0.05, legend=c("Male effect", "Female effect", "Interaction", "MM.EM", "EE.ME"), col =c(cbblue,cborange,cbgreen, cbred, cbgreen), pch=c(20,20,20,4,3), cex=.8 ) 

par(mar=c(5,5,4,3))
plot(ovdata$logFC.EE.MEvsEM.MM, ovdata$logFC.EE.EMvsMM.ME, xlab="logFC EE.MEvsEM.MM (male effect)", ylab="logFC EE.EMvsMM.ME (female effect)", main="", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
points(ovdata$logFC.EE.MEvsEM.MM[ovdata$FDR.EE.MEvsEM.MM<0.05], ovdata$logFC.EE.EMvsMM.ME[ovdata$FDR.EE.MEvsEM.MM<0.05], pch=20, col=cbblue)
points(ovdata$logFC.EE.MEvsEM.MM[ovdata$FDR.EE.EMvsMM.ME<0.05], ovdata$logFC.EE.EMvsMM.ME[ovdata$FDR.EE.EMvsMM.ME<0.05], pch=20, col=cborange)
points(ovdata$logFC.EE.MEvsEM.MM[ovdata$FDR.EE.MMvsEM.ME<0.05], ovdata$logFC.EE.EMvsMM.ME[ovdata$FDR.EE.MMvsEM.ME<0.05], pch=20, col=cbgreen)
points(ovdata$logFC.EE.MEvsEM.MM[ovdata$FDR.MMvsEM<0.05], ovdata$logFC.EE.EMvsMM.ME[ovdata$FDR.MMvsEM<0.05], pch=4, col=cbred)
points(ovdata$logFC.EE.MEvsEM.MM[ovdata$FDR.EEvsME<0.05], ovdata$logFC.EE.EMvsMM.ME[ovdata$FDR.EEvsME<0.05], pch=3, col=cbgreen)
abline(v=0, lty=2, lw=1, col=1)
abline(h=0, lty=2, lw=1, col=1)
dev.off()




plot(ovdata$logFC.MMvsEM, ovdata$logFC.EEvsME, xlab="logFC MMvsEM", ylab="logFC EEvsME", main="", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
points(ovdata$logFC.MMvsEM[ovdata$FDR.MMvsEM<0.05], ovdata$logFC.EEvsME[ovdata$FDR.MMvsEM<0.05], pch=20, col=cbblue)
points(ovdata$logFC.MMvsEM[ovdata$FDR.EEvsME<0.05], ovdata$logFC.EEvsME[ovdata$FDR.EEvsME<0.05], pch=20, col=cbgreen)
points(ovdata$logFC.MMvsEM[ovdata$FDR.EEvsEM<0.05], ovdata$logFC.EEvsME[ovdata$FDR.EEvsEM<0.05], pch=4, col=cbred)
points(ovdata$logFC.MMvsEM[ovdata$FDR.EE.EMvsEM.MM<0.05], ovdata$logFC.EEvsME[ovdata$FDR.EE.EMvsEM.MM<0.05], pch=3, col=cborange) # do this
legend("bottomleft", inset=0.05, legend=c("MMvsEM", "EEvsME", "EEvsEM"), col =c(cbblue,cbgreen,cbred), pch=c(20,20,4), cex=.8 ) 



pdf(file.path(outpath,"maleFemalePlot2.pdf"), width=12, height=6)
par(mfrow=c(1,2)) 
par(mar=c(5,5,4,3))

plot(frtdata$logFC.EEvsEM, frtdata$logFC.MMvsME, xlab="logFC EEvsEM (E female effect)", ylab="logFC MMvsME (M female effect)", main="", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
points(frtdata$logFC.EEvsEM[frtdata$FDR.EE.MEvsEM.MM<0.05], frtdata$logFC.MMvsME[frtdata$FDR.EE.MEvsEM.MM<0.05], pch=20, col=cbblue)
points(frtdata$logFC.EEvsEM[frtdata$FDR.EE.EMvsMM.ME<0.05], frtdata$logFC.MMvsME[frtdata$FDR.EE.EMvsMM.ME<0.05], pch=20, col=cborange)
points(frtdata$logFC.EEvsEM[frtdata$FDR.EE.MMvsEM.ME<0.05], frtdata$logFC.MMvsME[frtdata$FDR.EE.MMvsEM.ME<0.05], pch=20, col=cbgreen)
abline(v=0, lty=2, lw=1, col=1)
abline(h=0, lty=2, lw=1, col=1)
legend("topleft", inset=0.05, legend=c("Male effect", "Female effect", "Interaction"), col =c(cbblue,cborange,cbgreen), pch=c(20), cex=.8 ) 

par(mar=c(5,5,4,3))
plot(ovdata$logFC.EE.MEvsEM.MM, ovdata$logFC.EE.EMvsMM.ME, xlab="logFC EEvsEM (E female effect)", ylab="logFC MMvsME (M female effect)", main="", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
points(ovdata$logFC.EEvsEM[ovdata$FDR.EE.MEvsEM.MM<0.05], ovdata$logFC.MMvsME[ovdata$FDR.EE.MEvsEM.MM<0.05], pch=20, col=cbblue)
points(ovdata$logFC.EEvsEM[ovdata$FDR.EE.EMvsMM.ME<0.05], ovdata$logFC.MMvsME[ovdata$FDR.EE.EMvsMM.ME<0.05], pch=20, col=cborange)
points(ovdata$logFC.EEvsEM[ovdata$FDR.EE.MMvsEM.ME<0.05], ovdata$logFC.MMvsME[ovdata$FDR.EE.MMvsEM.ME<0.05], pch=20, col=cbgreen)
abline(v=0, lty=2, lw=1, col=1)
abline(h=0, lty=2, lw=1, col=1)
dev.off()



pdf(file.path(outpath,"maleFemalePlot3.pdf"), width=12, height=6)
par(mfrow=c(1,2)) 
par(mar=c(5,5,4,3))

plot(frtdata$logFC.EEvsME, frtdata$logFC.MMvsEM, xlab="logFC EEvsME (E male effect)", ylab="logFC MMvsEM (M male effect)", main="", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
points(frtdata$logFC.EEvsME[frtdata$FDR.EE.MEvsEM.MM<0.05], frtdata$logFC.MMvsEM[frtdata$FDR.EE.MEvsEM.MM<0.05], pch=20, col=cbblue)
points(frtdata$logFC.EEvsME[frtdata$FDR.EE.EMvsMM.ME<0.05], frtdata$logFC.MMvsEM[frtdata$FDR.EE.EMvsMM.ME<0.05], pch=20, col=cborange)
points(frtdata$logFC.EEvsME[frtdata$FDR.EE.MMvsEM.ME<0.05], frtdata$logFC.MMvsEM[frtdata$FDR.EE.MMvsEM.ME<0.05], pch=20, col=cbgreen)
abline(v=0, lty=2, lw=1, col=1)
abline(h=0, lty=2, lw=1, col=1)
legend("bottomleft", inset=0.05, legend=c("Male effect", "Female effect", "Interaction"), col =c(cbblue,cborange,cbgreen), pch=c(20), cex=.8 ) 

par(mar=c(5,5,4,3))
plot(ovdata$logFC.EEvsME, ovdata$logFC.MMvsEM, xlab="logFC EEvsME (E male effect)", ylab="logFC MMvsEM (M male effect)", main="", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
points(ovdata$logFC.EEvsME[ovdata$FDR.EE.MEvsEM.MM<0.05], ovdata$logFC.MMvsEM[ovdata$FDR.EE.MEvsEM.MM<0.05], pch=20, col=cbblue)
points(ovdata$logFC.EEvsME[ovdata$FDR.EE.EMvsMM.ME<0.05], ovdata$logFC.MMvsEM[ovdata$FDR.EE.EMvsMM.ME<0.05], pch=20, col=cborange)
points(ovdata$logFC.EEvsME[ovdata$FDR.EE.MMvsEM.ME<0.05], ovdata$logFC.MMvsEM[ovdata$FDR.EE.MMvsEM.ME<0.05], pch=20, col=cbgreen)
abline(v=0, lty=2, lw=1, col=1)
abline(h=0, lty=2, lw=1, col=1)
dev.off()

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







pdf(file.path(outpath,"4libs_frt.pdf"), width=9, height=9)

par(mfrow=c(2,2)) 

par(mar=c(5,5,4,3))
plot(frtdata$logFC.EEvsMM, frtdata$logFC.EEvsEM, xlab="logFC EE vs MM", ylab="logFC EE vs EM", main="", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
points(frtdata$logFC.EEvsMM[frtdata$FDR.EEvsMM<0.05], frtdata$logFC.EEvsEM[frtdata$FDR.EEvsMM<0.05],  pch=4, col=cborange)
points(frtdata$logFC.EEvsMM[frtdata$FDR.EEvsEM<0.05], frtdata$logFC.EEvsEM[frtdata$FDR.EEvsEM<0.05],  pch=3, col=cbblue)
points(frtdata$logFC.EEvsMM[frtdata$FDR.MMvsEM<0.05], frtdata$logFC.EEvsEM[frtdata$FDR.MMvsEM<0.05],  pch=2, col=cbred)
points(frtdata$logFC.EEvsMM[frtdata$FDR.MMvsME<0.05], frtdata$logFC.EEvsEM[frtdata$FDR.MMvsME<0.05],  pch=6, col=cbpurple)
abline(v=0, lty=2, lw=1, col=1)
abline(h=0, lty=2, lw=1, col=1)
legend("topleft", inset=0.05, legend=c("EE vs MM (113)", "EE vs EM (0)"), col =c(cborange,cbblue), pch=c(4,3), cex=.8 ) 

par(mar=c(5,5,4,3))
plot(frtdata$logFC.EEvsMM, frtdata$logFC.MMvsME, xlab="logFC EE vs MM", ylab="logFC MM vs ME", main="", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
points(frtdata$logFC.EEvsMM[frtdata$FDR.EEvsMM<0.05], frtdata$logFC.MMvsME[frtdata$FDR.EEvsMM<0.05],  pch=4, col=cborange)
points(frtdata$logFC.EEvsMM[frtdata$FDR.MMvsME<0.05], frtdata$logFC.MMvsME[frtdata$FDR.MMvsME<0.05],  pch=3, col=cbblue)
points(frtdata$logFC.EEvsMM[frtdata$FDR.MMvsEM<0.05], frtdata$logFC.MMvsME[frtdata$FDR.MMvsEM<0.05],  pch=2, col=cbred)
abline(v=0, lty=2, lw=1, col=1)
abline(h=0, lty=2, lw=1, col=1)
legend("topleft", inset=0.05, legend=c("EE vs MM (113)", "MM vs ME (0)"), col =c(cborange,cbblue), pch=c(4,3), cex=.8 ) 

par(mar=c(5,5,4,3))
plot(frtdata$logFC.EEvsMM, frtdata$logFC.EEvsME, xlab="logFC EE vs MM", ylab="logFC EE vs ME", main="", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
points(frtdata$logFC.EEvsMM[frtdata$FDR.EEvsMM<0.05], frtdata$logFC.EEvsME[frtdata$FDR.EEvsMM<0.05],  pch=4, col=cborange)
points(frtdata$logFC.EEvsMM[frtdata$FDR.EEvsME<0.05], frtdata$logFC.EEvsME[frtdata$FDR.EEvsME<0.05],  pch=3, col=cbblue)
abline(v=0, lty=2, lw=1, col=1)
abline(h=0, lty=2, lw=1, col=1)
legend("bottomright", inset=0.05, legend=c("EE vs MM (113)", "EE vs ME (43)"), col =c(cborange,cbblue), pch=c(4,3), cex=.8 ) 


par(mar=c(5,5,4,3))
plot(frtdata$logFC.EEvsMM, frtdata$logFC.MMvsEM, xlab="logFC EE vs MM", ylab="logFC MM vs EM", main="", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
points(frtdata$logFC.EEvsMM[frtdata$FDR.EEvsMM<0.05], frtdata$logFC.MMvsEM[frtdata$FDR.EEvsMM<0.05],  pch=4, col=cborange)
points(frtdata$logFC.EEvsMM[frtdata$FDR.MMvsEM<0.05], frtdata$logFC.MMvsEM[frtdata$FDR.MMvsEM<0.05],  pch=3, col=cbblue)
abline(v=0, lty=2, lw=1, col=1)
abline(h=0, lty=2, lw=1, col=1)
legend("topright", inset=0.05, legend=c("EE vs MM (113)", "MM vs EM (289)"), col =c(cborange,cbblue), pch=c(4,3), cex=.8 ) 

dev.off()



pdf(file.path(outpath,"4libs_ov.pdf"), width=9, height=9)

par(mfrow=c(2,2)) 

par(mar=c(5,5,4,3))
plot(ovdata$logFC.EEvsMM, ovdata$logFC.EEvsEM, xlab="logFC EE vs MM", ylab="logFC EE vs EM", main="", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
points(ovdata$logFC.EEvsMM[ovdata$FDR.EEvsMM<0.05], ovdata$logFC.EEvsEM[ovdata$FDR.EEvsMM<0.05],  pch=4, col=cborange)
points(ovdata$logFC.EEvsMM[ovdata$FDR.EEvsEM<0.05], ovdata$logFC.EEvsEM[ovdata$FDR.EEvsEM<0.05],  pch=3, col=cbblue)
abline(v=0, lty=2, lw=1, col=1)
abline(h=0, lty=2, lw=1, col=1)
legend("bottomleft", inset=0.05, legend=c("EE vs MM (231)", "EE vs EM (27)"), col =c(cborange,cbblue), pch=c(4,3), cex=.8 ) 

par(mar=c(5,5,4,3))
plot(ovdata$logFC.EEvsMM, ovdata$logFC.MMvsME, xlab="logFC EE vs MM", ylab="logFC MM vs ME", main="", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
points(ovdata$logFC.EEvsMM[ovdata$FDR.EEvsMM<0.05], ovdata$logFC.MMvsME[ovdata$FDR.EEvsMM<0.05],  pch=4, col=cborange)
points(ovdata$logFC.EEvsMM[ovdata$FDR.MMvsME<0.05], ovdata$logFC.MMvsME[ovdata$FDR.MMvsME<0.05],  pch=3, col=cbblue)
abline(v=0, lty=2, lw=1, col=1)
abline(h=0, lty=2, lw=1, col=1)
legend("topleft", inset=0.05, legend=c("EE vs MM (231)", "MM vs ME (0)"), col =c(cborange,cbblue), pch=c(4,3), cex=.8 ) 

par(mar=c(5,5,4,3))
plot(ovdata$logFC.EEvsMM, ovdata$logFC.EEvsME, xlab="logFC EE vs MM", ylab="logFC EE vs ME", main="", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
points(ovdata$logFC.EEvsMM[ovdata$FDR.EEvsMM<0.05], ovdata$logFC.EEvsME[ovdata$FDR.EEvsMM<0.05],  pch=4, col=cborange)
points(ovdata$logFC.EEvsMM[ovdata$FDR.EEvsME<0.05], ovdata$logFC.EEvsME[ovdata$FDR.EEvsME<0.05],  pch=3, col=cbblue)
abline(v=0, lty=2, lw=1, col=1)
abline(h=0, lty=2, lw=1, col=1)
legend("topleft", inset=0.05, legend=c("EE vs MM (231)", "EE vs ME (391)"), col =c(cborange,cbblue), pch=c(4,3), cex=.8 ) 


par(mar=c(5,5,4,3))
plot(ovdata$logFC.EEvsMM, ovdata$logFC.MMvsEM, xlab="logFC EE vs MM", ylab="logFC MM vs EM", main="", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
points(ovdata$logFC.EEvsMM[ovdata$FDR.EEvsMM<0.05], ovdata$logFC.MMvsEM[ovdata$FDR.EEvsMM<0.05],  pch=4, col=cborange)
points(ovdata$logFC.EEvsMM[ovdata$FDR.MMvsEM<0.05], ovdata$logFC.MMvsEM[ovdata$FDR.MMvsEM<0.05],  pch=3, col=cbblue)
abline(v=0, lty=2, lw=1, col=1)
abline(h=0, lty=2, lw=1, col=1)
legend("topright", inset=0.05, legend=c("EE vs MM (231)", "MM vs EM (144)"), col =c(cborange,cbblue), pch=c(4,3), cex=.8 ) 

dev.off()




pdf(file.path(outpath,"6libs_frt.pdf"), width=9, height=9)

par(mfrow=c(2,2)) 

par(mar=c(5,5,4,3))
plot(frtdata$logFC.EvsEE, frtdata$logFC.EvsEM, xlab="logFC E vs EE", ylab="logFC E vs EM", main="", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
points(frtdata$logFC.EvsEE[frtdata$FDR.EvsEE<0.05], frtdata$logFC.EvsEM[frtdata$FDR.EvsEE<0.05],  pch=4, col=cborange)
points(frtdata$logFC.EvsEE[frtdata$FDR.EvsEM<0.05], frtdata$logFC.EvsEM[frtdata$FDR.EvsEM<0.05],  pch=3, col=cbblue)
abline(v=0, lty=2, lw=1, col=1)
abline(h=0, lty=2, lw=1, col=1)
legend("topleft", inset=0.05, legend=c("E vs EE (200)", "E vs EM (90)"), col =c(cborange,cbblue), pch=c(4,3), cex=.8 ) 


par(mar=c(5,5,4,3))
plot(frtdata$logFC.MvsMM, frtdata$logFC.MvsME, xlab="logFC M vs MM", ylab="logFC M vs ME", main="", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
points(frtdata$logFC.MvsMM[frtdata$FDR.MvsME<0.05], frtdata$logFC.MvsME[frtdata$FDR.MvsME<0.05],  pch=3, col=cbblue)
points(frtdata$logFC.MvsMM[frtdata$FDR.MvsMM<0.05], frtdata$logFC.MvsME[frtdata$FDR.MvsMM<0.05],  pch=4, col=cborange)
abline(v=0, lty=2, lw=1, col=1)
abline(h=0, lty=2, lw=1, col=1)
legend("bottomright", inset=0.05, legend=c("M vs MM (163)", "M vs ME (363)"), col =c(cborange,cbblue), pch=c(4,3), cex=.8 ) 


par(mar=c(5,5,4,3))
plot(frtdata$logFC.EvsEE, frtdata$logFC.EvsME, xlab="logFC E vs EE", ylab="logFC E vs ME", main="", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
points(frtdata$logFC.EvsEE[frtdata$FDR.EvsEE<0.05], frtdata$logFC.EvsME[frtdata$FDR.EvsEE<0.05],  pch=4, col=cborange)
points(frtdata$logFC.EvsEE[frtdata$FDR.EvsME<0.05], frtdata$logFC.EvsME[frtdata$FDR.EvsME<0.05],  pch=3, col=cbblue)
abline(v=0, lty=2, lw=1, col=1)
abline(h=0, lty=2, lw=1, col=1)
legend("topleft", inset=0.05, legend=c("E vs EE (200)", "E vs ME (152)"), col =c(cborange,cbblue), pch=c(4,3), cex=.8 ) 

par(mar=c(5,5,4,3))
plot(frtdata$logFC.MvsMM, frtdata$logFC.MvsEM, xlab="logFC M vs MM", ylab="logFC M vs EM", main="", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
points(frtdata$logFC.MvsMM[frtdata$FDR.MvsEM<0.05], frtdata$logFC.MvsEM[frtdata$FDR.MvsEM<0.05],  pch=3, col=cbblue)
points(frtdata$logFC.MvsMM[frtdata$FDR.MvsMM<0.05], frtdata$logFC.MvsEM[frtdata$FDR.MvsMM<0.05],  pch=4, col=cborange)
abline(v=0, lty=2, lw=1, col=1)
abline(h=0, lty=2, lw=1, col=1)
legend("bottomright", inset=0.05, legend=c("M vs MM (163)", "M vs ME (2387)"), col =c(cborange,cbblue), pch=c(4,3), cex=.8 ) 

dev.off()






pdf(file.path(outpath,"6libs_ov.pdf"), width=9, height=9)

par(mfrow=c(2,2)) 

par(mar=c(5,5,4,3))
plot(ovdata$logFC.EvsEE, ovdata$logFC.EvsEM, xlab="logFC E vs EE", ylab="logFC E vs EM", main="", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
points(ovdata$logFC.EvsEE[ovdata$FDR.EvsEE<0.05], ovdata$logFC.EvsEM[ovdata$FDR.EvsEE<0.05],  pch=4, col=cborange)
points(ovdata$logFC.EvsEE[ovdata$FDR.EvsEM<0.05], ovdata$logFC.EvsEM[ovdata$FDR.EvsEM<0.05],  pch=3, col=cbblue)
abline(v=0, lty=2, lw=1, col=1)
abline(h=0, lty=2, lw=1, col=1)
legend("topleft", inset=0.05, legend=c("E vs EE (6)", "E vs EM (0)"), col =c(cborange,cbblue), pch=c(4,3), cex=.8 ) 


par(mar=c(5,5,4,3))
plot(ovdata$logFC.MvsMM, ovdata$logFC.MvsME, xlab="logFC M vs MM", ylab="logFC M vs ME", main="", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
points(ovdata$logFC.MvsMM[ovdata$FDR.MvsME<0.05], ovdata$logFC.MvsME[ovdata$FDR.MvsME<0.05],  pch=3, col=cbblue)
points(ovdata$logFC.MvsMM[ovdata$FDR.MvsMM<0.05], ovdata$logFC.EvsME[ovdata$FDR.MvsMM<0.05],  pch=4, col=cborange)
abline(v=0, lty=2, lw=1, col=1)
abline(h=0, lty=2, lw=1, col=1)
legend("bottomright", inset=0.05, legend=c("M vs MM (1)", "M vs ME (3)"), col =c(cborange,cbblue), pch=c(4,3), cex=.8 ) 


par(mar=c(5,5,4,3))
plot(ovdata$logFC.EvsEE, ovdata$logFC.EvsME, xlab="logFC E vs EE", ylab="logFC E vs ME", main="", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
points(ovdata$logFC.EvsEE[ovdata$FDR.EvsEE<0.05], ovdata$logFC.EvsME[ovdata$FDR.EvsEE<0.05],  pch=4, col=cborange)
points(ovdata$logFC.EvsEE[ovdata$FDR.EvsME<0.05], ovdata$logFC.EvsME[ovdata$FDR.EvsME<0.05],  pch=3, col=cbblue)
abline(v=0, lty=2, lw=1, col=1)
abline(h=0, lty=2, lw=1, col=1)
legend("topleft", inset=0.05, legend=c("E vs EE (6)", "E vs ME (270)"), col =c(cborange,cbblue), pch=c(4,3), cex=.8 ) 

par(mar=c(5,5,4,3))
plot(ovdata$logFC.MvsMM, ovdata$logFC.MvsEM, xlab="logFC M vs MM", ylab="logFC M vs EM", main="", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
points(ovdata$logFC.MvsMM[ovdata$FDR.MvsEM<0.05], ovdata$logFC.MvsEM[ovdata$FDR.MvsEM<0.05],  pch=3, col=cbblue)
points(ovdata$logFC.MvsMM[ovdata$FDR.MvsMM<0.05], ovdata$logFC.MvsEM[ovdata$FDR.MvsMM<0.05],  pch=4, col=cborange)
abline(v=0, lty=2, lw=1, col=1)
abline(h=0, lty=2, lw=1, col=1)
legend("bottomright", inset=0.05, legend=c("M vs MM (1)", "M vs ME (174)"), col =c(cborange,cbblue), pch=c(4,3), cex=.8 ) 

dev.off()








pdf(file.path(outpath,"sexEffects.pdf"), width=9, height=9)

par(mfrow=c(2,2)) 

par(mar=c(5,5,4,3))
plot(frtdata$logFC.MMvsEM, frtdata$logFC.EEvsME, xlab="logFC MM vs EM", ylab="logFC EE vs ME", main="FRT", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
points(frtdata$logFC.MMvsEM[frtdata$FDR.MMvsEM<0.05], frtdata$logFC.EEvsME[frtdata$FDR.MMvsEM<0.05],  pch=4, col=cborange)
points(frtdata$logFC.MMvsEM[frtdata$FDR.EEvsME<0.05], frtdata$logFC.EEvsME[frtdata$FDR.EEvsME<0.05],  pch=3, col=cbblue)
points(frtdata$logFC.MMvsEM[frtdata$FDR.EEvsMM<0.05 & frtdata$logFC.EEvsMM>0], frtdata$logFC.EEvsME[frtdata$FDR.EEvsMM<0.05 & frtdata$logFC.EEvsMM>0],  pch=2, col=cbred)
points(frtdata$logFC.MMvsEM[frtdata$FDR.EEvsMM<0.05 & frtdata$logFC.EEvsMM<0], frtdata$logFC.EEvsME[frtdata$FDR.EEvsMM<0.05 & frtdata$logFC.EEvsMM<0],  pch=6, col=cbpurple)
abline(v=0, lty=2, lw=1, col=1)
abline(h=0, lty=2, lw=1, col=1)
legend("topright", inset=0.05, legend=c("MM vs EM (289)", "EE vs ME (43)", "EE up MM down (19)", "EE down MM up (94)"), col =c(cborange,cbblue,cbred,cbpurple), pch=c(4,3,2,6), cex=.7 ) 


par(mar=c(5,5,4,3))
plot(frtdata$logFC.MMvsME, frtdata$logFC.EEvsEM, xlab="logFC MM vs ME", ylab="logFC EE vs EM", main="FRT", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
points(frtdata$logFC.MMvsME[frtdata$FDR.MMvsME<0.05], frtdata$logFC.EEvsEM[frtdata$FDR.MMvsME<0.05],  pch=4, col=cborange)
points(frtdata$logFC.MMvsME[frtdata$FDR.EEvsEM<0.05], frtdata$logFC.EEvsME[frtdata$FDR.EEvsEM<0.05],  pch=3, col=cbblue)
points(frtdata$logFC.MMvsME[frtdata$FDR.EEvsMM<0.05 & frtdata$logFC.EEvsMM>0], frtdata$logFC.EEvsEM[frtdata$FDR.EEvsMM<0.05 & frtdata$logFC.EEvsMM>0],  pch=2, col=cbred)
points(frtdata$logFC.MMvsME[frtdata$FDR.EEvsMM<0.05 & frtdata$logFC.EEvsMM<0], frtdata$logFC.EEvsEM[frtdata$FDR.EEvsMM<0.05 & frtdata$logFC.EEvsMM<0],  pch=6, col=cbpurple)
abline(v=0, lty=2, lw=1, col=1)
abline(h=0, lty=2, lw=1, col=1)
legend("topleft", inset=0.05, legend=c("MM vs EM (0)", "EE vs ME (0)", "EE up MM down (19)", "EE down MM up (94)"), col =c(cborange,cbblue,cbred,cbpurple), pch=c(4,3,2,6), cex=.7 ) 



plot(ovdata$logFC.MMvsEM, ovdata$logFC.EEvsME, xlab="logFC MM vs EM", ylab="logFC EE vs ME", main="ov", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
points(ovdata$logFC.MMvsEM[ovdata$FDR.MMvsEM<0.05], ovdata$logFC.EEvsME[ovdata$FDR.MMvsEM<0.05],  pch=4, col=cborange)
points(ovdata$logFC.MMvsEM[ovdata$FDR.EEvsME<0.05], ovdata$logFC.EEvsME[ovdata$FDR.EEvsME<0.05],  pch=3, col=cbblue)
points(ovdata$logFC.MMvsEM[ovdata$FDR.EEvsMM<0.05 & ovdata$logFC.EEvsMM>0], ovdata$logFC.EEvsME[ovdata$FDR.EEvsMM<0.05 & ovdata$logFC.EEvsMM>0],  pch=2, col=cbred)
points(ovdata$logFC.MMvsEM[ovdata$FDR.EEvsMM<0.05 & ovdata$logFC.EEvsMM<0], ovdata$logFC.EEvsME[ovdata$FDR.EEvsMM<0.05 & ovdata$logFC.EEvsMM<0],  pch=6, col=cbpurple)
abline(v=0, lty=2, lw=1, col=1)
abline(h=0, lty=2, lw=1, col=1)
legend("topright", inset=0.05, legend=c("MM vs EM (114)", "EE vs ME (391)", "EE up MM down (108)", "EE down MM up (123)"), col =c(cborange,cbblue,cbred,cbpurple), pch=c(4,3,2,6), cex=.7 ) 


par(mar=c(5,5,4,3))
plot(ovdata$logFC.MMvsME, ovdata$logFC.EEvsEM, xlab="logFC MM vs ME", ylab="logFC EE vs EM", main="ov", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
points(ovdata$logFC.MMvsME[ovdata$FDR.MMvsME<0.05], ovdata$logFC.EEvsEM[ovdata$FDR.MMvsME<0.05],  pch=4, col=cborange)
points(ovdata$logFC.MMvsME[ovdata$FDR.EEvsEM<0.05], ovdata$logFC.EEvsME[ovdata$FDR.EEvsEM<0.05],  pch=3, col=cbblue)
points(ovdata$logFC.MMvsME[ovdata$FDR.EEvsMM<0.05 & ovdata$logFC.EEvsMM>0], ovdata$logFC.EEvsEM[ovdata$FDR.EEvsMM<0.05 & ovdata$logFC.EEvsMM>0],  pch=2, col=cbred)
points(ovdata$logFC.MMvsME[ovdata$FDR.EEvsMM<0.05 & ovdata$logFC.EEvsMM<0], ovdata$logFC.EEvsEM[ovdata$FDR.EEvsMM<0.05 & ovdata$logFC.EEvsMM<0],  pch=6, col=cbpurple)
abline(v=0, lty=2, lw=1, col=1)
abline(h=0, lty=2, lw=1, col=1)
legend("bottomright", inset=0.05, legend=c("MM vs EM (0)", "EE vs ME (27)", "EE up MM down (108)", "EE down MM up (123)"), col =c(cborange,cbblue,cbred,cbpurple), pch=c(4,3,2,6), cex=.7 ) 

dev.off()





par(mar=c(5,5,4,3))
plot(frtdata$logFC.EEvsMM, frtdata$logFC.EMvsME, xlab="logFC EE vs MM", ylab="logFC EM vs ME", main="frt", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
points(frtdata$logFC.EEvsMM[frtdata$FDR.EEvsMM<0.05], frtdata$logFC.EMvsME[frtdata$FDR.EEvsMM<0.05],  pch=3, col=cbblue)
points(frtdata$logFC.EEvsMM[frtdata$FDR.MMvsEM<0.05], frtdata$logFC.EMvsME[frtdata$FDR.MMvsEM<0.05],  pch=2, col=cbred)
points(frtdata$logFC.EEvsMM[frtdata$FDR.MMvsME<0.05], frtdata$logFC.EMvsME[frtdata$FDR.MMvsME<0.05],  pch=6, col=cbpurple)
points(frtdata$logFC.EEvsMM[frtdata$FDR.EEvsEM<0.05], frtdata$logFC.EMvsME[frtdata$FDR.EEvsEM<0.05],  pch=4, col=cbred)
points(frtdata$logFC.EEvsMM[frtdata$FDR.EEvsME<0.05], frtdata$logFC.EMvsME[frtdata$FDR.EEvsME<0.05],  pch=4, col=cbpurple)
abline(v=0, lty=2, lw=1, col=1)
abline(h=0, lty=2, lw=1, col=1)
legend("topleft", inset=0.05, legend=c("EE vs MM (113)", "MM vs EM (288)", "MM vs ME (0)", "EE vs EM (0)", "EE vs ME (43)"), col =c(cbblue, cbred, cbpurple, cbred, cbpurple), pch=c(3,2,6,4,4), cex=.8 ) 




par(mar=c(5,5,4,3))
plot(ovdata$logFC.EEvsMM, ovdata$logFC.EMvsME, xlab="logFC EE vs MM", ylab="logFC EM vs ME", main="ovaries", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
points(ovdata$logFC.EEvsMM[ovdata$FDR.EEvsMM<0.05], ovdata$logFC.EMvsME[ovdata$FDR.EEvsMM<0.05],  pch=3, col=cbblue)
points(ovdata$logFC.EEvsMM[ovdata$FDR.MMvsEM<0.05], ovdata$logFC.EMvsME[ovdata$FDR.MMvsEM<0.05],  pch=2, col=cbred)
points(ovdata$logFC.EEvsMM[ovdata$FDR.MMvsME<0.05], ovdata$logFC.EMvsME[ovdata$FDR.MMvsME<0.05],  pch=6, col=cbpurple)
points(ovdata$logFC.EEvsMM[ovdata$FDR.EEvsEM<0.05], ovdata$logFC.EMvsME[ovdata$FDR.EEvsEM<0.05],  pch=4, col=cbred)
points(ovdata$logFC.EEvsMM[ovdata$FDR.EEvsME<0.05], ovdata$logFC.EMvsME[ovdata$FDR.EEvsME<0.05],  pch=4, col=cbpurple)
abline(v=0, lty=2, lw=1, col=1)
abline(h=0, lty=2, lw=1, col=1)
legend("topleft", inset=0.05, legend=c("EE vs MM (238)", "MM vs EM (144)", "MM vs ME (0)", "EE vs EM (27)", "EE vs ME (391)"), col =c(cbblue, cbred, cbpurple, cbred, cbpurple), pch=c(3,2,6,4,4), cex=.8 ) 



Mave <- frtdata$M_rt_1 + frtdata$M_rt_2 + frtdata$M_rt_3 + frtdata$M_rt_4
Eave <- frtdata$E_rt_1 + frtdata$E_rt_2 + frtdata$E_rt_3 + frtdata$E_rt_4
MEave <- frtdata$ME_rt_1 + frtdata$ME_rt_2 + frtdata$ME_rt_3 + frtdata$ME_rt_4
EMave <- frtdata$EM_rt_1 + frtdata$EM_rt_2 + frtdata$EM_rt_3 + frtdata$EM_rt_4
MMave <- frtdata$MM_rt_1 + frtdata$MM_rt_2 + frtdata$MM_rt_3 + frtdata$MM_rt_4
EEave <- frtdata$EE_rt_1 + frtdata$EE_rt_2 + frtdata$EE_rt_3 + frtdata$EE_rt_4

head(data.frame(frtdata[order(frtdata$FDR.MMvsEM, decreasing=F),]))

sdata <- subset(frtdata, frtdata$gene==“FBgn0078687”)

listGenes <- c("FBgn0078687", "FBgn0248682")

listE <- list()
listM <- list()
listEE <- list()
listMM <- list()

for (i in listGenes) {
listEE[[i]] <- EEave[frtdata$gene==i]
listMM[[i]] <- MMave[frtdata$gene==i]
listEM[[i]] <- EMave[frtdata$gene==i]
listME[[i]] <- MEave[frtdata$gene==i]

# listse[[i]] <- as.matrix(cbind(	c(Mave[frtdata$gene==i], MMave[frtdata$gene==i]), 
	# 							c(Eave[frtdata$gene==i], EEave[frtdata$gene==i])))
}


gene <- "FBgn0248682" # EE vs ME
gene <- "FBgn0250047" # MM vs EM
gene <- "FBgn0250659" # only EE vs ME

frtdata$FDR.EEvsME[frtdata$gene==gene]
frtdata$FDR.MMvsEM[frtdata$gene==gene]

tail(subset(frtdata, frtdata$FDR.EEvsME<0.05 & frtdata$FDR.MMvsEM>0.05))

neo <- as.matrix(cbind(	c(MMave[frtdata$gene==gene], EEave[frtdata$gene==gene]), 
					c(MEave[frtdata$gene==gene], EMave[frtdata$gene==gene])))
								
colnames(neo) <- c("Within line", "Between line")
row.names(neo) <- c("M", "E")
neo


par(mfrow=c(2,2)) 

plot( density(frtdata$logFC.EvsEE), col=cbred, main="E females")
lines( density(frtdata$logFC.EvsEE), col=cbred, main="E females", lw=2)
lines(density(frtdata$logFC.EvsEM), col=cbpink, lw=2)
legend("topright", inset=0.05, legend=c("EvsEE", "EvsEM"), pch =c(19,19), col=c(cbred,cbpink) )
abline(v=0, lty=2, lw=1, col=1)

plot(density(frtdata$logFC.MvsMM), col=cbblue, main="M females")
lines(density(frtdata$logFC.MvsMM), col=cbblue, main="M females", lw=2)
lines(density(frtdata$logFC.MvsME), col=cbpurple, lw=2)
legend("topright", inset=0.05, legend=c("MvsMM", "MvsME"), pch =c(19,19), col=c(cbblue,cbpurple) )
abline(v=0, lty=2, lw=1, col=1)

plot( density(frtdata$logFC.EvsEE[frtdata$FDR.EvsEE<0.05]), col=cbred, main="E females")
lines( density(frtdata$logFC.EvsEE[frtdata$FDR.EvsEE<0.05]), col=cbred, main="E females", lw=2)
lines(density(frtdata$logFC.EvsEM[frtdata$FDR.EvsEM<0.05]), col=cbpink, lw=2)
legend("topright", inset=0.05, legend=c("EvsEE", "EvsEM"), pch =c(19,19), col=c(cbred,cbpink) )
abline(v=0, lty=2, lw=1, col=1)

plot(density(frtdata$logFC.MvsME[frtdata$FDR.MvsME<0.05]), col=cbpurple, main="M females")
lines(density(frtdata$logFC.MvsME[frtdata$FDR.MvsME<0.05]), col=cbpurple, main="M females", lw=2)
lines(density(frtdata$logFC.MvsMM[frtdata$FDR.MvsMM<0.05]), col=cbblue, main="M females", lw=2)
legend("topright", inset=0.05, legend=c("MvsMM", "MvsME"), pch =c(19,19), col=c(cbblue,cbpurple) )
abline(v=0, lty=2, lw=1, col=1)

plot( density(frtdata$logFC.MvsMM), xlim=range( c(frtdata$logFC.MvsMM, frtdata$logFC.MvsME) ) )
lines(density(frtdata$logFC.MvsME), col=2)
