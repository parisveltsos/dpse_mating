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

ovpath <- "~/Library/Mobile Documents/com~apple~CloudDocs/MyDocuments/Work/Projects/gitOut/dpseMating/dedata/ovary4"
ovdata <- read.table(file.path(ovpath, "LogCPM_0.05_ovary4.txt"), header=T)
str(ovdata)

frtpath <- "~/Library/Mobile Documents/com~apple~CloudDocs/MyDocuments/Work/Projects/gitOut/dpseMating/dedata/rtract4"
frtdata <- read.table(file.path(frtpath, "LogCPM_0.05_rtract4.txt"), header=T)
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
