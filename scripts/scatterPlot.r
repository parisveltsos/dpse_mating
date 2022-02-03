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

ovpath <- "~/git/dpse_mating/output/dedata/ovary"
ovdata <- read.table(file.path(ovpath, "LogCPM_0.05_ovary.txt"), header=T)
str(ovdata)

frtpath <- "~/git/dpse_mating/output/dedata/rtract"
frtdata <- read.table(file.path(frtpath, "LogCPM_0.05_rtract.txt"), header=T)
str(frtdata)

outpath <- '~/git/dpse_mating/output/scatterplot'

pdf(file.path(outpath,"virginMatedFRTPlot.pdf"), width=6, height=6)
par(mar=c(5,5,4,3))
plot(frtdata$logFC.MvsMM*(-1), frtdata$logFC.EvsEE*(-1), xlab="logFC M vs MM", ylab="logFC E vs EE", main="", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
points(frtdata$logFC.MvsMM[frtdata$FDR.E.MvsEE.MM<0.05]*(-1), frtdata$logFC.EvsEE[frtdata$FDR.E.MvsEE.MM<0.05]*(-1), pch=2, col=cbblue)
# points(frtdata$logFC.MvsMM[frtdata$FDR.allEvsM<0.05 & frtdata$logFC.allEvsM>0]*(-1), frtdata$logFC.EvsEE[frtdata$FDR.allEvsM<0.05 & frtdata$logFC.allEvsM>0]*(-1), pch=4, col=cbred)
# points(frtdata$logFC.MvsMM[frtdata$FDR.allEvsM<0.05 & frtdata$logFC.allEvsM<0]*(-1), frtdata$logFC.EvsEE[frtdata$FDR.allEvsM<0.05 & frtdata$logFC.allEvsM<0]*(-1), pch=3, col=cborange)
points(frtdata$logFC.MvsMM[frtdata$FDR.E.MMvsM.EE<0.05]*(-1), frtdata$logFC.EvsEE[frtdata$FDR.E.MMvsM.EE<0.05]*(-1), pch=4, col=cbgreen)
abline(v=0, lty=2, lw=1, col=1)
abline(h=0, lty=2, lw=1, col=1)

legend("topleft", inset=0.05, legend=c("Virgin vs mated (244)", "Interaction (146)"), col =c(cbblue,cbgreen), pch=c(2,4), cex=.8 ) 
dev.off()


pdf(file.path(outpath,"maleEffect.pdf"), width=12, height=6)

par(mfrow=c(1,2)) 
par(mar=c(5,5,4,3))

plot(frtdata$logFC.MvsMM*(-1), frtdata$logFC.MvsME*(-1), xlab="logFC M vs MM", ylab="logFC M vs ME", main="FRT M", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
points(frtdata$logFC.MvsMM[frtdata$FDR.MvsMM<0.05 & frtdata$FDR.MvsME>0.05]*(-1), frtdata$logFC.MvsME[frtdata$FDR.MvsMM<0.05 & frtdata$FDR.MvsME>0.05]*(-1), pch=3, col=cbblue)
points(frtdata$logFC.MvsMM[frtdata$FDR.MvsME<0.05 & frtdata$FDR.MvsMM>0.05]*(-1), frtdata$logFC.MvsME[frtdata$FDR.MvsME<0.05 & frtdata$FDR.MvsMM>0.05]*(-1), pch=4, col=cbred)
points(frtdata$logFC.MvsMM[frtdata$FDR.MvsME<0.05 & frtdata$FDR.MvsMM<0.05]*(-1), frtdata$logFC.MvsME[frtdata$FDR.MvsME<0.05 & frtdata$FDR.MvsMM<0.05]*(-1), pch=20, col=cbpurple)
legend("topleft", inset=0.05, legend=c("Within-line only", "Between-line only", "Both"), pch =c(3,4,20), col=c(cbblue, cbred, cbpurple), cex=.7 )
# legend("topleft", inset=0.05, legend=c("Within-line only (3 M-biased, 55 MM-biased)", "Between-line only (103 M-biased, 155 ME-biased) ", "Both (54 virgin-biased, 51 mated-biased)"), pch =c(3,4,20), col=c(cbblue, cbred, cbpurple), cex=.7 )
abline(v=0, lty=2, lw=1, col=1)
abline(h=0, lty=2, lw=1, col=1)

plot(frtdata$logFC.EvsEE*(-1), frtdata$logFC.EvsEM*(-1), xlab="logFC E vs EE", ylab="logFC E vs EM", main="FRT E", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
points(frtdata$logFC.EvsEE[frtdata$FDR.EvsEE<0.05 & frtdata$FDR.EvsEM>0.05]*(-1), frtdata$logFC.EvsEM[frtdata$FDR.EvsEE<0.05 & frtdata$FDR.EvsEM>0.05]*(-1), pch=3, col=cbblue)
points(frtdata$logFC.EvsEE[frtdata$FDR.EvsEM<0.05 & frtdata$FDR.EvsEE>0.05]*(-1), frtdata$logFC.EvsEM[frtdata$FDR.EvsEM<0.05 & frtdata$FDR.EvsEE>0.05]*(-1), pch=4, col=cbred)
points(frtdata$logFC.EvsEE[frtdata$FDR.EvsEM<0.05 & frtdata$FDR.EvsEE<0.05]*(-1), frtdata$logFC.EvsEM[frtdata$FDR.EvsEM<0.05 & frtdata$FDR.EvsEE<0.05]*(-1), pch=20, col=cbpurple)
legend("topleft", inset=0.05, legend=c("Within line only", "Between line only", "Both"), pch =c(3,4,20), col=c(cbblue, cbred, cbpurple), cex=.7 )
# legend("topleft", inset=0.05, legend=c("Within line only (74 E-biased, 66 EE-biased)", "Between line only (25 E-biased, 5 EM-biased)", "Both (26 virgin-biased, 34 mated-biased)"), pch =c(3,4,20), col=c(cbblue, cbred, cbpurple), cex=.7 )
abline(v=0, lty=2, lw=1, col=1)
abline(h=0, lty=2, lw=1, col=1)

# plot(ovdata$logFC.MvsMM*(-1), ovdata$logFC.MvsME*(-1), xlab="logFC M vs MM", ylab="logFC M vs ME", main="OV M", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
# points(ovdata$logFC.MvsMM[ovdata$FDR.MvsMM<0.05 & ovdata$FDR.MvsME>0.05]*(-1), ovdata$logFC.MvsME[ovdata$FDR.MvsMM<0.05 & ovdata$FDR.MvsME>0.05]*(-1), pch=3, col=cbblue)
# points(ovdata$logFC.MvsMM[ovdata$FDR.MvsME<0.05 & ovdata$FDR.MvsMM>0.05]*(-1), ovdata$logFC.MvsME[ovdata$FDR.MvsME<0.05 & ovdata$FDR.MvsMM>0.05]*(-1), pch=4, col=cbred)
# points(ovdata$logFC.MvsMM[ovdata$FDR.MvsME<0.05 & ovdata$FDR.MvsMM<0.05]*(-1), ovdata$logFC.MvsME[ovdata$FDR.MvsME<0.05 & ovdata$FDR.MvsMM<0.05]*(-1), pch=20, col=cbpurple)
# abline(v=0, lty=2, lw=1, col=1)
# abline(h=0, lty=2, lw=1, col=1)
# 
# plot(ovdata$logFC.EvsEE*(-1), ovdata$logFC.EvsEM*(-1), xlab="logFC E vs EE", ylab="logFC E vs EM", main="OV E", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
# points(ovdata$logFC.EvsEE[ovdata$FDR.EvsEE<0.05 & ovdata$FDR.EvsEM>0.05]*(-1), ovdata$logFC.EvsEM[ovdata$FDR.EvsEE<0.05 & ovdata$FDR.EvsEM>0.05]*(-1), pch=3, col=cbblue)
# points(ovdata$logFC.EvsEE[ovdata$FDR.EvsEM<0.05 & ovdata$FDR.EvsEE>0.05]*(-1), ovdata$logFC.EvsEM[ovdata$FDR.EvsEM<0.05 & ovdata$FDR.EvsEE>0.05]*(-1), pch=4, col=cbred)
# points(ovdata$logFC.EvsEE[ovdata$FDR.EvsEM<0.05 & ovdata$FDR.EvsEE<0.05]*(-1), ovdata$logFC.EvsEM[ovdata$FDR.EvsEM<0.05 & ovdata$FDR.EvsEE<0.05]*(-1), pch=20, col=cbpurple)
# abline(v=0, lty=2, lw=1, col=1)
# abline(h=0, lty=2, lw=1, col=1)

dev.off()




pdf(file.path(outpath,"virginMatedOVPlot.pdf"), width=6, height=6)
par(mar=c(5,5,4,3))

plot(ovdata$logFC.allEvsM, ovdata$logFC.E.MvsEE.MM*(-1), xlab="logFC E vs M", ylab="logFC virgin vs mated", main="", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
points(ovdata$logFC.allEvsM[ovdata$FDR.allEvsM<0.05], ovdata$logFC.E.MvsEE.MM[ovdata$FDR.allEvsM<0.05]*(-1),  pch=5, col=cborange)
points(ovdata$logFC.allEvsM[ovdata$FDR.E.MvsEE.MM<0.05], ovdata$logFC.E.MvsEE.MM[ovdata$FDR.E.MvsEE.MM<0.05]*(-1),  pch=2, col=cbblue)

# points(ovdata$logFC.EvsM[ovdata$FDR.E.MMvsM.EE<0.05], ovdata$logFC.E.MvsEE.MM[ovdata$FDR.E.MMvsM.EE<0.05],  pch=4, col=cbgreen)
legend("topleft", inset=0.05, legend=c("Virgin vs mated (43)", "E vs M (604)"), col =c(cbblue, cborange), pch=c(2,5), cex=.8 ) 
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



# generation of 'venn' GO data fig 4a,c coloured points are a venn diagram, output here as txt for GO analysis 
# replace `frt` with `ov` in code below to generate the ovary data (not used)
#

outpath <- '~/git/dpse_mating/output/venn'
dir.create(file.path(outpath, "frt.datagen"))

# venn32 made manually through venny

MvsMMonly.down <- data.frame("gene"=frtdata$gene[frtdata$FDR.MvsMM < 0.05 & frtdata$FDR.MvsME > 0.05 & frtdata$logFC.MvsMM < 0])
MvsMMonly.down$FDR.MvsMMonly.down <- 0.001
MvsMMonly.up <- data.frame("gene"=frtdata$gene[frtdata$FDR.MvsMM < 0.05 & frtdata$FDR.MvsME > 0.05 & frtdata$logFC.MvsMM > 0])
MvsMMonly.up$FDR.MvsMMonly.up <- 0.001

MvsMEonly.down <- data.frame("gene"=frtdata$gene[frtdata$FDR.MvsMM > 0.05 & frtdata$FDR.MvsME < 0.05 & frtdata$logFC.MvsMM < 0])
MvsMEonly.down$FDR.MvsMEonly.down <- 0.001
MvsMEonly.up <- data.frame("gene"=frtdata$gene[frtdata$FDR.MvsMM > 0.05 & frtdata$FDR.MvsME < 0.05 & frtdata$logFC.MvsMM > 0])
MvsMEonly.up$FDR.MvsMEonly.up <- 0.001

Mboth.down <- data.frame("gene"=frtdata$gene[frtdata$FDR.MvsMM < 0.05 & frtdata$FDR.MvsME < 0.05 & frtdata$logFC.MvsMM < 0])
Mboth.down$FDR.Mboth.down <- 0.001
Mboth.up <- data.frame("gene"=frtdata$gene[frtdata$FDR.MvsMM < 0.05 & frtdata$FDR.MvsME < 0.05 & frtdata$logFC.MvsMM > 0])
Mboth.up$FDR.Mboth.up <- 0.001

frtsmall <- data.frame("gene"=frtdata$gene, "FDR.MvsME"=frtdata$FDR.MvsME, "logFC.MvsME"=frtdata$logFC.MvsME, "FDR.MvsMM"=frtdata$FDR.MvsMM, "logFC.MvsMM"=frtdata$logFC.MvsMM)

merged1 <- merge(frtsmall, MvsMMonly.down,  by.x='gene', by.y='gene', all=T)
merged2 <- merge(merged1, MvsMMonly.up,  by.x='gene', by.y='gene', all=T)
merged3 <- merge(merged2, MvsMEonly.down,  by.x='gene', by.y='gene', all=T)
merged4 <- merge(merged3, MvsMEonly.up,  by.x='gene', by.y='gene', all=T)
merged5 <- merge(merged4, Mboth.up,  by.x='gene', by.y='gene', all=T)
merged6 <- merge(merged5, Mboth.down,  by.x='gene', by.y='gene', all=T)

merged6$FDR.MvsMMonly.down[is.na(merged6$FDR.MvsMMonly.down)] <- 1
merged6$FDR.MvsMMonly.up[is.na(merged6$FDR.MvsMMonly.up)] <- 1
merged6$FDR.MvsMEonly.down[is.na(merged6$FDR.MvsMEonly.down)] <- 1
merged6$FDR.MvsMEonly.up[is.na(merged6$FDR.MvsMEonly.up)] <- 1
merged6$FDR.Mboth.up[is.na(merged6$FDR.Mboth.up)] <- 1
merged6$FDR.Mboth.down[is.na(merged6$FDR.Mboth.down)] <- 1

dir.create(file.path(outpath, "frt.MvsMMonly.up"))
dir.create(file.path(outpath, "frt.MvsMMonly.down"))
dir.create(file.path(outpath, "frt.MvsMEonly.up"))
dir.create(file.path(outpath, "frt.MvsMEonly.down"))
dir.create(file.path(outpath, "frt.Mboth.up"))
dir.create(file.path(outpath, "frt.Mboth.down"))


write.table(data.frame(merged6$gene, merged6$FDR.MvsMMonly.down), file=file.path(outpath, "frt.MvsMMonly.down/GO_pvalues.txt"), quote=F, row.names=F, sep="\t")
write.table(data.frame(merged6$gene, merged6$FDR.MvsMMonly.up), file=file.path(outpath, "frt.MvsMMonly.up/GO_pvalues.txt"), quote=F, row.names=F, sep="\t")
write.table(data.frame(merged6$gene, merged6$FDR.MvsMEonly.down), file=file.path(outpath, "frt.MvsMEonly.down/GO_pvalues.txt"), quote=F, row.names=F, sep="\t")
write.table(data.frame(merged6$gene, merged6$FDR.MvsMEonly.up), file=file.path(outpath, "frt.MvsMEonly.up/GO_pvalues.txt"), quote=F, row.names=F, sep="\t")
write.table(data.frame(merged6$gene, merged6$FDR.Mboth.down), file=file.path(outpath, "frt.Mboth.down/GO_pvalues.txt"), quote=F, row.names=F, sep="\t")
write.table(data.frame(merged6$gene, merged6$FDR.Mboth.up), file=file.path(outpath, "frt.Mboth.up/GO_pvalues.txt"), quote=F, row.names=F, sep="\t")







# Show interaction (not used)
# par(mar=c(5,5,4,3))
# plot(frtdata$logFC.allEvsM, frtdata$logFC.E.MvsEE.MM, xlab="logFC E vs M", ylab="logFC virgin vs mated", main="", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
# points(frtdata$logFC.allEvsM[frtdata$FDR.allEvsM<0.05], frtdata$logFC.E.MvsEE.MM[frtdata$FDR.allEvsM<0.05],  pch=4, col=cborange)
# points(frtdata$logFC.allEvsM[frtdata$FDR.E.MvsEE.MM<0.05], frtdata$logFC.E.MvsEE.MM[frtdata$FDR.E.MvsEE.MM<0.05],  pch=3, col=cbblue)
# points(frtdata$logFC.allEvsM[frtdata$FDR.E.MMvsM.EE<0.05], frtdata$logFC.E.MvsEE.MM[frtdata$FDR.E.MMvsM.EE<0.05],  pch=1, col=cbgreen)
# abline(v=0, lty=2, lw=1, col=1)
# abline(h=0, lty=2, lw=1, col=1)
# legend("bottomright", inset=0.05, legend=c("Virgin vs mated (362)", "E vs M (205)", "Interaction (169)"), col =c(cbblue,cborange,cbgreen), pch=c(3,4,1), cex=.8 ) 

# Plot male vs female effect (not used)
# pdf(file.path(outpath,"maleFemalePlot.pdf"), width=12, height=6)
# par(mfrow=c(1,2)) 
# par(mar=c(5,5,4,3))
# 
# plot(frtdata$logFC.EE.MEvsEM.MM, frtdata$logFC.EE.EMvsMM.ME, xlab="logFC EE.MEvsEM.MM (male effect)", ylab="logFC EE.EMvsMM.ME (female effect)", main="", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
# points(frtdata$logFC.EE.MEvsEM.MM[frtdata$FDR.EE.MEvsEM.MM<0.05], frtdata$logFC.EE.EMvsMM.ME[frtdata$FDR.EE.MEvsEM.MM<0.05], pch=20, col=cbblue)
# points(frtdata$logFC.EE.MEvsEM.MM[frtdata$FDR.EE.EMvsMM.ME<0.05], frtdata$logFC.EE.EMvsMM.ME[frtdata$FDR.EE.EMvsMM.ME<0.05], pch=20, col=cborange)
# points(frtdata$logFC.EE.MEvsEM.MM[frtdata$FDR.EE.MMvsEM.ME<0.05], frtdata$logFC.EE.EMvsMM.ME[frtdata$FDR.EE.MMvsEM.ME<0.05], pch=20, col=cbgreen)
# points(frtdata$logFC.EE.MEvsEM.MM[frtdata$FDR.MMvsEM<0.05], frtdata$logFC.EE.EMvsMM.ME[frtdata$FDR.MMvsEM<0.05], pch=4, col=cbred)
# points(frtdata$logFC.EE.MEvsEM.MM[frtdata$FDR.EEvsME<0.05], frtdata$logFC.EE.EMvsMM.ME[frtdata$FDR.EEvsME<0.05], pch=3, col=cbgreen)
# abline(v=0, lty=2, lw=1, col=1)
# abline(h=0, lty=2, lw=1, col=1)
# legend("topleft", inset=0.05, legend=c("Male effect", "Female effect", "Interaction"), col =c(cbblue,cborange,cbgreen), pch=c(20), cex=.8 ) 
# legend("bottomright", inset=0.05, legend=c("Male effect", "Female effect", "Interaction", "MM.EM", "EE.ME"), col =c(cbblue,cborange,cbgreen, cbred, cbgreen), pch=c(20,20,20,4,3), cex=.8 ) 
# 
# par(mar=c(5,5,4,3))
# plot(ovdata$logFC.EE.MEvsEM.MM, ovdata$logFC.EE.EMvsMM.ME, xlab="logFC EE.MEvsEM.MM (male effect)", ylab="logFC EE.EMvsMM.ME (female effect)", main="", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
# points(ovdata$logFC.EE.MEvsEM.MM[ovdata$FDR.EE.MEvsEM.MM<0.05], ovdata$logFC.EE.EMvsMM.ME[ovdata$FDR.EE.MEvsEM.MM<0.05], pch=20, col=cbblue)
# points(ovdata$logFC.EE.MEvsEM.MM[ovdata$FDR.EE.EMvsMM.ME<0.05], ovdata$logFC.EE.EMvsMM.ME[ovdata$FDR.EE.EMvsMM.ME<0.05], pch=20, col=cborange)
# points(ovdata$logFC.EE.MEvsEM.MM[ovdata$FDR.EE.MMvsEM.ME<0.05], ovdata$logFC.EE.EMvsMM.ME[ovdata$FDR.EE.MMvsEM.ME<0.05], pch=20, col=cbgreen)
# points(ovdata$logFC.EE.MEvsEM.MM[ovdata$FDR.MMvsEM<0.05], ovdata$logFC.EE.EMvsMM.ME[ovdata$FDR.MMvsEM<0.05], pch=4, col=cbred)
# points(ovdata$logFC.EE.MEvsEM.MM[ovdata$FDR.EEvsME<0.05], ovdata$logFC.EE.EMvsMM.ME[ovdata$FDR.EEvsME<0.05], pch=3, col=cbgreen)
# abline(v=0, lty=2, lw=1, col=1)
# abline(h=0, lty=2, lw=1, col=1)
# dev.off()

# Alternative male vs female effect plot - not used
# pdf(file.path(outpath,"maleFemalePlot2.pdf"), width=12, height=6)
# par(mfrow=c(1,2)) 
# par(mar=c(5,5,4,3))
# 
# plot(frtdata$logFC.EEvsEM, frtdata$logFC.MMvsME, xlab="logFC EEvsEM (E female effect)", ylab="logFC MMvsME (M female effect)", main="", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
# points(frtdata$logFC.EEvsEM[frtdata$FDR.EE.MEvsEM.MM<0.05], frtdata$logFC.MMvsME[frtdata$FDR.EE.MEvsEM.MM<0.05], pch=20, col=cbblue)
# points(frtdata$logFC.EEvsEM[frtdata$FDR.EE.EMvsMM.ME<0.05], frtdata$logFC.MMvsME[frtdata$FDR.EE.EMvsMM.ME<0.05], pch=20, col=cborange)
# points(frtdata$logFC.EEvsEM[frtdata$FDR.EE.MMvsEM.ME<0.05], frtdata$logFC.MMvsME[frtdata$FDR.EE.MMvsEM.ME<0.05], pch=20, col=cbgreen)
# abline(v=0, lty=2, lw=1, col=1)
# abline(h=0, lty=2, lw=1, col=1)
# legend("topleft", inset=0.05, legend=c("Male effect", "Female effect", "Interaction"), col =c(cbblue,cborange,cbgreen), pch=c(20), cex=.8 ) 
# 
# par(mar=c(5,5,4,3))
# plot(ovdata$logFC.EE.MEvsEM.MM, ovdata$logFC.EE.EMvsMM.ME, xlab="logFC EEvsEM (E female effect)", ylab="logFC MMvsME (M female effect)", main="", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
# points(ovdata$logFC.EEvsEM[ovdata$FDR.EE.MEvsEM.MM<0.05], ovdata$logFC.MMvsME[ovdata$FDR.EE.MEvsEM.MM<0.05], pch=20, col=cbblue)
# points(ovdata$logFC.EEvsEM[ovdata$FDR.EE.EMvsMM.ME<0.05], ovdata$logFC.MMvsME[ovdata$FDR.EE.EMvsMM.ME<0.05], pch=20, col=cborange)
# points(ovdata$logFC.EEvsEM[ovdata$FDR.EE.MMvsEM.ME<0.05], ovdata$logFC.MMvsME[ovdata$FDR.EE.MMvsEM.ME<0.05], pch=20, col=cbgreen)
# abline(v=0, lty=2, lw=1, col=1)
# abline(h=0, lty=2, lw=1, col=1)
# dev.off()


# Alternative male vs female plot - not used
# pdf(file.path(outpath,"maleFemalePlot3.pdf"), width=12, height=6)
# par(mfrow=c(1,2)) 
# par(mar=c(5,5,4,3))
# 
# plot(frtdata$logFC.EEvsME, frtdata$logFC.MMvsEM, xlab="logFC EEvsME (E male effect)", ylab="logFC MMvsEM (M male effect)", main="", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
# points(frtdata$logFC.EEvsME[frtdata$FDR.EE.MEvsEM.MM<0.05], frtdata$logFC.MMvsEM[frtdata$FDR.EE.MEvsEM.MM<0.05], pch=20, col=cbblue)
# points(frtdata$logFC.EEvsME[frtdata$FDR.EE.EMvsMM.ME<0.05], frtdata$logFC.MMvsEM[frtdata$FDR.EE.EMvsMM.ME<0.05], pch=20, col=cborange)
# points(frtdata$logFC.EEvsME[frtdata$FDR.EE.MMvsEM.ME<0.05], frtdata$logFC.MMvsEM[frtdata$FDR.EE.MMvsEM.ME<0.05], pch=20, col=cbgreen)
# abline(v=0, lty=2, lw=1, col=1)
# abline(h=0, lty=2, lw=1, col=1)
# legend("bottomleft", inset=0.05, legend=c("Male effect", "Female effect", "Interaction"), col =c(cbblue,cborange,cbgreen), pch=c(20), cex=.8 ) 
# 
# par(mar=c(5,5,4,3))
# plot(ovdata$logFC.EEvsME, ovdata$logFC.MMvsEM, xlab="logFC EEvsME (E male effect)", ylab="logFC MMvsEM (M male effect)", main="", cex.main=1.8, cex.lab=1.3, col=sgrey, pch=20)
# points(ovdata$logFC.EEvsME[ovdata$FDR.EE.MEvsEM.MM<0.05], ovdata$logFC.MMvsEM[ovdata$FDR.EE.MEvsEM.MM<0.05], pch=20, col=cbblue)
# points(ovdata$logFC.EEvsME[ovdata$FDR.EE.EMvsMM.ME<0.05], ovdata$logFC.MMvsEM[ovdata$FDR.EE.EMvsMM.ME<0.05], pch=20, col=cborange)
# points(ovdata$logFC.EEvsME[ovdata$FDR.EE.MMvsEM.ME<0.05], ovdata$logFC.MMvsEM[ovdata$FDR.EE.MMvsEM.ME<0.05], pch=20, col=cbgreen)
# abline(v=0, lty=2, lw=1, col=1)
# abline(h=0, lty=2, lw=1, col=1)
# dev.off()


# Density plots v vs mated - not used
# par(mfrow=c(2,2)) 
# 
# plot( density(frtdata$logFC.EvsEE), col=cbred, main="E females")
# lines( density(frtdata$logFC.EvsEE), col=cbred, main="E females", lw=2)
# lines(density(frtdata$logFC.EvsEM), col=cbpink, lw=2)
# legend("topright", inset=0.05, legend=c("EvsEE", "EvsEM"), pch =c(19,19), col=c(cbred,cbpink) )
# abline(v=0, lty=2, lw=1, col=1)
# 
# plot(density(frtdata$logFC.MvsMM), col=cbblue, main="M females")
# lines(density(frtdata$logFC.MvsMM), col=cbblue, main="M females", lw=2)
# lines(density(frtdata$logFC.MvsME), col=cbpurple, lw=2)
# legend("topright", inset=0.05, legend=c("MvsMM", "MvsME"), pch =c(19,19), col=c(cbblue,cbpurple) )
# abline(v=0, lty=2, lw=1, col=1)
# 
# plot( density(frtdata$logFC.EvsEE[frtdata$FDR.EvsEE<0.05]), col=cbred, main="E females")
# lines( density(frtdata$logFC.EvsEE[frtdata$FDR.EvsEE<0.05]), col=cbred, main="E females", lw=2)
# lines(density(frtdata$logFC.EvsEM[frtdata$FDR.EvsEM<0.05]), col=cbpink, lw=2)
# legend("topright", inset=0.05, legend=c("EvsEE", "EvsEM"), pch =c(19,19), col=c(cbred,cbpink) )
# abline(v=0, lty=2, lw=1, col=1)
# 
# plot(density(frtdata$logFC.MvsME[frtdata$FDR.MvsME<0.05]), col=cbpurple, main="M females")
# lines(density(frtdata$logFC.MvsME[frtdata$FDR.MvsME<0.05]), col=cbpurple, main="M females", lw=2)
# lines(density(frtdata$logFC.MvsMM[frtdata$FDR.MvsMM<0.05]), col=cbblue, main="M females", lw=2)
# legend("topright", inset=0.05, legend=c("MvsMM", "MvsME"), pch =c(19,19), col=c(cbblue,cbpurple) )
# abline(v=0, lty=2, lw=1, col=1)
