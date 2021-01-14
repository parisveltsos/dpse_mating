# merge all data from supplementary tables of Karr 2019 https://doi.org/10.1074/mcp.RA118.001098 
datapath <- '~/git/dpse_mating/input/agcp/'
outpath <- '~/git/dpse_mating/output/agcp/'
dir.create(file.path(outpath))

alldata <- read.table(file.path(datapath, 'all.txt'), header=T) # Table S1
str(alldata) 
agcpsdata <- read.table(file.path(datapath, 'agcps.txt'), header=T) # Table S4
str(agcpsdata) 
orthodata <- read.table(file.path(datapath, 'ortho.txt'), header=T) # Table S3
str(orthodata)

mergedA <- merge(alldata, agcpsdata, by.x="id", by.y="id", all=T )
agcpsdata <- merge(mergedA, orthodata, by.x="id", by.y="id", all=T )

# findlaydata <- read.table('~/git/dpse_mating/input/agcp/not used/findlay.txt', header=T) # Findlay 2014 https://doi.org/10.1371/journal.pgen.1004108
# str(findlaydata)

agdata <- read.table('~/git/dpse_mating/output/dedata/agland/LogCPM_0.05_agland.txt', header=T) # data from DE analysis
str(agdata)

# merged2 <- merge(agdata, findlaydata, by.x="gene", by.y="id", all=F ) # could be used to compare to Findlay 2014 genes

merged1 <- merge(agdata, agcpsdata, by.x="gene", by.y="id", all=F )
summary(merged1)
merged1$de <- merged1$FDR.maleVirgins_E.M < 0.05  
merged1$defdr1 <- merged1$FDR.maleVirgins_E.M < 0.05 & merged1$logFC.maleVirgins_E.M > 1
merged1$defdr2 <- merged1$FDR.maleVirgins_E.M < 0.05 & merged1$logFC.maleVirgins_E.M > 2

deonly <- subset(merged1, merged1$FDR.maleVirgins_E.M < 0.05 & merged1$all=="all") 
summary(deonly)


pdf(file.path(outpath,paste('agcpPlot.pdf', sep="")), width=12, height=6) # Fig 2
par(mfrow=c(1,2)) 
plot(density(merged1$logFC.maleVirgins_E.M), xlim=range(-2,2), main='AG E vs M (all)', col=2)
lines(density(merged1$logFC.maleVirgins_E.M), col=2, lw=2)
lines(density(agdata$logFC.maleVirgins_E.M), col=1, lw=2)
lines(density(merged1$logFC.maleVirgins_E.M[merged1$agcpS=="agcpS"], na.rm=T), col=3, lw=2)
lines(density(merged1$logFC.maleVirgins_E.M[merged1$ortho=="ortho"], na.rm=T), col=4, lw=2)
rug(agdata$logFC.maleVirgins_E.M, col=1, ticksize=0.01, line=2.5)
rug(merged1$logFC.maleVirgins_E.M, col=2, ticksize=0.01, line=3.0)
# legend('topleft', inset=0.05, legend=c('All genes', 'Findlay homologues (76/89/138)'), pch =c(19,19), col=c(1,2), cex=0.8 )
# legend('topright', inset=0.05, legend=c('W=48552, p=0.19'),  cex=0.8 )
legend('topleft', inset=0.05, legend=c('All RNAseq', 'All proteome', 'Secretome', 'SFP'), pch =c(19,19), col=c(1,2,3,4), cex=0.8 )
# legend('topright', inset=0.05, legend=c('W=21345000, p=<0.001'),  cex=0.8 )
abline(v=0, lty=2)
wilcox.test(agdata$logFC.maleVirgins_E.M, merged1$logFC.maleVirgins_E.M) # significant
wilcox.test(merged1$logFC.maleVirgins_E.M, merged1$logFC.maleVirgins_E.M[merged1$agcpS=="agcpS"]) # NS
wilcox.test(merged1$logFC.maleVirgins_E.M, merged1$logFC.maleVirgins_E.M[merged1$ortho=="ortho"]) # NS
wilcox.test(agdata$logFC.maleVirgins_E.M, merged1$logFC.maleVirgins_E.M[merged1$agcpS=="agcpS"]) # significant
wilcox.test(agdata$logFC.maleVirgins_E.M, merged1$logFC.maleVirgins_E.M[merged1$ortho=="ortho"]) # NS
wilcox.test(merged1$logFC.maleVirgins_E.M[merged1$agcpS=="agcpS"], merged1$logFC.maleVirgins_E.M[merged1$ortho=="ortho"]) # NS

plot(density(merged1$logFC.maleVirgins_E.M[merged1$FDR.maleVirgins_E.M<0.05]), main='AG E vs M (FDR<0.05)', xlim=range(-6,6), col=1)
lines(density(merged1$logFC.maleVirgins_E.M[merged1$FDR.maleVirgins_E.M<0.05]), col=1, lw=2)
lines(density(agdata$logFC.maleVirgins_E.M[agdata$FDR.maleVirgins_E.M<0.05]), col=2, lw=2)
# lines(density(merged1$logFC.maleVirgins_E.M[merged1$FDR.maleVirgins_E.M<0.05 & merged1$agcpS=="agcpS"], na.rm=T), col=3)
# lines(density(merged1$logFC.maleVirgins_E.M[merged1$FDR.maleVirgins_E.M<0.05 & merged1$ortho=="ortho"], na.rm=T), col=4)
rug(agdata$logFC.maleVirgins_E.M, col=1, ticksize=0.01, line=2.5)
rug(merged1$logFC.maleVirgins_E.M, col=2, ticksize=0.01, line=3.0)
# legend('topleft', inset=0.05, legend=c('DE genes', 'Findlay homologues'), pch =c(19,19), col=c(1,2), cex=0.8 )
legend('topleft', inset=0.05, legend=c('All RNAseq', 'All proteome'), pch =c(19,19), col=c(1,2), cex=0.8 )

dev.off()

# testis omnibus genes
# testdata <- read.table('~/git/dpse_mating/output/dedata/testis/LogCPM_0.05_testis.txt', header=T)
# str(testdata)
# 
# omnidatatest <- read.table('~/git/dpse_mating/input/agcp/not used/testis_genes_dpse_omnibus.txt', header=T)
# str(omnidatatest)
# 
# merged1 <- merge(testdata, omnidatatest, by.x="gene", by.y="gene", all=F )
# 
# plot(density(testdata$logFC.maleVirgins_E.M), main='Testis E vs M (all)')
# lines(density(merged1$logFC.maleVirgins_E.M), col=2)
# rug(testdata$logFC.maleVirgins_E.M, col=1, ticksize=0.01, line=2.5)
# rug(merged1$logFC.maleVirgins_E.M, col=2, ticksize=0.01, line=3.0)
# legend('topleft', inset=0.05, legend=c('All genes', 'Omnibus testis (3276)'), pch =c(19,19), col=c(1,2), cex=0.8 )
# legend('topright', inset=0.05, legend=c('W=20793000, p<0.001'),  cex=0.8 )
# wilcox.test(testdata$logFC.maleVirgins_E.M, merged1$logFC.maleVirgins_E.M)
# 
# plot(density(merged1$logFC.maleVirgins_E.M[merged1$FDR.maleVirgins_E.M<0.05]), xlim=range(-6,6), main='Testis E vs M (FDR<0.05)', col=2)
# lines(density(testdata$logFC.maleVirgins_E.M[testdata$FDR.maleVirgins_E.M<0.05]), col=1)
# rug(testdata$logFC.maleVirgins_E.M, col=2, ticksize=0.01, line=2.5)
# rug(merged1$logFC.maleVirgins_E.M, col=1, ticksize=0.01, line=3.0)
# legend('topleft', inset=0.05, legend=c('DE genes', 'Omnibus testis'), pch =c(19,19), col=c(1,2), cex=0.8 )
# 
# dev.off()
