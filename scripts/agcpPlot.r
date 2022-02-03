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

agdata <- read.table('/Users/pveltsos/Library/Mobile Documents/com~apple~CloudDocs/MyDocuments/Work/Projects/gitOut/dpseMating/dedata/agland/LogCPM_0.05_agland.txt', header=T) # data from DE analysis
str(agdata)

# merged2 <- merge(agdata, findlaydata, by.x="gene", by.y="id", all=F ) # could be used to compare to Findlay 2014 genes

merged1 <- merge(agdata, agcpsdata, by.x="gene", by.y="id", all=F )
summary(merged1)
merged1$de <- merged1$FDR.maleVirgins_E.M < 0.05  
merged1$defdr1 <- merged1$FDR.maleVirgins_E.M < 0.05 & merged1$logFC.maleVirgins_E.M > 1
merged1$defdr2 <- merged1$FDR.maleVirgins_E.M < 0.05 & merged1$logFC.maleVirgins_E.M > 2

merged1$box <- 'NA'
merged1$box[merged1]

deonly <- subset(merged1, merged1$FDR.maleVirgins_E.M < 0.05 & merged1$all=="all") 
summary(deonly)


par(mfrow=c(1,2)) 
par(mar=c(5,5,4,3), oma=c(0,0,2,0))
# boxplot(agdata$logFC.maleVirgins_E.M, merged1$logFC.maleVirgins_E.M[merged1$all=="all"], merged1$logFC.maleVirgins_E.M[merged1$agcpS=="agcpS"], merged1$logFC.maleVirgins_E.M[merged1$ortho=="ortho"], names=list('All RNAseq', 'All proteome', 'Secretome', 'SFP'), col=c(2,'gray50',3,4), ylab='logFC', cex.lab=1.3)


agdataSig <- subset(agdata, agdata$FDR.maleVirgins_E.M < 0.05)
agdataSig$box <- agdataSig$logFC.maleVirgins_E.M < 0

agdataSig$box <- factor(agdataSig$box)
mylevels <- levels(agdataSig$box)
levelProportions <- summary(agdataSig$box)/nrow(agdataSig)


boxplot(agdataSig$logFC.maleVirgins_E.M[agdataSig$FDR.maleVirgins_E.M<0.05 & agdataSig$logFC.maleVirgins_E.M>0],
		abs(agdataSig$logFC.maleVirgins_E.M[agdataSig$FDR.maleVirgins_E.M<0.05 & agdataSig$logFC.maleVirgins_E.M<0]),
		NULL, NULL, 
		names=list('All RNAseq E-up', 'All RNAseq M-up', 'Proteome E-up', 'Proteome M-up'), col=c('gray60','gray30','pink','red'), ylab='abs(logFC)', cex.lab=1.3)
		


for(i in 1:length(mylevels)){
  thislevel <- mylevels[i]
  thisvalues <- abs(agdataSig[agdataSig$box==thislevel, "logFC.maleVirgins_E.M"])
   
  # take the x-axis indices and add a jitter, proportional to the N in each level
#  myjitter <- jitter(rep(i, length(thisvalues)), amount=levelProportions[i]/2)
  myjitter <- jitter(rep(i, length(thisvalues)), amount=.25)
  points(myjitter, thisvalues, pch=20, cex=0.6, col=rgb(0,0,0,.8))    
}

m1Sig <- subset(merged1, merged1$FDR.maleVirgins_E.M < 0.05 & merged1$all=="all")
m1Sig$box <- m1Sig$logFC.maleVirgins_E.M < 0 # determines the scatter points


boxplot(NULL, NULL,
		m1Sig$logFC.maleVirgins_E.M[m1Sig$FDR.maleVirgins_E.M<0.05 & m1Sig$logFC.maleVirgins_E.M>0 & m1Sig$all=="all"],
		abs(m1Sig$logFC.maleVirgins_E.M[m1Sig$FDR.maleVirgins_E.M<0.05 & m1Sig$logFC.maleVirgins_E.M<0 & m1Sig$all=="all"]), add=T, xlab='null',
		col=c('gray60','gray30','pink','red'))

m1Sig$box <- factor(m1Sig$box)
mylevels <- levels(m1Sig$box)
levelProportions <- summary(m1Sig$box)/nrow(m1Sig)

for(i in 1:length(mylevels)){
  thislevel <- mylevels[i]
  thisvalues <- abs(m1Sig[m1Sig$box==thislevel, "logFC.maleVirgins_E.M"])
   
  # take the x-axis indices and add a jitter, proportional to the N in each level
#  myjitter <- jitter(rep(i, length(thisvalues)), amount=levelProportions[i]/2)
  myjitter <- jitter(rep(i+2, length(thisvalues)), amount=.25)
  points(myjitter, thisvalues, pch=20, cex=0.6, col=rgb(0,0,0,.8))    
}



pdf(file.path(outpath,paste('agcpPlot.pdf', sep="")), width=12, height=6) # Fig 2
par(mfrow=c(1,2)) 
plot(density(merged1$logFC.maleVirgins_E.M), xlim=range(-2,2), main='AG E vs M (all)', col=2)
lines(density(merged1$logFC.maleVirgins_E.M), col=2, lw=2)
lines(density(agdata$logFC.maleVirgins_E.M), col=1, lw=2)
lines(density(merged1$logFC.maleVirgins_E.M[merged1$agcpS=="agcpS"], na.rm=T), col=3, lw=2)
lines(density(merged1$logFC.maleVirgins_E.M[merged1$ortho=="ortho"], na.rm=T), col=4, lw=2)
# rug(agdata$logFC.maleVirgins_E.M, col=1, ticksize=0.01, line=2.5)
# rug(merged1$logFC.maleVirgins_E.M, col=2, ticksize=0.01, line=3.0)
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
