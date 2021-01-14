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

par(mfrow=c(2,3)) 

xname <- read.table("~/git/dpse_mating/output/dedata/rtract_matingE/LogCPM_0.05_rtract_matingE.txt", header=T)
str(xname)
xname$logFC.matingE_E.EEInv <- xname$logFC.matingE_E.EE*(-1)

yname <- read.table("~/git/dpse_mating/output/dedata/rtract_matingM/LogCPM_0.05_rtract_matingM.txt", header=T)
str(yname)
yname$logFC.matingM_M.MMInv <- yname$logFC.matingM_M.MM*(-1)

merged1 = merge(xname, yname, by.x="gene", by.y="gene", all=F )
head(merged1)
logFC.y <- merged1$logFC.matingE_E.EEInv
FDR.y <- merged1$FDR.matingE_E.EE
logFC.x <- merged1$logFC.matingM_M.MMInv
FDR.x <- merged1$FDR.matingM_M.MM


par(mar=c(5,5,4,3))
plot(logFC.x, logFC.y, type="n", xlab="", ylab="", main="Within line effect of mating", cex.main=1.8, cex.lab=1.3)
# plot(logFC.x, logFC.y, type="n", xlab="logFC MM/M", ylab="logFC EE/E", main="Within line effect of mating", cex.main=1.8, cex.lab=1.3)
points(logFC.x, logFC.y, pch=20, col=sgrey, cex=.8)
points(logFC.x[FDR.x < 0.05 & FDR.y > 0.05], logFC.y[FDR.x < 0.05 & FDR.y > 0.05], pch=20, col=cborange, cex=1.3)
points(logFC.x[FDR.x > 0.05 & FDR.y < 0.05], logFC.y[FDR.x > 0.05 & FDR.y < 0.05], pch=20, col=cbblue, cex=1.3)
points(logFC.x[FDR.x < 0.05 & FDR.y < 0.05], logFC.y[FDR.x < 0.05 & FDR.y < 0.05], pch=18, col=cbgreen, cex=1.3 )
legend('topleft', inset=0.05, legend=c('M/MM FDR < 0.05', 'E/EE FDR < 0.05', 'E & M FDR < 0.05'), pch =c(20,20,18), col=c(cborange, cbblue, cbgreen), cex=0.8 )
abline(v=0, lty=2, lw=1, col='grey')
abline(h=0, lty=2, lw=1, col='grey')


# nrow(subset(merged1, merged1$logFC.matingE_E.EEInv < 0.05 & merged5$FDR.matingE_E.EE < 0.05)) 				


aname <- read.table("~/git/dpse_mating/output/dedata/rtract_maleEffectOnE/LogCPM_0.05_rtract_maleEffectOnE.txt", header=T)
str(aname)

bname <- read.table("~/git/dpse_mating/output/dedata/rtract_maleEffectOnM/LogCPM_0.05_rtract_maleEffectOnM.txt", header=T)
str(bname)

merged1 = merge(xname, yname, by.x="gene", by.y="gene", all=T )
merged2 = merge(aname, bname, by.x="gene", by.y="gene", all=T )
merged3 = merge(merged1, merged2, by.x="gene", by.y="gene", all=T )

logFC.a <- merged3$logFC.maleEffectOnE_EE.EM
FDR.a <- merged3$FDR.maleEffectOnE_EE.EM
logFC.b <- merged3$logFC.maleEffectOnM_MM.ME
FDR.b <- merged3$FDR.maleEffectOnM_MM.ME
logFC.y2 <- merged3$logFC.matingE_E.EEInv
FDR.y2 <- merged3$FDR.matingE_E.EE
logFC.x2 <- merged3$logFC.matingM_M.MMInv
FDR.x2 <- merged3$FDR.matingM_M.MM



par(mar=c(5,5,4,3))
plot(logFC.x, logFC.y, type="n", xlab="", ylab="", main="Within line effect of male", cex.main=1.8, cex.lab=1.3)
# plot(logFC.x, logFC.y, type="n", xlab="logFC MM/M", ylab="logFC EE/E", main="Within line effect of male", cex.main=1.8, cex.lab=1.3)
points(logFC.x2, logFC.y2, pch=20, col=sgrey, cex=.8)
points(logFC.x2[FDR.a < 0.05 & logFC.a > 0], logFC.y2[FDR.a < 0.05 & logFC.a > 0], pch=2, col=cbpurple, cex=1.3)
points(logFC.x2[FDR.a < 0.05 & logFC.a < 0], logFC.y2[FDR.a < 0.05 & logFC.a < 0], pch=6, col=cborange, cex=1.3)
points(logFC.x2[FDR.b < 0.05 & logFC.b > 0], logFC.y2[FDR.b < 0.05 & logFC.b > 0], pch=3, col=cblblue, cex=1.3)
points(logFC.x2[FDR.b < 0.05 & logFC.b < 0], logFC.y2[FDR.b < 0.05 & logFC.b < 0], pch=4, col=cbred, cex=1.3)
legend('topleft', inset=0.07, legend=c('EE up - EM down', 'EE down - EM up', 'MM up - ME down', 'MM down - ME up'), pch=c(2,6,3,4), col=c(cbpurple, cborange, cblblue, cbred), cex=0.8)
abline(v=0, lty=2, lw=1, col='grey')
abline(h=0, lty=2, lw=1, col='grey')



cname <- read.table("~/git/dpse_mating/output/dedata/rtract_maleEffectOfE/LogCPM_0.05_rtract_maleEffectOfE.txt", header=T)
str(cname)
dname <- read.table("~/git/dpse_mating/output/dedata/rtract_maleEffectOfM/LogCPM_0.05_rtract_maleEffectOfM.txt", header=T)
str(dname)

# merged1 = merge(xname, yname, by.x="gene", by.y="gene", all=F )
merged4 = merge(cname, dname, by.x="gene", by.y="gene", all=T )
merged5 = merge(merged1, merged4, by.x="gene", by.y="gene", all=T )

# head(merged5)

# subset.merged5 <- subset(merged5, logFC.y2<0 & FDR.y2<0.05 & logFC.c<0 & FDR.c<0.05)

logFC.y2 <- merged5$logFC.matingE_E.EEInv
FDR.y2 <- merged5$FDR.matingE_E.EE
logFC.x2 <- merged5$logFC.matingM_M.MMInv
FDR.x2 <- merged5$FDR.matingM_M.MM
logFC.c <- merged5$logFC.maleEffectOfE_EE.ME
FDR.c <- merged5$FDR.maleEffectOfE_EE.ME
logFC.d <- merged5$logFC.maleEffectOfM_MM.EM
FDR.d <- merged5$FDR.maleEffectOfM_MM.EM


par(mar=c(5,5,4,3))
plot(logFC.x, logFC.y, type="n", xlab="", ylab="", main="Between line effect of male", cex.main=1.8, cex.lab=1.3)
# plot(logFC.x, logFC.y, type="n", xlab="logFC MM/M", ylab="logFC EE/E", main="Between line effect of male", cex.main=1.8, cex.lab=1.3)
points(logFC.x, logFC.y, pch=20, col=sgrey, cex=.8)
points(logFC.x2[FDR.d < 0.05 & logFC.d > 0], logFC.y2[FDR.d < 0.05 & logFC.d > 0], pch=3, col=cbpurple, cex=1.3)
points(logFC.x2[FDR.d < 0.05 & logFC.d < 0], logFC.y2[FDR.d < 0.05 & logFC.d < 0], pch=4, col=cborange, cex=1.3)
points(logFC.x2[FDR.c < 0.05 & logFC.c > 0], logFC.y2[FDR.c < 0.05 & logFC.c > 0], pch=2, col=cblblue, cex=1.3)
points(logFC.x2[FDR.c < 0.05 & logFC.c < 0], logFC.y2[FDR.c < 0.05 & logFC.c < 0], pch=6, col=cbred, cex=1.3)
legend('topleft', inset=0.07, legend=c('EE up - ME down', 'EE down - ME up', 'MM up - EM down', 'MM down - EM up'), pch=c(2,6,3,4), col=c(cblblue, cbred, cbpurple, cborange), cex=0.8)
abline(v=0, lty=2, lw=1, col='grey')
abline(h=0, lty=2, lw=1, col='grey')

																													  # fig-a	fig-c
nrow(subset(merged5, merged5$FDR.matingE_E.EE < 0.05 & merged5$FDR.matingM_M.MM < 0.05)) 							  # green 143
nrow(subset(merged5, merged5$FDR.matingE_E.EE < 0.05 & merged5$FDR.matingM_M.MM < 0.05 & FDR.d < 0.05 & logFC.d > 0)) # green purple 27
nrow(subset(merged5, merged5$FDR.matingE_E.EE < 0.05 & merged5$FDR.matingM_M.MM < 0.05 & FDR.d < 0.05 & logFC.d < 0)) # green orange 1
nrow(subset(merged5, merged5$FDR.matingE_E.EE < 0.05 & merged5$FDR.matingM_M.MM < 0.05 & FDR.c < 0.05 & logFC.c > 0)) # green lblue 0
nrow(subset(merged5, merged5$FDR.matingE_E.EE < 0.05 & merged5$FDR.matingM_M.MM < 0.05 & FDR.c < 0.05 & logFC.c < 0)) # green red 7

nrow(subset(merged5, merged5$FDR.matingE_E.EE < 0.05 & merged5$FDR.matingM_M.MM > 0.05)) 							  # blue 248
nrow(subset(merged5, merged5$FDR.matingE_E.EE < 0.05 & merged5$FDR.matingM_M.MM > 0.05 & FDR.d < 0.05 & logFC.d > 0)) # blue purple 29
nrow(subset(merged5, merged5$FDR.matingE_E.EE < 0.05 & merged5$FDR.matingM_M.MM > 0.05 & FDR.d < 0.05 & logFC.d < 0)) # blue orange 1
nrow(subset(merged5, merged5$FDR.matingE_E.EE < 0.05 & merged5$FDR.matingM_M.MM > 0.05 & FDR.c < 0.05 & logFC.c > 0)) # blue lblue 1
nrow(subset(merged5, merged5$FDR.matingE_E.EE < 0.05 & merged5$FDR.matingM_M.MM > 0.05 & FDR.c < 0.05 & logFC.c < 0)) # blue red 8

nrow(subset(merged5, merged5$FDR.matingE_E.EE > 0.05 & merged5$FDR.matingM_M.MM < 0.05)) 							  # yellow 202
nrow(subset(merged5, merged5$FDR.matingE_E.EE > 0.05 & merged5$FDR.matingM_M.MM < 0.05 & FDR.d < 0.05 & logFC.d > 0)) # yellow purple 19
nrow(subset(merged5, merged5$FDR.matingE_E.EE > 0.05 & merged5$FDR.matingM_M.MM < 0.05 & FDR.d < 0.05 & logFC.d < 0)) # yellow orange 1
nrow(subset(merged5, merged5$FDR.matingE_E.EE > 0.05 & merged5$FDR.matingM_M.MM < 0.05 & FDR.c < 0.05 & logFC.c > 0)) # yellow lblue 2
nrow(subset(merged5, merged5$FDR.matingE_E.EE > 0.05 & merged5$FDR.matingM_M.MM < 0.05 & FDR.c < 0.05 & logFC.c < 0)) # yellow red 1

nrow(subset(merged5, merged5$FDR.matingE_E.EE > 0.05 & merged5$FDR.matingM_M.MM > 0.05 & FDR.d < 0.05 & logFC.d > 0)) # gray purple 106
nrow(subset(merged5, merged5$FDR.matingE_E.EE > 0.05 & merged5$FDR.matingM_M.MM > 0.05 & FDR.d < 0.05 & logFC.d < 0)) # gray orange 213
nrow(subset(merged5, merged5$FDR.matingE_E.EE > 0.05 & merged5$FDR.matingM_M.MM > 0.05 & FDR.c < 0.05 & logFC.c > 0)) # gray lblue 11
nrow(subset(merged5, merged5$FDR.matingE_E.EE > 0.05 & merged5$FDR.matingM_M.MM > 0.05 & FDR.c < 0.05 & logFC.c < 0)) # gray red 7




points(logFC.x2[merged5$logFC.matingE_E.EE < 0 & merged5$logFC.matingM_M.MM > 0 & logFC.d < 0 & FDR.d < 0.05], logFC.y2[merged5$logFC.matingE_E.EE < 0 & merged5$logFC.matingM_M.MM > 0 & logFC.d < 0 & FDR.d < 0.05], pch=4, col=cbpurple, cex=1.3)
points(logFC.x2[merged5$logFC.matingE_E.EE > 0 & merged5$logFC.matingM_M.MM < 0 & logFC.d < 0 & FDR.d < 0.05], logFC.y2[merged5$logFC.matingE_E.EE > 0 & merged5$logFC.matingM_M.MM < 0 & logFC.d < 0 & FDR.d < 0.05], pch=4, col=cbpurple, cex=1.3)

nrow(subset(merged5, merged5$logFC.matingE_E.EE < 0 & merged5$logFC.matingM_M.MM > 0 & logFC.d > 0 & FDR.d < 0.05)) # purple top left (incongr) 10
nrow(subset(merged5, merged5$logFC.matingE_E.EE > 0 & merged5$logFC.matingM_M.MM < 0 & logFC.d > 0 & FDR.d < 0.05)) # purple bottom right (incongr) 112
nrow(subset(merged5, merged5$logFC.matingE_E.EE > 0 & merged5$logFC.matingM_M.MM > 0 & logFC.d > 0 & FDR.d < 0.05)) # purple bottom left (congr) 10
nrow(subset(merged5, merged5$logFC.matingE_E.EE < 0 & merged5$logFC.matingM_M.MM < 0 & logFC.d > 0 & FDR.d < 0.05)) # purple top right (congr) 49

prop.test(x = c(10+112, 10+49), n = c(nrow(merged5), nrow(merged5))) # purple incongr vs congr 0.012 vs 0.006, x2=21.42, df=1, p<0.001

nrow(subset(merged5, merged5$logFC.matingE_E.EE < 0 & merged5$logFC.matingM_M.MM > 0 & logFC.d < 0 & FDR.d < 0.05)) # orange top left (incongr) 5
nrow(subset(merged5, merged5$logFC.matingE_E.EE > 0 & merged5$logFC.matingM_M.MM < 0 & logFC.d < 0 & FDR.d < 0.05)) # orange bottom right (incongr) 138
nrow(subset(merged5, merged5$logFC.matingE_E.EE > 0 & merged5$logFC.matingM_M.MM > 0 & logFC.d < 0 & FDR.d < 0.05)) # orange bottom left (congr) 18
nrow(subset(merged5, merged5$logFC.matingE_E.EE < 0 & merged5$logFC.matingM_M.MM < 0 & logFC.d < 0 & FDR.d < 0.05)) # orange top right (congr) 55

prop.test(x = c(5+138, 18+55), n = c(nrow(merged5), nrow(merged5))) # orange incongr vs congr 0.014 vs 0.007, x2=22.27, df=1, p<0.001

prop.test(x = c(5+138+10+112, 10+49+18+55), n = c(nrow(merged5), nrow(merged5))) # orange + purple (EM effect) incongr vs congr 0.025 vs 0.013, x2=44.73, df=1, p<0.001

nrow(subset(merged5, merged5$logFC.matingE_E.EE < 0 & merged5$logFC.matingM_M.MM > 0 & logFC.c < 0 & FDR.c < 0.05)) # red top left (incongr) 2
nrow(subset(merged5, merged5$logFC.matingE_E.EE > 0 & merged5$logFC.matingM_M.MM < 0 & logFC.c < 0 & FDR.c < 0.05)) # red bottom right (incongr) 12
nrow(subset(merged5, merged5$logFC.matingE_E.EE > 0 & merged5$logFC.matingM_M.MM > 0 & logFC.c < 0 & FDR.c < 0.05)) # red bottom left (congr) 8
nrow(subset(merged5, merged5$logFC.matingE_E.EE < 0 & merged5$logFC.matingM_M.MM < 0 & logFC.c < 0 & FDR.c < 0.05)) # red top right (congr) 1

prop.test(x = c(2+12, 8+1), n = c(nrow(merged5), nrow(merged5))) # NS

nrow(subset(merged5, merged5$logFC.matingE_E.EE < 0 & merged5$logFC.matingM_M.MM > 0 & logFC.c > 0 & FDR.c < 0.05)) # lblue top left (incongr) 7
nrow(subset(merged5, merged5$logFC.matingE_E.EE > 0 & merged5$logFC.matingM_M.MM < 0 & logFC.c > 0 & FDR.c < 0.05)) # lblue bottom right (incongr) 1
nrow(subset(merged5, merged5$logFC.matingE_E.EE > 0 & merged5$logFC.matingM_M.MM > 0 & logFC.c > 0 & FDR.c < 0.05)) # lblue bottom left (congr) 3
nrow(subset(merged5, merged5$logFC.matingE_E.EE < 0 & merged5$logFC.matingM_M.MM < 0 & logFC.c > 0 & FDR.c < 0.05)) # lblue top right (congr) 3

prop.test(x = c(7+1, 3+3), n = c(nrow(merged5), nrow(merged5))) # NS
prop.test(x = c(7+1+2+12, 3+3+8+1), n = c(nrow(merged5), nrow(merged5))) # NS

### Make subsets of diagonals for GO analysis

#### output folder
outpath <- '~/git/dpse_mating/output/scatterplot/'
dir.create(file.path(outpath))

#### EE.ME
#### subset to only non na data
frt.EE.ME.all <- data.frame(merged5$gene, merged5$logFC.maleEffectOfE_EE.ME, merged5$FDR.maleEffectOfE_EE.ME, merged5$logFC.matingE_E.EEInv, merged5$logFC.matingM_M.MMInv)
names(frt.EE.ME.all) <- c('gene', 'logFC.maleEffectOfE_EE.ME', 'FDR.maleEffectOfE_EE.ME', 'logFC.matingE_E.EEInv', 'logFC.matingM_M.MMInv')
frt.EE.ME.all <- subset(frt.EE.ME.all, !is.na(frt.EE.ME.all$FDR.maleEffectOfE_EE.ME) & !is.na(frt.EE.ME.all$logFC.matingE_E.EEInv) & !is.na(frt.EE.ME.all$logFC.matingM_M.MMInv))

##### keep only the significant values on diagonal
frt.EE.ME.all$FDR.EE.ME.DiagExpEE <- frt.EE.ME.all$FDR.maleEffectOfE_EE.ME
frt.EE.ME.all$FDR.EE.ME.DiagExpME <- frt.EE.ME.all$FDR.maleEffectOfE_EE.ME
# points(frt.EE.ME.all$logFC.matingM_M.MMInv[frt.EE.ME.all$FDR.maleEffectOfE_EE.ME < 0.05], frt.EE.ME.all$logFC.matingE_E.EEInv[frt.EE.ME.all$FDR.maleEffectOfE_EE.ME < 0.05], pch=5, col=1, cex=1.3)
# points(frt.EE.ME.all$logFC.matingM_M.MMInv[frt.EE.ME.all$FDR.EE.ME.DiagExp < 0.05], frt.EE.ME.all$logFC.matingE_E.EEInv[frt.EE.ME.all$FDR.EE.ME.DiagExp < 0.05], pch=5, col=2, cex=1.3)

frt.EE.ME.all$FDR.EE.ME.DiagExpEE[(frt.EE.ME.all$logFC.maleEffectOfE_EE.ME > 0 & frt.EE.ME.all$logFC.matingE_E.EEInv > 0 & frt.EE.ME.all$logFC.matingM_M.MMInv < 0 ) | (frt.EE.ME.all$logFC.maleEffectOfE_EE.ME > 0 & frt.EE.ME.all$logFC.matingE_E.EEInv < 0 & frt.EE.ME.all$logFC.matingM_M.MMInv > 0)] <- 1 

##### Make diagonal data frame for GO output
frt.EE.ME.DiagExpEE <- data.frame(frt.EE.ME.all$gene, frt.EE.ME.all$FDR.EE.ME.DiagExpEE)
names(frt.EE.ME.DiagExpEE) <- c('gene', 'FDR.EE.ME.DiagExp')

gopathDiagExpEE <- paste(outpath, 'frtEEMEDiagExpEE/', sep="")
dir.create(file.path(gopathDiagExpEE))
write.table(frt.EE.ME.DiagExpEE, file=file.path(gopathDiagExpEE, "GO_pvalues.txt"), quote=F, row.names=F, sep="\t")

frt.EE.ME.all$FDR.EE.ME.DiagExpME[(frt.EE.ME.all$logFC.maleEffectOfE_EE.ME < 0 & frt.EE.ME.all$logFC.matingE_E.EEInv > 0 & frt.EE.ME.all$logFC.matingM_M.MMInv < 0 ) | (frt.EE.ME.all$logFC.maleEffectOfE_EE.ME < 0 & frt.EE.ME.all$logFC.matingE_E.EEInv < 0 & frt.EE.ME.all$logFC.matingM_M.MMInv > 0)] <- 1 

##### Make diagonal data frame for GO output
frt.EE.ME.DiagExpME <- data.frame(frt.EE.ME.all$gene, frt.EE.ME.all$FDR.EE.ME.DiagExpME)
names(frt.EE.ME.DiagExpME) <- c('gene', 'FDR.EE.ME.DiagExp')

gopathDiagExpME <- paste(outpath, 'frtEEMEDiagExpME/', sep="")
dir.create(file.path(gopathDiagExpME))
write.table(frt.EE.ME.DiagExpME, file=file.path(gopathDiagExpME, "GO_pvalues.txt"), quote=F, row.names=F, sep="\t")


##### keep only the significant values on NOT exp diagonal
frt.EE.ME.all$FDR.EE.ME.DiagNOTExpEE <- frt.EE.ME.all$FDR.maleEffectOfE_EE.ME
frt.EE.ME.all$FDR.EE.ME.DiagNOTExpME <- frt.EE.ME.all$FDR.maleEffectOfE_EE.ME

frt.EE.ME.all$FDR.EE.ME.DiagNOTExpEE[(frt.EE.ME.all$logFC.maleEffectOfE_EE.ME > 0 & frt.EE.ME.all$logFC.matingE_E.EEInv > 0 & frt.EE.ME.all$logFC.matingM_M.MMInv > 0 ) | (frt.EE.ME.all$logFC.maleEffectOfE_EE.ME > 0 & frt.EE.ME.all$logFC.matingE_E.EEInv < 0 & frt.EE.ME.all$logFC.matingM_M.MMInv < 0)] <- 1 

##### Make diagonal data frame for GO output
frt.EE.ME.DiagNOTExpEE <- data.frame(frt.EE.ME.all$gene, frt.EE.ME.all$FDR.EE.ME.DiagNOTExpEE)
names(frt.EE.ME.DiagNOTExpEE) <- c('gene', 'FDR.EE.ME.DiagNOTExp')

gopathDiagNOTExpEE <- paste(outpath, 'frtEEMEDiagNOTExpEE/', sep="")
dir.create(file.path(gopathDiagNOTExpEE))
write.table(frt.EE.ME.DiagNOTExpEE, file=file.path(gopathDiagNOTExpEE, "GO_pvalues.txt"), quote=F, row.names=F, sep="\t")

frt.EE.ME.all$FDR.EE.ME.DiagNOTExpEE[(frt.EE.ME.all$logFC.maleEffectOfE_EE.ME < 0 & frt.EE.ME.all$logFC.matingE_E.EEInv > 0 & frt.EE.ME.all$logFC.matingM_M.MMInv > 0 ) | (frt.EE.ME.all$logFC.maleEffectOfE_EE.ME < 0 & frt.EE.ME.all$logFC.matingE_E.EEInv < 0 & frt.EE.ME.all$logFC.matingM_M.MMInv < 0)] <- 1 


##### Make diagonal data frame for GO output
frt.EE.ME.DiagNOTExpME <- data.frame(frt.EE.ME.all$gene, frt.EE.ME.all$FDR.EE.ME.DiagNOTExpME)
names(frt.EE.ME.DiagNOTExpME) <- c('gene', 'FDR.EE.ME.DiagNOTExp')

gopathDiagNOTExpME <- paste(outpath, 'frtEEMEDiagNOTExpME/', sep="")
dir.create(file.path(gopathDiagNOTExpME))
write.table(frt.EE.ME.DiagNOTExpME, file=file.path(gopathDiagNOTExpME, "GO_pvalues.txt"), quote=F, row.names=F, sep="\t")



#### MM.EM
#### subset to only non na data
frt.MM.EM.all <- data.frame(merged5$gene, merged5$logFC.maleEffectOfM_MM.EM, merged5$FDR.maleEffectOfM_MM.EM, merged5$logFC.matingE_E.EEInv, merged5$logFC.matingM_M.MMInv)
names(frt.MM.EM.all) <- c('gene', 'logFC.maleEffectOfM_MM.EM', 'FDR.maleEffectOfM_MM.EM', 'logFC.matingE_E.EEInv', 'logFC.matingM_M.MMInv')
frt.MM.EM.all <- subset(frt.MM.EM.all, !is.na(frt.MM.EM.all$FDR.maleEffectOfM_MM.EM) & !is.na(frt.MM.EM.all$logFC.matingE_E.EEInv) & !is.na(frt.MM.EM.all$logFC.matingM_M.MMInv))

##### keep only the significant values on diagonal
frt.MM.EM.all$FDR.MM.EM.DiagExpEE <- frt.MM.EM.all$FDR.maleEffectOfM_MM.EM
frt.MM.EM.all$FDR.MM.EM.DiagExpME <- frt.MM.EM.all$FDR.maleEffectOfM_MM.EM
# points(frt.MM.EM.all$logFC.matingM_M.MMInv[frt.MM.EM.all$FDR.maleEffectOfM_MM.EM < 0.05], frt.MM.EM.all$logFC.matingE_E.EEInv[frt.MM.EM.all$FDR.maleEffectOfM_MM.EM < 0.05], pch=5, col=1, cex=1.3)
# points(frt.MM.EM.all$logFC.matingM_M.MMInv[frt.MM.EM.all$FDR.MM.EM.DiagExp < 0.05], frt.MM.EM.all$logFC.matingE_E.EEInv[frt.MM.EM.all$FDR.MM.EM.DiagExp < 0.05], pch=5, col=2, cex=1.3)

frt.MM.EM.all$FDR.MM.EM.DiagExpEE[(frt.MM.EM.all$logFC.maleEffectOfM_MM.EM > 0 & frt.MM.EM.all$logFC.matingE_E.EEInv > 0 & frt.MM.EM.all$logFC.matingM_M.MMInv < 0 ) | (frt.MM.EM.all$logFC.maleEffectOfM_MM.EM > 0 & frt.MM.EM.all$logFC.matingE_E.EEInv < 0 & frt.MM.EM.all$logFC.matingM_M.MMInv > 0)] <- 1 

##### Make diagonal data frame for GO output
frt.MM.EM.DiagExpEE <- data.frame(frt.MM.EM.all$gene, frt.MM.EM.all$FDR.MM.EM.DiagExpEE)
names(frt.MM.EM.DiagExpEE) <- c('gene', 'FDR.MM.EM.DiagExp')

gopathDiagExpEE <- paste(outpath, 'frtMMEMDiagExpEE/', sep="")
dir.create(file.path(gopathDiagExpEE))
write.table(frt.MM.EM.DiagExpEE, file=file.path(gopathDiagExpEE, "GO_pvalues.txt"), quote=F, row.names=F, sep="\t")

frt.MM.EM.all$FDR.MM.EM.DiagExpME[(frt.MM.EM.all$logFC.maleEffectOfM_MM.EM < 0 & frt.MM.EM.all$logFC.matingE_E.EEInv > 0 & frt.MM.EM.all$logFC.matingM_M.MMInv < 0 ) | (frt.MM.EM.all$logFC.maleEffectOfM_MM.EM < 0 & frt.MM.EM.all$logFC.matingE_E.EEInv < 0 & frt.MM.EM.all$logFC.matingM_M.MMInv > 0)] <- 1 

##### Make diagonal data frame for GO output
frt.MM.EM.DiagExpME <- data.frame(frt.MM.EM.all$gene, frt.MM.EM.all$FDR.MM.EM.DiagExpME)
names(frt.MM.EM.DiagExpME) <- c('gene', 'FDR.MM.EM.DiagExp')

gopathDiagExpME <- paste(outpath, 'frtMMEMDiagExpME/', sep="")
dir.create(file.path(gopathDiagExpME))
write.table(frt.MM.EM.DiagExpME, file=file.path(gopathDiagExpME, "GO_pvalues.txt"), quote=F, row.names=F, sep="\t")


##### keep only the significant values on NOT exp diagonal
frt.MM.EM.all$FDR.MM.EM.DiagNOTExpEE <- frt.MM.EM.all$FDR.maleEffectOfM_MM.EM
frt.MM.EM.all$FDR.MM.EM.DiagNOTExpME <- frt.MM.EM.all$FDR.maleEffectOfM_MM.EM

frt.MM.EM.all$FDR.MM.EM.DiagNOTExpEE[(frt.MM.EM.all$logFC.maleEffectOfM_MM.EM > 0 & frt.MM.EM.all$logFC.matingE_E.EEInv > 0 & frt.MM.EM.all$logFC.matingM_M.MMInv > 0 ) | (frt.MM.EM.all$logFC.maleEffectOfM_MM.EM > 0 & frt.MM.EM.all$logFC.matingE_E.EEInv < 0 & frt.MM.EM.all$logFC.matingM_M.MMInv < 0)] <- 1 

##### Make diagonal data frame for GO output
frt.MM.EM.DiagNOTExpEE <- data.frame(frt.MM.EM.all$gene, frt.MM.EM.all$FDR.MM.EM.DiagNOTExpEE)
names(frt.MM.EM.DiagNOTExpEE) <- c('gene', 'FDR.MM.EM.DiagNOTExp')

gopathDiagNOTExpEE <- paste(outpath, 'frtMMEMDiagNOTExpEE/', sep="")
dir.create(file.path(gopathDiagNOTExpEE))
write.table(frt.MM.EM.DiagNOTExpEE, file=file.path(gopathDiagNOTExpEE, "GO_pvalues.txt"), quote=F, row.names=F, sep="\t")

frt.MM.EM.all$FDR.MM.EM.DiagNOTExpEE[(frt.MM.EM.all$logFC.maleEffectOfM_MM.EM < 0 & frt.MM.EM.all$logFC.matingE_E.EEInv > 0 & frt.MM.EM.all$logFC.matingM_M.MMInv > 0 ) | (frt.MM.EM.all$logFC.maleEffectOfM_MM.EM < 0 & frt.MM.EM.all$logFC.matingE_E.EEInv < 0 & frt.MM.EM.all$logFC.matingM_M.MMInv < 0)] <- 1 


##### Make diagonal data frame for GO output
frt.MM.EM.DiagNOTExpME <- data.frame(frt.MM.EM.all$gene, frt.MM.EM.all$FDR.MM.EM.DiagNOTExpME)
names(frt.MM.EM.DiagNOTExpME) <- c('gene', 'FDR.MM.EM.DiagNOTExp')

gopathDiagNOTExpME <- paste(outpath, 'frtMMEMDiagNOTExpME/', sep="")
dir.create(file.path(gopathDiagNOTExpME))
write.table(frt.MM.EM.DiagNOTExpME, file=file.path(gopathDiagNOTExpME, "GO_pvalues.txt"), quote=F, row.names=F, sep="\t")






## ovaries

xname <- read.table("~/git/dpse_mating/output/dedata/ovary_matingE_rmEEov4/LogCPM_0.05_ovary_matingE_rmEEov4.txt", header=T)
str(xname)
xname$logFC.matingE_E.EEInv <- xname$logFC.matingE_E.EE*(-1)
yname <- read.table("~/git/dpse_mating/output/dedata/ovary_matingM/LogCPM_0.05_ovary_matingM.txt", header=T)
str(yname)
yname$logFC.matingM_M.MMInv <- yname$logFC.matingM_M.MM*(-1)

merged1 = merge(xname, yname, by.x="gene", by.y="gene", all=F )
logFC.y <- merged1$logFC.matingE_E.EEInv
FDR.y <- merged1$FDR.matingE_E.EE
logFC.x <- merged1$logFC.matingM_M.MMInv
FDR.x <- merged1$FDR.matingM_M.MM

par(mar=c(5,5,4,3))
plot(logFC.x, logFC.y, type="n", xlab="", ylab="", main="Within line effect of mating", cex.main=1.8, cex.lab=1.3)
# plot(logFC.x, logFC.y, type="n", xlab="logFC MM/M", ylab="logFC EE/E", main="Within line effect of mating", cex.main=1.8, cex.lab=1.3)
points(logFC.x[FDR.x > 0.05 & FDR.x > 0.05], logFC.y[FDR.x > 0.05 & FDR.x > 0.05], pch=20, col=sgrey, cex=.8)
points(logFC.x[FDR.x < 0.05 & FDR.y > 0.05], logFC.y[FDR.x < 0.05 & FDR.y > 0.05], pch=20, col=cborange, cex=1.3)
points(logFC.x[FDR.x > 0.05 & FDR.y < 0.05], logFC.y[FDR.x > 0.05 & FDR.y < 0.05], pch=20, col=cbblue, cex=1.3)
points(logFC.x[FDR.x < 0.05 & FDR.y < 0.05], logFC.y[FDR.x < 0.05 & FDR.y < 0.05], pch=18, col=cbgreen, cex=1.3 )
legend('topleft', inset=0.05, legend=c('FDR MM/M < 0.05', 'FDR EE/E < 0.05', 'FDR < 0.05 (both)'), pch =c(20,20,18), col=c(cborange, cbblue, cbgreen), cex=0.8 )
abline(v=0, lty=2, lw=1, col='grey')
abline(h=0, lty=2, lw=1, col='grey')



aname <- read.table("~/git/dpse_mating/output/dedata/ovary_maleEffectOnE/LogCPM_0.05_ovary_maleEffectOnE.txt", header=T)
str(aname)

bname <- read.table("~/git/dpse_mating/output/dedata/ovary_maleEffectOnM/LogCPM_0.05_ovary_maleEffectOnM.txt", header=T)
str(bname)

merged1 = merge(xname, yname, by.x="gene", by.y="gene", all=T )
merged2 = merge(aname, bname, by.x="gene", by.y="gene", all=T )
merged3 = merge(merged1, merged2, by.x="gene", by.y="gene", all=T)

logFC.a <- merged3$logFC.maleEffectOnE_EE.EM
FDR.a <- merged3$FDR.maleEffectOnE_EE.EM
logFC.b <- merged3$logFC.maleEffectOnM_MM.ME
FDR.b <- merged3$FDR.maleEffectOnM_MM.ME
logFC.y2 <- merged3$logFC.matingE_E.EEInv
FDR.y2 <- merged3$FDR.matingE_E.EE
logFC.x2 <- merged3$logFC.matingM_M.MMInv
FDR.x2 <- merged3$FDR.matingM_M.MM



par(mar=c(5,5,4,3))
plot(logFC.x, logFC.y, type="n", xlab="", ylab="", main="Within line effect of male", cex.main=1.8, cex.lab=1.3)
# plot(logFC.x, logFC.y, type="n", xlab="logFC MM/M", ylab="logFC EE/E", main="Within line effect of male", cex.main=1.8, cex.lab=1.3)
points(logFC.x, logFC.y, pch=20, col=sgrey, cex=.8)
points(logFC.x2[FDR.a < 0.05 & logFC.a > 0], logFC.y2[FDR.a < 0.05 & logFC.a > 0], pch=2, col=cbpurple, cex=1.3)
points(logFC.x2[FDR.a < 0.05 & logFC.a < 0], logFC.y2[FDR.a < 0.05 & logFC.a < 0], pch=6, col=cborange, cex=1.3)
points(logFC.x2[FDR.b < 0.05 & logFC.b > 0], logFC.y2[FDR.b < 0.05 & logFC.b > 0], pch=3, col=cblblue, cex=1.3)
points(logFC.x2[FDR.b < 0.05 & logFC.b < 0], logFC.y2[FDR.b < 0.05 & logFC.b < 0], pch=4, col=cbred, cex=1.3)
legend('topleft', inset=0.07, legend=c('EE up - EM down', 'EE down - EM up', 'MM up - ME down', 'MM down - ME up'), pch=c(2,6,3,4), col=c(cbpurple, cborange, cblblue, cbred), cex=0.8)
abline(v=0, lty=2, lw=1, col='grey')
abline(h=0, lty=2, lw=1, col='grey')





cname <- read.table("~/git/dpse_mating/output/dedata/ovary_maleEffectOfE/LogCPM_0.05_ovary_maleEffectOfE.txt", header=T)
str(cname)
dname <- read.table("~/git/dpse_mating/output/dedata/ovary_maleEffectOfM/LogCPM_0.05_ovary_maleEffectOfM.txt", header=T)
str(dname)

merged1 = merge(xname, yname, by.x="gene", by.y="gene", all=T )
merged4 = merge(cname, dname, by.x="gene", by.y="gene", all=T )
merged5 = merge(merged1, merged4, by.x="gene", by.y="gene", all=T )

logFC.y2 <- merged5$logFC.matingE_E.EEInv
FDR.y2 <- merged5$FDR.matingE_E.EE
logFC.x2 <- merged5$logFC.matingM_M.MMInv
FDR.x2 <- merged5$FDR.matingM_M.MM
logFC.c <- merged5$logFC.maleEffectOfE_EE.ME
FDR.c <- merged5$FDR.maleEffectOfE_EE.ME
logFC.d <- merged5$logFC.maleEffectOfM_MM.EM
FDR.d <- merged5$FDR.maleEffectOfM_MM.EM


par(mar=c(5,5,4,3))
plot(logFC.x, logFC.y, type="n", xlab="", ylab="", main="Between line effect of male", cex.main=1.8, cex.lab=1.3)
# plot(logFC.x, logFC.y, type="n", xlab="logFC MM/M", ylab="logFC EE/E", main="Between line effect of male", cex.main=1.8, cex.lab=1.3)
points(logFC.x, logFC.y, pch=20, col=sgrey, cex=.8)
points(logFC.x2[FDR.c < 0.05 & logFC.c > 0], logFC.y2[FDR.c < 0.05 & logFC.c > 0], pch=2, col=cblblue, cex=1.3)
points(logFC.x2[FDR.c < 0.05 & logFC.c < 0], logFC.y2[FDR.c < 0.05 & logFC.c < 0], pch=6, col=cbred, cex=1.3)
points(logFC.x2[FDR.d < 0.05 & logFC.d > 0], logFC.y2[FDR.d < 0.05 & logFC.d > 0], pch=3, col=cbpurple, cex=1.3)
points(logFC.x2[FDR.d < 0.05 & logFC.d < 0], logFC.y2[FDR.d < 0.05 & logFC.d < 0], pch=4, col=cborange, cex=1.3)
legend('topleft', inset=0.07, legend=c('EE up - ME down', 'EE down - ME up', 'MM up - EM down', 'MM down - EM up'), pch=c(2,6,3,4), col=c(cblblue, cbred, cbpurple, cborange), cex=0.8)
abline(v=0, lty=2, lw=1, col='grey')
abline(h=0, lty=2, lw=1, col='grey')

dev.copy(pdf,file.path(paste(outpath, "fig3.pdf", sep="")), width=15, height=10)
dev.off()



# Calculate number of congruent and incongruent genes
points(logFC.x2[merged5$logFC.matingE_E.EE < 0 & merged5$logFC.matingM_M.MM > 0 & logFC.d < 0 & FDR.d < 0.05], logFC.y2[merged5$logFC.matingE_E.EE < 0 & merged5$logFC.matingM_M.MM > 0 & logFC.d < 0 & FDR.d < 0.05], pch=4, col=cbpurple, cex=1.3)
points(logFC.x2[merged5$logFC.matingE_E.EE > 0 & merged5$logFC.matingM_M.MM < 0 & logFC.d < 0 & FDR.d < 0.05], logFC.y2[merged5$logFC.matingE_E.EE > 0 & merged5$logFC.matingM_M.MM < 0 & logFC.d < 0 & FDR.d < 0.05], pch=4, col=cbpurple, cex=1.3)

nrow(subset(merged5, merged5$logFC.matingE_E.EE < 0 & merged5$logFC.matingM_M.MM > 0 & logFC.d > 0 & FDR.d < 0.05)) # purple top left (incongr) 6
nrow(subset(merged5, merged5$logFC.matingE_E.EE > 0 & merged5$logFC.matingM_M.MM < 0 & logFC.d > 0 & FDR.d < 0.05)) # purple bottom right (incongr) 13
nrow(subset(merged5, merged5$logFC.matingE_E.EE > 0 & merged5$logFC.matingM_M.MM > 0 & logFC.d > 0 & FDR.d < 0.05)) # purple bottom left (congr) 9
nrow(subset(merged5, merged5$logFC.matingE_E.EE < 0 & merged5$logFC.matingM_M.MM < 0 & logFC.d > 0 & FDR.d < 0.05)) # purple top right (congr) 30

prop.test(x = c(6+13, 9+30), n = c(nrow(merged5), nrow(merged5))) # purple incongr vs congr 0.002 vs 0.004, x2=6.25, df=1, p=0.013

nrow(subset(merged5, merged5$logFC.matingE_E.EE < 0 & merged5$logFC.matingM_M.MM > 0 & logFC.d < 0 & FDR.d < 0.05)) # orange top left (incongr) 12
nrow(subset(merged5, merged5$logFC.matingE_E.EE > 0 & merged5$logFC.matingM_M.MM < 0 & logFC.d < 0 & FDR.d < 0.05)) # orange bottom right (incongr) 11
nrow(subset(merged5, merged5$logFC.matingE_E.EE > 0 & merged5$logFC.matingM_M.MM > 0 & logFC.d < 0 & FDR.d < 0.05)) # orange bottom left (congr) 30
nrow(subset(merged5, merged5$logFC.matingE_E.EE < 0 & merged5$logFC.matingM_M.MM < 0 & logFC.d < 0 & FDR.d < 0.05)) # orange top right (congr) 9

prop.test(x = c(12+11, 30+9), n = c(nrow(merged5), nrow(merged5))) # orange incongr vs congr 0.0026 vs 0.0045, x2=3.64, df=1, p<0.056

prop.test(x = c(6+13+12+11, 9+30+9+30), n = c(nrow(merged5), nrow(merged5))) # orange + purple (EM effect) incongr vs congr 0.005 vs 0.009, x2=10.28, df=1, p<0.001

nrow(subset(merged5, merged5$logFC.matingE_E.EE < 0 & merged5$logFC.matingM_M.MM > 0 & logFC.c < 0 & FDR.c < 0.05)) # red top left (incongr) 2
nrow(subset(merged5, merged5$logFC.matingE_E.EE > 0 & merged5$logFC.matingM_M.MM < 0 & logFC.c < 0 & FDR.c < 0.05)) # red bottom right (incongr) 42
nrow(subset(merged5, merged5$logFC.matingE_E.EE > 0 & merged5$logFC.matingM_M.MM > 0 & logFC.c < 0 & FDR.c < 0.05)) # red bottom left (congr) 46
nrow(subset(merged5, merged5$logFC.matingE_E.EE < 0 & merged5$logFC.matingM_M.MM < 0 & logFC.c < 0 & FDR.c < 0.05)) # red top right (congr) 10

prop.test(x = c(2+42, 46+10), n = c(nrow(merged5), nrow(merged5))) # NS

nrow(subset(merged5, merged5$logFC.matingE_E.EE < 0 & merged5$logFC.matingM_M.MM > 0 & logFC.c > 0 & FDR.c < 0.05)) # lblue top left (incongr) 62
nrow(subset(merged5, merged5$logFC.matingE_E.EE > 0 & merged5$logFC.matingM_M.MM < 0 & logFC.c > 0 & FDR.c < 0.05)) # lblue bottom right (incongr) 9
nrow(subset(merged5, merged5$logFC.matingE_E.EE > 0 & merged5$logFC.matingM_M.MM > 0 & logFC.c > 0 & FDR.c < 0.05)) # lblue bottom left (congr) 14
nrow(subset(merged5, merged5$logFC.matingE_E.EE < 0 & merged5$logFC.matingM_M.MM < 0 & logFC.c > 0 & FDR.c < 0.05)) # lblue top right (congr) 100

prop.test(x = c(62+9, 14+100), n = c(nrow(merged5), nrow(merged5))) # lblue incongr vs congr 0.008 vs 0.13, x2=9.64, p=0.002
prop.test(x = c(2+42+62+9, 46+10+14+100), n = c(nrow(merged5), nrow(merged5))) # red + lblue incongr vs congr 0.013 vs 0.019, x2=10.4, p=0.0012



																													  # fig-d	fig-f
nrow(subset(merged5, merged5$FDR.matingE_E.EE < 0.05 & merged5$FDR.matingM_M.MM < 0.05)) 							  # green 1
nrow(subset(merged5, merged5$FDR.matingE_E.EE < 0.05 & merged5$FDR.matingM_M.MM < 0.05 & FDR.d < 0.05 & logFC.d > 0)) # green purple 0
nrow(subset(merged5, merged5$FDR.matingE_E.EE < 0.05 & merged5$FDR.matingM_M.MM < 0.05 & FDR.d < 0.05 & logFC.d < 0)) # green orange 1
nrow(subset(merged5, merged5$FDR.matingE_E.EE < 0.05 & merged5$FDR.matingM_M.MM < 0.05 & FDR.c < 0.05 & logFC.c > 0)) # green lblue 0
nrow(subset(merged5, merged5$FDR.matingE_E.EE < 0.05 & merged5$FDR.matingM_M.MM < 0.05 & FDR.c < 0.05 & logFC.c < 0)) # green red 1

nrow(subset(merged5, merged5$FDR.matingE_E.EE < 0.05 & merged5$FDR.matingM_M.MM > 0.05)) 							  # blue 21
nrow(subset(merged5, merged5$FDR.matingE_E.EE < 0.05 & merged5$FDR.matingM_M.MM > 0.05 & FDR.d < 0.05 & logFC.d > 0)) # blue purple 4
nrow(subset(merged5, merged5$FDR.matingE_E.EE < 0.05 & merged5$FDR.matingM_M.MM > 0.05 & FDR.d < 0.05 & logFC.d < 0)) # blue orange 0
nrow(subset(merged5, merged5$FDR.matingE_E.EE < 0.05 & merged5$FDR.matingM_M.MM > 0.05 & FDR.c < 0.05 & logFC.c > 0)) # blue lblue 0
nrow(subset(merged5, merged5$FDR.matingE_E.EE < 0.05 & merged5$FDR.matingM_M.MM > 0.05 & FDR.c < 0.05 & logFC.c < 0)) # blue red 0

nrow(subset(merged5, merged5$FDR.matingE_E.EE > 0.05 & merged5$FDR.matingM_M.MM < 0.05)) 							  # yellow 0
nrow(subset(merged5, merged5$FDR.matingE_E.EE > 0.05 & merged5$FDR.matingM_M.MM < 0.05 & FDR.d < 0.05 & logFC.d > 0)) # yellow purple 0
nrow(subset(merged5, merged5$FDR.matingE_E.EE > 0.05 & merged5$FDR.matingM_M.MM < 0.05 & FDR.d < 0.05 & logFC.d < 0)) # yellow orange 0
nrow(subset(merged5, merged5$FDR.matingE_E.EE > 0.05 & merged5$FDR.matingM_M.MM < 0.05 & FDR.c < 0.05 & logFC.c > 0)) # yellow lblue 0
nrow(subset(merged5, merged5$FDR.matingE_E.EE > 0.05 & merged5$FDR.matingM_M.MM < 0.05 & FDR.c < 0.05 & logFC.c < 0)) # yellow red 0

nrow(subset(merged5, merged5$FDR.matingE_E.EE > 0.05 & merged5$FDR.matingM_M.MM > 0.05 & FDR.d < 0.05 & logFC.d > 0)) # gray purple 54
nrow(subset(merged5, merged5$FDR.matingE_E.EE > 0.05 & merged5$FDR.matingM_M.MM > 0.05 & FDR.d < 0.05 & logFC.d < 0)) # gray orange 62
nrow(subset(merged5, merged5$FDR.matingE_E.EE > 0.05 & merged5$FDR.matingM_M.MM > 0.05 & FDR.c < 0.05 & logFC.c > 0)) # gray lblue 185
nrow(subset(merged5, merged5$FDR.matingE_E.EE > 0.05 & merged5$FDR.matingM_M.MM > 0.05 & FDR.c < 0.05 & logFC.c < 0)) # gray red 99







#### EE.ME
#### subset to only non na data
ov.EE.ME.all <- data.frame(merged5$gene, merged5$logFC.maleEffectOfE_EE.ME, merged5$FDR.maleEffectOfE_EE.ME, merged5$logFC.matingE_E.EEInv, merged5$logFC.matingM_M.MMInv)
names(ov.EE.ME.all) <- c('gene', 'logFC.maleEffectOfE_EE.ME', 'FDR.maleEffectOfE_EE.ME', 'logFC.matingE_E.EEInv', 'logFC.matingM_M.MMInv')
ov.EE.ME.all <- subset(ov.EE.ME.all, !is.na(ov.EE.ME.all$FDR.maleEffectOfE_EE.ME) & !is.na(ov.EE.ME.all$logFC.matingE_E.EEInv) & !is.na(ov.EE.ME.all$logFC.matingM_M.MMInv))

##### keep only the significant values on diagonal
ov.EE.ME.all$FDR.EE.ME.DiagExpEE <- ov.EE.ME.all$FDR.maleEffectOfE_EE.ME
ov.EE.ME.all$FDR.EE.ME.DiagExpME <- ov.EE.ME.all$FDR.maleEffectOfE_EE.ME
# points(ov.EE.ME.all$logFC.matingM_M.MMInv[ov.EE.ME.all$FDR.maleEffectOfE_EE.ME < 0.05], ov.EE.ME.all$logFC.matingE_E.EEInv[ov.EE.ME.all$FDR.maleEffectOfE_EE.ME < 0.05], pch=5, col=1, cex=1.3)
# points(ov.EE.ME.all$logFC.matingM_M.MMInv[ov.EE.ME.all$FDR.EE.ME.DiagExp < 0.05], ov.EE.ME.all$logFC.matingE_E.EEInv[ov.EE.ME.all$FDR.EE.ME.DiagExp < 0.05], pch=5, col=2, cex=1.3)

ov.EE.ME.all$FDR.EE.ME.DiagExpEE[(ov.EE.ME.all$logFC.maleEffectOfE_EE.ME > 0 & ov.EE.ME.all$logFC.matingE_E.EEInv > 0 & ov.EE.ME.all$logFC.matingM_M.MMInv < 0 ) | (ov.EE.ME.all$logFC.maleEffectOfE_EE.ME > 0 & ov.EE.ME.all$logFC.matingE_E.EEInv < 0 & ov.EE.ME.all$logFC.matingM_M.MMInv > 0)] <- 1 

##### Make diagonal data frame for GO output
ov.EE.ME.DiagExpEE <- data.frame(ov.EE.ME.all$gene, ov.EE.ME.all$FDR.EE.ME.DiagExpEE)
names(ov.EE.ME.DiagExpEE) <- c('gene', 'FDR.EE.ME.DiagExp')

gopathDiagExpEE <- paste(outpath, 'ovEEMEDiagExpEE/', sep="")
dir.create(file.path(gopathDiagExpEE))
write.table(ov.EE.ME.DiagExpEE, file=file.path(gopathDiagExpEE, "GO_pvalues.txt"), quote=F, row.names=F, sep="\t")

ov.EE.ME.all$FDR.EE.ME.DiagExpME[(ov.EE.ME.all$logFC.maleEffectOfE_EE.ME < 0 & ov.EE.ME.all$logFC.matingE_E.EEInv > 0 & ov.EE.ME.all$logFC.matingM_M.MMInv < 0 ) | (ov.EE.ME.all$logFC.maleEffectOfE_EE.ME < 0 & ov.EE.ME.all$logFC.matingE_E.EEInv < 0 & ov.EE.ME.all$logFC.matingM_M.MMInv > 0)] <- 1 

##### Make diagonal data frame for GO output
ov.EE.ME.DiagExpME <- data.frame(ov.EE.ME.all$gene, ov.EE.ME.all$FDR.EE.ME.DiagExpME)
names(ov.EE.ME.DiagExpME) <- c('gene', 'FDR.EE.ME.DiagExp')

gopathDiagExpME <- paste(outpath, 'ovEEMEDiagExpME/', sep="")
dir.create(file.path(gopathDiagExpME))
write.table(ov.EE.ME.DiagExpME, file=file.path(gopathDiagExpME, "GO_pvalues.txt"), quote=F, row.names=F, sep="\t")


##### keep only the significant values on NOT exp diagonal
ov.EE.ME.all$FDR.EE.ME.DiagNOTExpEE <- ov.EE.ME.all$FDR.maleEffectOfE_EE.ME
ov.EE.ME.all$FDR.EE.ME.DiagNOTExpME <- ov.EE.ME.all$FDR.maleEffectOfE_EE.ME

ov.EE.ME.all$FDR.EE.ME.DiagNOTExpEE[(ov.EE.ME.all$logFC.maleEffectOfE_EE.ME > 0 & ov.EE.ME.all$logFC.matingE_E.EEInv > 0 & ov.EE.ME.all$logFC.matingM_M.MMInv > 0 ) | (ov.EE.ME.all$logFC.maleEffectOfE_EE.ME > 0 & ov.EE.ME.all$logFC.matingE_E.EEInv < 0 & ov.EE.ME.all$logFC.matingM_M.MMInv < 0)] <- 1 

##### Make diagonal data frame for GO output
ov.EE.ME.DiagNOTExpEE <- data.frame(ov.EE.ME.all$gene, ov.EE.ME.all$FDR.EE.ME.DiagNOTExpEE)
names(ov.EE.ME.DiagNOTExpEE) <- c('gene', 'FDR.EE.ME.DiagNOTExp')

gopathDiagNOTExpEE <- paste(outpath, 'ovEEMEDiagNOTExpEE/', sep="")
dir.create(file.path(gopathDiagNOTExpEE))
write.table(ov.EE.ME.DiagNOTExpEE, file=file.path(gopathDiagNOTExpEE, "GO_pvalues.txt"), quote=F, row.names=F, sep="\t")

ov.EE.ME.all$FDR.EE.ME.DiagNOTExpEE[(ov.EE.ME.all$logFC.maleEffectOfE_EE.ME < 0 & ov.EE.ME.all$logFC.matingE_E.EEInv > 0 & ov.EE.ME.all$logFC.matingM_M.MMInv > 0 ) | (ov.EE.ME.all$logFC.maleEffectOfE_EE.ME < 0 & ov.EE.ME.all$logFC.matingE_E.EEInv < 0 & ov.EE.ME.all$logFC.matingM_M.MMInv < 0)] <- 1 


##### Make diagonal data frame for GO output
ov.EE.ME.DiagNOTExpME <- data.frame(ov.EE.ME.all$gene, ov.EE.ME.all$FDR.EE.ME.DiagNOTExpME)
names(ov.EE.ME.DiagNOTExpME) <- c('gene', 'FDR.EE.ME.DiagNOTExp')

gopathDiagNOTExpME <- paste(outpath, 'ovEEMEDiagNOTExpME/', sep="")
dir.create(file.path(gopathDiagNOTExpME))
write.table(ov.EE.ME.DiagNOTExpME, file=file.path(gopathDiagNOTExpME, "GO_pvalues.txt"), quote=F, row.names=F, sep="\t")



#### MM.EM
#### subset to only non na data
ov.MM.EM.all <- data.frame(merged5$gene, merged5$logFC.maleEffectOfM_MM.EM, merged5$FDR.maleEffectOfM_MM.EM, merged5$logFC.matingE_E.EEInv, merged5$logFC.matingM_M.MMInv)
names(ov.MM.EM.all) <- c('gene', 'logFC.maleEffectOfM_MM.EM', 'FDR.maleEffectOfM_MM.EM', 'logFC.matingE_E.EEInv', 'logFC.matingM_M.MMInv')
ov.MM.EM.all <- subset(ov.MM.EM.all, !is.na(ov.MM.EM.all$FDR.maleEffectOfM_MM.EM) & !is.na(ov.MM.EM.all$logFC.matingE_E.EEInv) & !is.na(ov.MM.EM.all$logFC.matingM_M.MMInv))

##### keep only the significant values on diagonal
ov.MM.EM.all$FDR.MM.EM.DiagExpEE <- ov.MM.EM.all$FDR.maleEffectOfM_MM.EM
ov.MM.EM.all$FDR.MM.EM.DiagExpME <- ov.MM.EM.all$FDR.maleEffectOfM_MM.EM
# points(ov.MM.EM.all$logFC.matingM_M.MMInv[ov.MM.EM.all$FDR.maleEffectOfM_MM.EM < 0.05], ov.MM.EM.all$logFC.matingE_E.EEInv[ov.MM.EM.all$FDR.maleEffectOfM_MM.EM < 0.05], pch=5, col=1, cex=1.3)
# points(ov.MM.EM.all$logFC.matingM_M.MMInv[ov.MM.EM.all$FDR.MM.EM.DiagExp < 0.05], ov.MM.EM.all$logFC.matingE_E.EEInv[ov.MM.EM.all$FDR.MM.EM.DiagExp < 0.05], pch=5, col=2, cex=1.3)

ov.MM.EM.all$FDR.MM.EM.DiagExpEE[(ov.MM.EM.all$logFC.maleEffectOfM_MM.EM > 0 & ov.MM.EM.all$logFC.matingE_E.EEInv > 0 & ov.MM.EM.all$logFC.matingM_M.MMInv < 0 ) | (ov.MM.EM.all$logFC.maleEffectOfM_MM.EM > 0 & ov.MM.EM.all$logFC.matingE_E.EEInv < 0 & ov.MM.EM.all$logFC.matingM_M.MMInv > 0)] <- 1 

##### Make diagonal data frame for GO output
ov.MM.EM.DiagExpEE <- data.frame(ov.MM.EM.all$gene, ov.MM.EM.all$FDR.MM.EM.DiagExpEE)
names(ov.MM.EM.DiagExpEE) <- c('gene', 'FDR.MM.EM.DiagExp')

gopathDiagExpEE <- paste(outpath, 'ovMMEMDiagExpEE/', sep="")
dir.create(file.path(gopathDiagExpEE))
write.table(ov.MM.EM.DiagExpEE, file=file.path(gopathDiagExpEE, "GO_pvalues.txt"), quote=F, row.names=F, sep="\t")

ov.MM.EM.all$FDR.MM.EM.DiagExpME[(ov.MM.EM.all$logFC.maleEffectOfM_MM.EM < 0 & ov.MM.EM.all$logFC.matingE_E.EEInv > 0 & ov.MM.EM.all$logFC.matingM_M.MMInv < 0 ) | (ov.MM.EM.all$logFC.maleEffectOfM_MM.EM < 0 & ov.MM.EM.all$logFC.matingE_E.EEInv < 0 & ov.MM.EM.all$logFC.matingM_M.MMInv > 0)] <- 1 

##### Make diagonal data frame for GO output
ov.MM.EM.DiagExpME <- data.frame(ov.MM.EM.all$gene, ov.MM.EM.all$FDR.MM.EM.DiagExpME)
names(ov.MM.EM.DiagExpME) <- c('gene', 'FDR.MM.EM.DiagExp')

gopathDiagExpME <- paste(outpath, 'ovMMEMDiagExpME/', sep="")
dir.create(file.path(gopathDiagExpME))
write.table(ov.MM.EM.DiagExpME, file=file.path(gopathDiagExpME, "GO_pvalues.txt"), quote=F, row.names=F, sep="\t")


##### keep only the significant values on NOT exp diagonal
ov.MM.EM.all$FDR.MM.EM.DiagNOTExpEE <- ov.MM.EM.all$FDR.maleEffectOfM_MM.EM
ov.MM.EM.all$FDR.MM.EM.DiagNOTExpME <- ov.MM.EM.all$FDR.maleEffectOfM_MM.EM

ov.MM.EM.all$FDR.MM.EM.DiagNOTExpEE[(ov.MM.EM.all$logFC.maleEffectOfM_MM.EM > 0 & ov.MM.EM.all$logFC.matingE_E.EEInv > 0 & ov.MM.EM.all$logFC.matingM_M.MMInv > 0 ) | (ov.MM.EM.all$logFC.maleEffectOfM_MM.EM > 0 & ov.MM.EM.all$logFC.matingE_E.EEInv < 0 & ov.MM.EM.all$logFC.matingM_M.MMInv < 0)] <- 1 

##### Make diagonal data frame for GO output
ov.MM.EM.DiagNOTExpEE <- data.frame(ov.MM.EM.all$gene, ov.MM.EM.all$FDR.MM.EM.DiagNOTExpEE)
names(ov.MM.EM.DiagNOTExpEE) <- c('gene', 'FDR.MM.EM.DiagNOTExp')

gopathDiagNOTExpEE <- paste(outpath, 'ovMMEMDiagNOTExpEE/', sep="")
dir.create(file.path(gopathDiagNOTExpEE))
write.table(ov.MM.EM.DiagNOTExpEE, file=file.path(gopathDiagNOTExpEE, "GO_pvalues.txt"), quote=F, row.names=F, sep="\t")

ov.MM.EM.all$FDR.MM.EM.DiagNOTExpEE[(ov.MM.EM.all$logFC.maleEffectOfM_MM.EM < 0 & ov.MM.EM.all$logFC.matingE_E.EEInv > 0 & ov.MM.EM.all$logFC.matingM_M.MMInv > 0 ) | (ov.MM.EM.all$logFC.maleEffectOfM_MM.EM < 0 & ov.MM.EM.all$logFC.matingE_E.EEInv < 0 & ov.MM.EM.all$logFC.matingM_M.MMInv < 0)] <- 1 


##### Make diagonal data frame for GO output
ov.MM.EM.DiagNOTExpME <- data.frame(ov.MM.EM.all$gene, ov.MM.EM.all$FDR.MM.EM.DiagNOTExpME)
names(ov.MM.EM.DiagNOTExpME) <- c('gene', 'FDR.MM.EM.DiagNOTExp')

gopathDiagNOTExpME <- paste(outpath, 'ovMMEMDiagNOTExpME/', sep="")
dir.create(file.path(gopathDiagNOTExpME))
write.table(ov.MM.EM.DiagNOTExpME, file=file.path(gopathDiagNOTExpME, "GO_pvalues.txt"), quote=F, row.names=F, sep="\t")
