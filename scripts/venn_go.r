# mating FRT

matingE <- read.table('~/git/dpse_mating/input/venn/mating/GO_rt_venn_mating/FRT_GO_pvalues_matingE.txt', header=T)

matingM <- read.table('~/git/dpse_mating/input/venn/mating/GO_rt_venn_mating/FRT_GO_pvalues_matingM.txt', header=T)

upE <- read.table('~/git/dpse_mating/input/venn/mating/GO_rt_venn_mating/upE.txt', header=T)
upE$FDR_upE <- 0.001
upEE <- read.table('~/git/dpse_mating/input/venn/mating/GO_rt_venn_mating/upEE.txt', header=T)
upEE$FDR_upEE <- 0.001
upM_upE <- read.table('~/git/dpse_mating/input/venn/mating/GO_rt_venn_mating/upM_upE.txt', header=T)
upM_upE$FDR_upM_upE <- 0.001
upM <- read.table('~/git/dpse_mating/input/venn/mating/GO_rt_venn_mating/upM.txt', header=T)
upM$FDR_upM <- 0.001
upMM_upE <- read.table('~/git/dpse_mating/input/venn/mating/GO_rt_venn_mating/upMM_upE.txt', header=T)
upMM_upE$FDR_upMM_upE <- 0.001
upMM_upEE <- read.table('~/git/dpse_mating/input/venn/mating/GO_rt_venn_mating/upMM_upEE.txt', header=T)
upMM_upEE$FDR_upMM_upEE <- 0.001
upMM <- read.table('~/git/dpse_mating/input/venn/mating/GO_rt_venn_mating/upMM.txt', header=T)
upMM$FDR_upMM <- 0.001

merged1 <- merge(matingE, matingM, by.x='gene', by.y='gene', all=T)
merged2 <- merge(merged1, upE, by.x='gene', by.y='gene', all=T)
merged3 <- merge(merged2, upEE, by.x='gene', by.y='gene', all=T)
merged4 <- merge(merged3, upM_upE, by.x='gene', by.y='gene', all=T)
merged5 <- merge(merged4, upM, by.x='gene', by.y='gene', all=T)
merged6 <- merge(merged5, upMM_upE, by.x='gene', by.y='gene', all=T)
merged7 <- merge(merged6, upMM_upEE, by.x='gene', by.y='gene', all=T)
merged8 <- merge(merged7, upMM, by.x='gene', by.y='gene', all=T)

merged8$FDR_upE[is.na(merged8$FDR_upE)] <- 1
merged8$FDR_upEE[is.na(merged8$FDR_upEE)] <- 1
merged8$FDR_upM[is.na(merged8$FDR_upM)] <- 1
merged8$FDR_upM_upE[is.na(merged8$FDR_upM_upE)] <- 1
merged8$FDR_upMM_upE[is.na(merged8$FDR_upMM_upE)] <- 1
merged8$FDR_upMM_upEE[is.na(merged8$FDR_upMM_upEE)] <- 1
merged8$FDR_upMM[is.na(merged8$FDR_upMM)] <- 1

outpath_upE <- '~/git/dpse_mating/output/venn/GO/mating/upE/'
dir.create(outpath_upE, recursive=T)
write.table(data.frame(merged8$gene, merged8$FDR_upE), file=file.path(outpath_upE, "GO_pvalues.txt"), quote=F, row.names=F, sep="\t")

outpath_upEE <- '~/git/dpse_mating/output/venn/GO/mating/upEE'
outpath_upM <- '~/git/dpse_mating/output/venn/GO/mating/upM'
outpath_upM_upE <- '~/git/dpse_mating/output/venn/GO/mating/upM_upE'
outpath_upMM_upE <- '~/git/dpse_mating/output/venn/GO/mating/upMM_upE'
outpath_upMM_upEE <- '~/git/dpse_mating/output/venn/GO/mating/upMM_upEE'
outpath_upMM <- '~/git/dpse_mating/output/venn/GO/mating/upMM'

dir.create(outpath_upEE)
dir.create(outpath_upM)
dir.create(outpath_upM_upE)
dir.create(outpath_upMM_upE)
dir.create(outpath_upMM_upEE)
dir.create(outpath_upMM)

write.table(data.frame(merged8$gene, merged8$FDR_upEE), file=file.path(outpath_upEE, "GO_pvalues.txt"), quote=F, row.names=F, sep="\t")
write.table(data.frame(merged8$gene, merged8$FDR_upM), file=file.path(outpath_upM, "GO_pvalues.txt"), quote=F, row.names=F, sep="\t")
write.table(data.frame(merged8$gene, merged8$FDR_upM_upE), file=file.path(outpath_upM_upE, "GO_pvalues.txt"), quote=F, row.names=F, sep="\t")
write.table(data.frame(merged8$gene, merged8$FDR_upMM_upE), file=file.path(outpath_upMM_upE, "GO_pvalues.txt"), quote=F, row.names=F, sep="\t")
write.table(data.frame(merged8$gene, merged8$FDR_upMM_upEE), file=file.path(outpath_upMM_upEE, "GO_pvalues.txt"), quote=F, row.names=F, sep="\t")
write.table(data.frame(merged8$gene, merged8$FDR_upMM), file=file.path(outpath_upMM, "GO_pvalues.txt"), quote=F, row.names=F, sep="\t")


# of FRT

ofE <- read.table('~/git/dpse_mating/input/venn/FRT_of/GO_rt_venn_of/FRT_GO_pvalues_ofE.txt', header=T)

ofM <- read.table('~/git/dpse_mating/input/venn/FRT_of/GO_rt_venn_of/FRT_GO_pvalues_ofM.txt', header=T)

onE_upEE_onM_upEM <- read.table('~/git/dpse_mating/input/venn/FRT_of/GO_rt_venn_of/onE_upEE_onM_upEM.txt', header=T)
onE_upEE_onM_upEM$FDR_onE_upEE_onM_upEM <- 0.001
onE_upEE <- read.table('~/git/dpse_mating/input/venn/FRT_of/GO_rt_venn_of/onE_upEE.txt', header=T)
onE_upEE$FDR_onE_upEE <- 0.001
onE_upEM <- read.table('~/git/dpse_mating/input/venn/FRT_of/GO_rt_venn_of/onE_upEM.txt', header=T)
onE_upEM$FDR_onE_upEM <- 0.001
onE_upEM_onM_upMM <- read.table('~/git/dpse_mating/input/venn/FRT_of/GO_rt_venn_of/onE_upEM_onM_upMM.txt', header=T)
onE_upEM_onM_upMM$FDR_onE_upEM_onM_upMM <- 0.001
onM_upMM <- read.table('~/git/dpse_mating/input/venn/FRT_of/GO_rt_venn_of/onM_upMM.txt', header=T)
onM_upMM$FDR_onM_upMM <- 0.001
onM_upEM <- read.table('~/git/dpse_mating/input/venn/FRT_of/GO_rt_venn_of/onM_upEM.txt', header=T)
onM_upEM$FDR_onM_upEM <- 0.001

merged1 <- merge(ofE, ofM, by.x='gene', by.y='gene', all=T)
merged2 <- merge(merged1, onE_upEE_onM_upEM, by.x='gene', by.y='gene', all=T)
merged3 <- merge(merged2, onE_upEE, by.x='gene', by.y='gene', all=T)
merged4 <- merge(merged3, onE_upEM, by.x='gene', by.y='gene', all=T)
merged5 <- merge(merged4, onE_upEM_onM_upMM, by.x='gene', by.y='gene', all=T)
merged6 <- merge(merged5, onM_upMM, by.x='gene', by.y='gene', all=T)
merged7 <- merge(merged6, onM_upEM, by.x='gene', by.y='gene', all=T)
merged8<- merged7

merged8$FDR_onE_upEE_onM_upEM[is.na(merged8$FDR_onE_upEE_onM_upEM)] <- 1
merged8$FDR_onE_upEE[is.na(merged8$FDR_onE_upEE)] <- 1
merged8$FDR_onE_upEM[is.na(merged8$FDR_onE_upEM)] <- 1
merged8$FDR_onE_upEM_onM_upMM[is.na(merged8$FDR_onE_upEM_onM_upMM)] <- 1
merged8$FDR_onM_upMM[is.na(merged8$FDR_onM_upMM)] <- 1
merged8$FDR_onM_upEM[is.na(merged8$FDR_onM_upEM)] <- 1

outpath_onE_upEE_onM_upEM <- '~/git/dpse_mating/output/venn/GO/FRT_of/onE_upEE_onM_upEM'
outpath_onE_upEE <- '~/git/dpse_mating/output/venn/GO/FRT_of/onE_upEE'
outpath_onE_upEM <- '~/git/dpse_mating/output/venn/GO/FRT_of/onE_upEM'
outpath_onE_upEM_onM_upMM <- '~/git/dpse_mating/output/venn/GO/FRT_of/onE_upEM_onM_upMM'
outpath_onM_upMM <- '~/git/dpse_mating/output/venn/GO/FRT_of/onM_upMM'
outpath_onM_upEM <- '~/git/dpse_mating/output/venn/GO/FRT_of/onM_upEM'

dir.create(outpath_onE_upEE_onM_upEM, recursive=T)
dir.create(outpath_onE_upEE)
dir.create(outpath_onE_upEM)
dir.create(outpath_onE_upEM_onM_upMM)
dir.create(outpath_onM_upMM)
dir.create(outpath_onM_upEM)

write.table(data.frame(merged8$gene, merged8$FDR_onE_upEE_onM_upEM), file=file.path(outpath_onE_upEE_onM_upEM, "GO_pvalues.txt"), quote=F, row.names=F, sep="\t")
write.table(data.frame(merged8$gene, merged8$FDR_onE_upEE), file=file.path(outpath_onE_upEE, "GO_pvalues.txt"), quote=F, row.names=F, sep="\t")
write.table(data.frame(merged8$gene, merged8$FDR_onE_upEM), file=file.path(outpath_onE_upEM, "GO_pvalues.txt"), quote=F, row.names=F, sep="\t")
write.table(data.frame(merged8$gene, merged8$FDR_onE_upEM_onM_upMM), file=file.path(outpath_onE_upEM_onM_upMM, "GO_pvalues.txt"), quote=F, row.names=F, sep="\t")
write.table(data.frame(merged8$gene, merged8$FDR_onM_upMM), file=file.path(outpath_onM_upMM, "GO_pvalues.txt"), quote=F, row.names=F, sep="\t")
write.table(data.frame(merged8$gene, merged8$FDR_onM_upEM), file=file.path(outpath_onM_upEM, "GO_pvalues.txt"), quote=F, row.names=F, sep="\t")

# ov of
					/git/dpse_mating/input/venn/OV_of/GO_ov_venn_of/
ofE <- read.table('~/git/dpse_mating/input/venn/OV_of/GO_ov_venn_of/ov_GO_pvalues_ofE.txt', header=T)

ofM <- read.table('~/git/dpse_mating/input/venn/OV_of/GO_ov_venn_of/ov_GO_pvalues_ofM.txt', header=T)

onE_upEE_onM_upEM <- read.table('~/git/dpse_mating/input/venn/OV_of/GO_ov_venn_of/onE_upEE_onM_upEM.txt', header=T)
onE_upEE_onM_upEM$FDR_onE_upEE_onM_upEM <- 0.001
onE_upEE <- read.table('~/git/dpse_mating/input/venn/OV_of/GO_ov_venn_of/onE_upEE.txt', header=T)
onE_upEE$FDR_onE_upEE <- 0.001
onE_upEM <- read.table('~/git/dpse_mating/input/venn/OV_of/GO_ov_venn_of/onE_upEM.txt', header=T)
onE_upEM$FDR_onE_upEM <- 0.001
onE_upEM_onM_upMM <- read.table('~/git/dpse_mating/input/venn/OV_of/GO_ov_venn_of/onE_upEM_onM_upMM.txt', header=T)
onE_upEM_onM_upMM$FDR_onE_upEM_onM_upMM <- 0.001
onM_upMM <- read.table('~/git/dpse_mating/input/venn/OV_of/GO_ov_venn_of/onM_upMM.txt', header=T)
onM_upMM$FDR_onM_upMM <- 0.001
onM_upEM <- read.table('~/git/dpse_mating/input/venn/OV_of/GO_ov_venn_of/onM_upEM.txt', header=T)
onM_upEM$FDR_onM_upEM <- 0.001

merged1 <- merge(ofE, ofM, by.x='gene', by.y='gene', all=T)
merged2 <- merge(merged1, onE_upEE_onM_upEM, by.x='gene', by.y='gene', all=T)
merged3 <- merge(merged2, onE_upEE, by.x='gene', by.y='gene', all=T)
merged4 <- merge(merged3, onE_upEM, by.x='gene', by.y='gene', all=T)
merged5 <- merge(merged4, onE_upEM_onM_upMM, by.x='gene', by.y='gene', all=T)
merged6 <- merge(merged5, onM_upMM, by.x='gene', by.y='gene', all=T)
merged7 <- merge(merged6, onM_upEM, by.x='gene', by.y='gene', all=T)
merged8 <- merged7

merged8$FDR_onE_upEE_onM_upEM[is.na(merged8$FDR_onE_upEE_onM_upEM)] <- 1
merged8$FDR_onE_upEE[is.na(merged8$FDR_onE_upEE)] <- 1
merged8$FDR_onE_upEM[is.na(merged8$FDR_onE_upEM)] <- 1
merged8$FDR_onE_upEM_onM_upMM[is.na(merged8$FDR_onE_upEM_onM_upMM)] <- 1
merged8$FDR_onM_upMM[is.na(merged8$FDR_onM_upMM)] <- 1
merged8$FDR_onM_upEM[is.na(merged8$FDR_onM_upEM)] <- 1
	
	
outpath_onE_upEE_onM_upEM <- '~/git/dpse_mating/output/venn/GO/OV_of/onE_upEE_onM_upEM'
outpath_onE_upEE <- '~/git/dpse_mating/output/venn/GO/OV_of/onE_upEE'
outpath_onE_upEM <- '~/git/dpse_mating/output/venn/GO/OV_of/onE_upEM'
outpath_onE_upEM_onM_upMM <- '~/git/dpse_mating/output/venn/GO/OV_of/onE_upEM_onM_upMM'
outpath_onM_upMM <- '~/git/dpse_mating/output/venn/GO/OV_of/onM_upMM'
outpath_onM_upEM <- '~/git/dpse_mating/output/venn/GO/OV_of/onM_upEM'

dir.create(outpath_onE_upEE_onM_upEM, recursive=T)
dir.create(outpath_onE_upEE)
dir.create(outpath_onE_upEM)
dir.create(outpath_onE_upEM_onM_upMM)
dir.create(outpath_onM_upMM)
dir.create(outpath_onM_upEM)

write.table(data.frame(merged8$gene, merged8$FDR_onE_upEE_onM_upEM), file=file.path(outpath_onE_upEE_onM_upEM, "GO_pvalues.txt"), quote=F, row.names=F, sep="\t")
write.table(data.frame(merged8$gene, merged8$FDR_onE_upEE), file=file.path(outpath_onE_upEE, "GO_pvalues.txt"), quote=F, row.names=F, sep="\t")
write.table(data.frame(merged8$gene, merged8$FDR_onE_upEM), file=file.path(outpath_onE_upEM, "GO_pvalues.txt"), quote=F, row.names=F, sep="\t")
write.table(data.frame(merged8$gene, merged8$FDR_onE_upEM_onM_upMM), file=file.path(outpath_onE_upEM_onM_upMM, "GO_pvalues.txt"), quote=F, row.names=F, sep="\t")
write.table(data.frame(merged8$gene, merged8$FDR_onM_upMM), file=file.path(outpath_onM_upMM, "GO_pvalues.txt"), quote=F, row.names=F, sep="\t")
write.table(data.frame(merged8$gene, merged8$FDR_onM_upEM), file=file.path(outpath_onM_upEM, "GO_pvalues.txt"), quote=F, row.names=F, sep="\t")
