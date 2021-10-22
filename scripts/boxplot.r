# Input data

datapath <- '~/git/dpse_mating/input/boxplots/'
outpath <- '~/git/dpse_mating/input/boxplots/'

male_data <- read.table(file.path(datapath, 'male_virgin_boxplotData.txt'), header=T)
str(male_data) 

female_data <- read.table(file.path(datapath, 'female_virgin_boxplotData.txt'), header=T)
str(female_data) 

dev.copy(pdf,file.path(outpath,'fig1.pdf'), width=12, height=6) # Fig S1
par(mfrow=c(1,2)) 
par(mar=c(5,5,4,3), oma=c(0,0,2,0))

boxplot(abs(male_data$logfc) ~ male_data$mating , col=c("gray80", "gray50"), outline=F, xlab="Library", ylab="abs(logFC)", cex.lab=1.3)
 
# Add data points
male_data$mating <- factor(male_data$mating)
mylevels <- levels(male_data$mating)
levelProportions <- summary(male_data$mating)/nrow(male_data)

for(i in 1:length(mylevels)){
  thislevel <- mylevels[i]
  thisvalues <- abs(male_data[male_data$mating==thislevel, "logfc"])
   
  # take the x-axis indices and add a jitter, proportional to the N in each level
#  myjitter <- jitter(rep(i, length(thisvalues)), amount=levelProportions[i]/2)
  myjitter <- jitter(rep(i, length(thisvalues)), amount=.25)
  points(myjitter, thisvalues, pch=20, cex=0.6, col=rgb(0,0,0,.8)) 
}  
  
  
boxplot(abs(female_data$logfc) ~ female_data$mating , col=c("gray80", "gray50"), outline=F, xlab="Library", ylab="abs(logFC)", cex.lab=1.3)
 
# Add data points
female_data$mating <- factor(female_data$mating)
mylevels <- levels(female_data$mating)
levelProportions <- summary(female_data$mating)/nrow(female_data)

for(i in 1:length(mylevels)){
  thislevel <- mylevels[i]
  thisvalues <- abs(female_data[female_data$mating==thislevel, "logfc"])
   
  # take the x-axis indices and add a jitter, proportional to the N in each level
#  myjitter <- jitter(rep(i, length(thisvalues)), amount=levelProportions[i]/2)
  myjitter <- jitter(rep(i, length(thisvalues)), amount=.25)
  points(myjitter, thisvalues, pch=20, cex=0.6, col=rgb(0,0,0,.8))    
}

dev.off()




frt_data <- read.table(file.path(datapath, 'frt_vm_data.txt'), header=T)
ov_data <- read.table(file.path(datapath, 'ov_vm_data.txt'), header=T)

dev.copy(pdf,file.path(outpath,'fig_vm.pdf'), width=12, height=6) # Fig S1
par(mfrow=c(1,2)) 
par(mar=c(5,5,4,3), oma=c(0,0,2,0))

boxplot(abs(frt_data$logfc) ~ frt_data$mating , col=c("gray80", "gray50"), outline=F, xlab="Library", ylab="abs(logFC)", cex.lab=1.3)
 
# Add data points
frt_data$mating <- factor(frt_data$mating)
mylevels <- levels(frt_data$mating)
levelProportions <- summary(frt_data$mating)/nrow(frt_data)

for(i in 1:length(mylevels)){
  thislevel <- mylevels[i]
  thisvalues <- abs(frt_data[frt_data$mating==thislevel, "logfc"])
   
  # take the x-axis indices and add a jitter, proportional to the N in each level
#  myjitter <- jitter(rep(i, length(thisvalues)), amount=levelProportions[i]/2)
  myjitter <- jitter(rep(i, length(thisvalues)), amount=.25)
  points(myjitter, thisvalues, pch=20, cex=0.6, col=rgb(0,0,0,.8)) 
}  
  
  
boxplot(abs(ov_data$logfc) ~ ov_data$mating , col=c("gray80", "gray50"), outline=F, xlab="Library", ylab="abs(logFC)", cex.lab=1.3)
 
# Add data points
ov_data$mating <- factor(ov_data$mating)
mylevels <- levels(ov_data$mating)
levelProportions <- summary(ov_data$mating)/nrow(ov_data)

for(i in 1:length(mylevels)){
  thislevel <- mylevels[i]
  thisvalues <- abs(ov_data[ov_data$mating==thislevel, "logfc"])
   
  # take the x-axis indices and add a jitter, proportional to the N in each level
#  myjitter <- jitter(rep(i, length(thisvalues)), amount=levelProportions[i]/2)
  myjitter <- jitter(rep(i, length(thisvalues)), amount=.25)
  points(myjitter, thisvalues, pch=20, cex=0.6, col=rgb(0,0,0,.8))    
}

dev.off()


wilcox.test(abs(frt_data$logfc[frt_data$mating=="frt_E.EEdown"]), frt_data$logfc[frt_data$mating=="frt_E.EEup"])
# W = 13085, p-value = 1.939e-07 FRT E.EE
wilcox.test(abs(frt_data$logfc[frt_data$mating=="frt_M.MMdown"]), frt_data$logfc[frt_data$mating=="frt_M.MMup"])
# W = 21188, p-value = 3.833e-13 FRT M.MM
wilcox.test(abs(ov_data$logfc[ov_data$mating=="ov_M.MMdown"]), ov_data$logfc[ov_data$mating=="ov_M.MMup"])
# W = 0, p-value = 1 OV M.MM
wilcox.test(abs(ov_data$logfc[ov_data$mating=="ov_E.EEdown"]), ov_data$logfc[ov_data$mating=="ov_E.EEup"])
# W = 11, p-value = 0.001176 OV E.EE



wilcox.test(abs(male_data$logfc[male_data$mating=="te_M"]), male_data$logfc[male_data$mating=="te_E"])