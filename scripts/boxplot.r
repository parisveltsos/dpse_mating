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