# Installation

Install the following libraries are available

	install.packages("gplots", dependencies=TRUE)

	install.packages("ggplot2", dependencies=TRUE)

	install.packages("sjPlot", dependencies=TRUE)

	install.packages("dynamicTreeCut", dependencies=TRUE)

Install topGO and edgeR from bioconductor

## Local installation

	mkdir ~/git

	cd ~/git

	git clone https://GitHub.com/parisveltsos/dpse_mating.git

	cd ~/git/dpse_mating

	mkdir output

# Analyses (in order)

## Main DE analysis

Note, more output than used in the paper is produced. The script 

	* performs DE analysis
	* plots MDS 
	* plots volcano plot
	* performs and plots chisq tests of DE gene number comparisons between different chromosomes, and produces the relevant tables with sjPlot
	* plots the chromosomal distribution of significant genes
	* makes violin plots showing differential expression magnitude in different chromosomes
	* makes heatmaps of the DE genes using all libraries
	
Run the EM contrast analyses for all tissues, with 5% FDR.

	source ~/git/dpse_mating/scripts/batch_analysis.sh

The input files have removed the rDNA genes, and library EEov4 which was an MDS plot outlier.

## GO analysis

Navigate to the folder of the relevant comparisons, and run the analysis. This needs to be done in the terminal for each contrast. eg for testis

	cd ~/git/dpse_mating/output/testis/GO_0.05_maleVirgins_E.M

	source ~/git/dpse_mating/scripts/GO_analysis_05.sh

The `Fisher.txt` files in the GO_0.05... and GO_DOWN and GO_UP folders make up the GO analysis results ("All DE" "Up" or "Down") summarized in File S1. These refer to the legend of each panel in Fig 3.

File S1 also shows GO terms for congruent and non-congruent genes from panel a. These were identified by constructing venn diagrams.

### Venn diagrams

	mkdir -p ~/git/dpse_mating/output/venn/GO/

DE genes from the analysis, eg `~/git/dpse_mating/output/dedata/ovary_matingE_rmEEov4/de_0.05_0_ovary_matingE_rmEEov4.txt` were used as input in [venny](https://bioinfogp.cnb.csic.es/tools/venny/) to generate the datasets in `/Users/pveltsos/git/dpse_mating/input/venn/`

	source ~/git/dpse_mating/scripts/venn_go.r
	
The script makes the GO universe of all genes in each venn diagram.

Conduct the GO analysis (only virgin/mated reported)

	cd ~/git/dpse_mating/output/venn/GO/mating/

	for i in $(ls); do cd $i; source ~/git/dpse_mating/scripts/GO_venn_analysis_05.sh; cd ..; done

	cd ~/git/dpse_mating/output/venn/GO/FRT_of/

	for i in $(ls); do cd $i; source ~/git/dpse_mating/scripts/GO_venn_analysis_05.sh; cd ..; done

	cd ~/git/dpse_mating/output/venn/GO/OV_of/

	for i in $(ls); do cd $i; source ~/git/dpse_mating/scripts/GO_venn_analysis_05.sh; cd ..; done

The results are written in `~/git/dpse_mating/output/venn/GO`

## Scatterplot

Fig 3 showing congruent and non congruent genes is made with

	source ~/git/dpse_mating/scripts/scatterPlot.r

The script also calculates the number of genes with different symbols in the plot, and produces files for GO analysis. This analysis is not presented.


### GO analysis of Fig 3 shapes

This analysis is not presented.

The subsets are produced by `~/git/dpse_mating/scripts/scatterPlot.r`

	cd /Users/pveltsos/git/dpse_mating/output/scatterplot/

## Agcp analysis

This analysis produces figure 2, comparing DE genes from this study with agcp genes from Karr 2019 https://doi.org/10.1074/mcp.RA118.001098. 

Run `~/git/dpse_mating/scripts/agcpPlot.r` See script for origin of input files.

## Ridge plots

The data in `input/ridgePlots` were made from the output of the DE analysis. For example the `input/ridgePlots/mating.txt` data whose second column reads `ov_E_mated` or `ov_E_virgin` show the logFC value of the DE genes in the ovary mated vs virgin contrast from file `~/git/dpse_mating/output/dedata/ovary_matingE_rmEEov4/de_0.05_0_ovary_matingE_rmEEov4.txt`.

Make the plots for Fig S1, Fig S2 and Fig S3, and the associated Wilcoxon tests with 

	source ~/git/dpse_mating/scripts/ridgeplots_wilcoxon_reproseq.r

The plots were manipulated for presentation in Affinity Designer, mostly to change the location of the gene numbers, add the Wilcoxon test asterisks and change the text of the y axis.

