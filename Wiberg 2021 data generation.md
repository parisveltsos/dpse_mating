The contrasts providing DE genes for Wiberg 2021 are

virgin E vs M per tissue

	Rscript ~/git/dpse_mating/scripts/de_analysis.r testis 0.05
	Rscript ~/git/dpse_mating/scripts/de_analysis.r agland 0.05
	Rscript ~/git/dpse_mating/scripts/de_analysis.r rtract_SSvirgins 0.05
	Rscript ~/git/dpse_mating/scripts/de_analysis.r ovary_SSvirgins 0.05

mated E vs M (within line) for female tissues

	Rscript ~/git/dpse_mating/scripts/de_analysis.r rtract_SSmated 0.05
	Rscript ~/git/dpse_mating/scripts/de_analysis.r ovary_SSmated_rmEEov4 0.05

These DE genes were merged, along with genes from virgin and courted heads and abdomens of each sex from https://doi.org/10.5281/zenodo.1044632 in file `/Users/pveltsos/git/dpse_mating/input/Wiberg_DE_genes.txt`.
 