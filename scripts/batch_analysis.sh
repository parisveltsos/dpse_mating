# Main
## males

Rscript ~/git/dpse_mating/scripts/de_analysis.r testis 0.05
Rscript ~/git/dpse_mating/scripts/de_analysis.r agland 0.05

## females - all 6 library types used in single model, EEov4 library is an outlier in MDS plot and was removed
Rscript ~/git/dpse_mating/scripts/de_analysis.r rtract 0.05
Rscript ~/git/dpse_mating/scripts/de_analysis.r ovary_rmEEov4 0.05


# Not used - analyses of libraries of specific contrasts only
# Rscript ~/git/dpse_mating/scripts/de_analysis.r rtract_SSvirgins 0.05
# Rscript ~/git/dpse_mating/scripts/de_analysis.r ovary_SSvirgins 0.05
# 
# ## within line virgin mated contrast
# Rscript ~/git/dpse_mating/scripts/de_analysis.r ovary_matingE_rmEEov4 0.05
# Rscript ~/git/dpse_mating/scripts/de_analysis.r ovary_matingM 0.05
# Rscript ~/git/dpse_mating/scripts/de_analysis.r rtract_matingE 0.05
# Rscript ~/git/dpse_mating/scripts/de_analysis.r rtract_matingM 0.05
# 
# ## Between line effect of male
# Rscript ~/git/dpse_mating/scripts/de_analysis.r rtract_maleEffectOnE 0.05
# Rscript ~/git/dpse_mating/scripts/de_analysis.r ovary_rmEEov4_maleEffectOnE 0.05
# Rscript ~/git/dpse_mating/scripts/de_analysis.r rtract_maleEffectOnM 0.05
# Rscript ~/git/dpse_mating/scripts/de_analysis.r ovary_maleEffectOnM 0.05
# 
# Rscript ~/git/dpse_mating/scripts/de_analysis.r rtract_maleEffectOfE 0.05
# Rscript ~/git/dpse_mating/scripts/de_analysis.r ovary_rmEEov4_maleEffectOfE 0.05
# Rscript ~/git/dpse_mating/scripts/de_analysis.r rtract_maleEffectOfM 0.05
# Rscript ~/git/dpse_mating/scripts/de_analysis.r ovary_maleEffectOfM 0.05




# Compare E and M mated - not shown
## The first two DE genes were used by Wiberg 2021
Rscript ~/git/dpse_mating/scripts/de_analysis.r rtract_SSmated 0.05
Rscript ~/git/dpse_mating/scripts/de_analysis.r ovary_SSmated_rmEEov4 0.05
# Rscript ~/git/dpse_mating/scripts/de_analysis.r rtract_SSmatedBothMales 0.05
# Rscript ~/git/dpse_mating/scripts/de_analysis.r ovary_SSmatedBothMales 0.05
# Rscript ~/git/dpse_mating/scripts/de_analysis.r rtract_matingWithinLine 0.05
# Rscript ~/git/dpse_mating/scripts/de_analysis.r ovary_rmEEov4_matingWithinLine 0.05
# Rscript ~/git/dpse_mating/scripts/de_analysis.r rtract_maleEffectOnBoth 0.05
# Rscript ~/git/dpse_mating/scripts/de_analysis.r ovary_rmEEov4_maleEffectOnBoth 0.05
# Rscript ~/git/dpse_mating/scripts/de_analysis.r rtract_maleEffectOfBoth 0.05
# Rscript ~/git/dpse_mating/scripts/de_analysis.r ovary_rmEEov4_maleEffectOfBoth 0.05

# effect of erroneous library, clear from MDS plot EEov4 is problematic
# Rscript ~/git/dpse_mating/scripts/de_analysis.r ovary 0.05
# Rscript ~/git/dpse_mating/scripts/de_analysis.r ovary_SSmated 0.05
# Rscript ~/git/dpse_mating/scripts/de_analysis.r ovary_matingE 0.05
# Rscript ~/git/dpse_mating/scripts/de_analysis.r ovary_matingWithinLine 0.05
# Rscript ~/git/dpse_mating/scripts/de_analysis.r ovary_maleEffectOnE 0.05
# Rscript ~/git/dpse_mating/scripts/de_analysis.r ovary_maleEffectOfBoth 0.05
# Rscript ~/git/dpse_mating/scripts/de_analysis.r ovary_maleEffectOnBoth 0.05
# Rscript ~/git/dpse_mating/scripts/de_analysis.r ovary_maleEffectOfE 0.05

