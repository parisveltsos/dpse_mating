sort -u ~/git/dpse_mating/input/Dpse_Dmel_GO.txt > GO_sorted_temp.txt

cut -f 1 GO_pvalues.txt | sort -u > GO_p_sorted_temp.txt

join GO_sorted_temp.txt GO_p_sorted_temp.txt | perl -pe 's/ /\t/' > GO_universe.txt

Rscript ~/git/dpse_mating/scripts/GO_BP.R GO_universe.txt GO_pvalues.txt out_BP 0.1
Rscript ~/git/dpse_mating/scripts/GO_CC.R GO_universe.txt GO_pvalues.txt out_CC 0.1
Rscript ~/git/dpse_mating/scripts/GO_MF.R GO_universe.txt GO_pvalues.txt out_MF 0.1

echo -e 'GO category\tGO\tTerm\tAnnotated\tSignificant\tExpected\ttopGO (Fisher)\tclassic (Fisher)' > header_temp.txt
cat table_Fisherout_BP | cut -f 1,2,3,4,5,7,9 | grep -v topgo |  perl -pe 's/^/BP\t/' > table_Fisherout_BP_temp
cat table_Fisherout_CC | cut -f 1,2,3,4,5,7,9 | grep -v topgo |  perl -pe 's/^/CC\t/' > table_Fisherout_CC_temp
cat table_Fisherout_MF | cut -f 1,2,3,4,5,7,9 | grep -v topgo |  perl -pe 's/^/MF\t/' > table_Fisherout_MF_temp

cat header_temp.txt table_Fisherout_BP_temp table_Fisherout_CC_temp table_Fisherout_MF_temp > Fisher.txt

rm *temp*

mkdir GO_UP
cp GO_pvalues_UP.txt GO_UP
cd GO_UP
Rscript ~/git/dpse_mating/scripts/GO_BP.R ../GO_universe.txt GO_pvalues_UP.txt out_BP 0.1
Rscript ~/git/dpse_mating/scripts/GO_CC.R ../GO_universe.txt GO_pvalues_UP.txt out_CC 0.1
Rscript ~/git/dpse_mating/scripts/GO_MF.R ../GO_universe.txt GO_pvalues_UP.txt out_MF 0.1
echo -e 'GO category\tGO\tTerm\tAnnotated\tSignificant\tExpected\ttopGO (Fisher)\tclassic (Fisher)' > header_temp.txt
cat table_Fisherout_BP | cut -f 1,2,3,4,5,7,9 | grep -v topgo |  perl -pe 's/^/BP\t/' > table_Fisherout_BP_temp
cat table_Fisherout_CC | cut -f 1,2,3,4,5,7,9 | grep -v topgo |  perl -pe 's/^/CC\t/' > table_Fisherout_CC_temp
cat table_Fisherout_MF | cut -f 1,2,3,4,5,7,9 | grep -v topgo |  perl -pe 's/^/MF\t/' > table_Fisherout_MF_temp
cat header_temp.txt table_Fisherout_BP_temp table_Fisherout_CC_temp table_Fisherout_MF_temp > Fisher.txt
rm *temp*
cd ..

mkdir GO_DOWN
cp GO_pvalues_DOWN.txt GO_DOWN
cd GO_DOWN
Rscript ~/git/dpse_mating/scripts/GO_BP.R ../GO_universe.txt GO_pvalues_DOWN.txt out_BP 0.1
Rscript ~/git/dpse_mating/scripts/GO_CC.R ../GO_universe.txt GO_pvalues_DOWN.txt out_CC 0.1
Rscript ~/git/dpse_mating/scripts/GO_MF.R ../GO_universe.txt GO_pvalues_DOWN.txt out_MF 0.1
echo -e 'GO category\tGO\tTerm\tAnnotated\tSignificant\tExpected\ttopGO (Fisher)\tclassic (Fisher)' > header_temp.txt
cat table_Fisherout_BP | cut -f 1,2,3,4,5,7,9 | grep -v topgo |  perl -pe 's/^/BP\t/' > table_Fisherout_BP_temp
cat table_Fisherout_CC | cut -f 1,2,3,4,5,7,9 | grep -v topgo |  perl -pe 's/^/CC\t/' > table_Fisherout_CC_temp
cat table_Fisherout_MF | cut -f 1,2,3,4,5,7,9 | grep -v topgo |  perl -pe 's/^/MF\t/' > table_Fisherout_MF_temp
cat header_temp.txt table_Fisherout_BP_temp table_Fisherout_CC_temp table_Fisherout_MF_temp > Fisher.txt
rm *temp*
cd ..


