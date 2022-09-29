HMAWK=/mnt/ws/sysbio/roygroup/shared/scripts/figscripts/Heatmap.awk

#obtain enrichments for example data 
export GOEXE=/mnt/ws/sysbio/roygroup/shared/programs/enrichanalyzer_Nongraph_Qval/enrichAnalyzer
export GODAT=/mnt/ws/sysbio/roygroup/shared/projects/EpigenomeComparison_5Species/data/go_data/ENSEMBL/go_inputs/homo_sapiens_
export GOBG=/mnt/ws/sysbio/roygroup/shared/projects/EpigenomeComparison_5Species/data/go_data/ENSEMBL/go_inputs/homo_sapiens_bg.txt
export KEGG=/mnt/ws/sysbio/roygroup/shared/projects/EpigenomeComparison_5Species/data/go_data/KEGG/KEGG_homo_sapiens_
export KEGGBG=/mnt/ws/sysbio/roygroup/shared/projects/EpigenomeComparison_5Species/data/go_data/KEGG/KEGG_homo_sapiens_bg.txt
#eval ${GOEXE} example_out/all_genesets.txt ${GOBG} ${GODAT} 1 example_out/GeneGroups_GO persg
#eval ${GOEXE} example_out/all_genesets.txt ${KEGGBG} ${KEGG} 1 example_out/GeneGroups_KEGG persg
export CISBPM=/mnt/ws/sysbio/roygroup/shared/projects/EpigenomeComparison_5Species/data/motifs/go_inputs/CisBP_Mammals_homo_sapiens_
#eval ${GOEXE} example_out/all_genesets.txt <(cut -f1 ${CISBPM}gotermap.txt | sort -u) ${CISBPM} 1  example_out/GeneGroups_CisBP_Mammals persg

#example enrichment summary plot for a single cluster
for c in `cut -f1 example_out/all_genesets.txt | grep Cluster100 `
do
	export L=${c}
	grep -F -w ${L} example_out/GeneGroups_CisBP_Mammals_details.txt | cut -f1-4 | sort -k4g | head -n 5 | awk '$4<=0.05{printf("%s\t%s\t%f|1|%.2f\n",$2,$1,-log($4)/log(10),-log($4)/log(10))}' > example_out/${L}.txt
	grep -F -w ${L} example_out/GeneGroups_GO_details.txt | cut -f1-4 | sort -k4g | head -n 5 | awk '$4<=0.05{printf("%s\t%s\t%f|2|%.2f\n",$2,$1,-log($4)/log(10),-log($4)/log(10))}' >> example_out/${L}.txt
	grep -F -w ${L} example_out/GeneGroups_KEGG_details.txt | cut -f1-4 | sort -k4g | head -n 5 | awk '$4<=0.05{printf("%s\t%s\t%f|3|%.2f\n",$2,$1,-log($4)/log(10),-log($4)/log(10))}' >> example_out/${L}.txt
	cat example_out/${L}.txt | ${HMAWK} -vStrokeC="-" -vStrokeSC="black" -vFontSize="8" -vL="CisBPMotifs GO KEGG" -vD="space" -vLegendNum="6" -vC="0:(255,255,255);5:(0,255,0) 0:(255,255,255);5:(0,0,255) 0:(255,255,255);5:(0,0,0)" > example_plots/${L}_enrichments.svg
done

#example summary plot

#prep GO information
export F=example_out/ordered_clusterset_means.txt
export GO=example_out/GeneGroups_GO_details.txt
export SGO=example_out/GO_Summary.txt
#prep motif information
export MO=example_out/GeneGroups_CisBP_Mammals_details.txt
export SM=example_out/CisBP_Mammals_Summary.txt
if [[ -e ${SGO} ]]
then
	rm ${SGO}
fi
if [[ -e ${SM} ]]
then
	rm ${SM}
fi
for C in `cat ${GO} ${MO} | awk '$4<=0.05{print}' | cut -f1 | sort -u`
do
	awk -vc="$C" '$1==c && $4<=0.05{printf("%s\t%s\t%f\n",$1,$2,-log($4)/log(10))}' ${GO} | sort -k4gr | head -n 3 >> ${SGO}
	awk -vc="$C" '$1==c && $4<=0.05{printf("%s\t%s\t%f\n",$1,$2,-log($4)/log(10))}' ${MO} | sort -k4gr | head -n 3 >> ${SM}
done
#prep ordering based on module ids in species
if [[ -e tmp_input.txt ]]
then
        rm tmp_input.txt
fi
#integrate all module assignment information available into a matrix and then re-sort the cluster ordeA
#grep -F "|1|" ${F} > tmp_input.txt  
grep -F "|1|" ${F} | sed 's/||8//g;s/||6//g' | cut -f1 -d"|" > tmp.txt
/mnt/ws/sysbio/roygroup/shared/scripts/figscripts/TT.awk -vIF=3 -vOF=m -vOFS=$'\t' tmp.txt | awk '{print}' > tmp2.txt
#now reorder based on sorting for specific species
awk 'NR>1{print}'  tmp2.txt | sort -k2,2n -k3,3n -k5,5n | cut -f1 > tmp_ordering.txt

#Prep cluster count information
paste <(cat example_out/all_genesets.txt | awk '{print $1}') <(cat example_out/all_genesets.txt | awk '{print $2}' | sed 's/[^#]//g' | awk '{print length+1}') > tmp_count.txt
#include module assignment heatmap information
for C in `cat tmp_ordering.txt | awk '{print $1"|"}' | sed 's/Cluster/clusterset/g'`
do
	grep -F ${C} ${F}| grep -F -v Exp | grep -F -v Spacer | sed 's/Anc3/Anc1/g;s/Anc4/Anc2/g' >> tmp_input.txt
done
#first vertical spacer
printf "|- myCSP1|Spacer\n" >> tmp_input.txt
#put gene found information
awk '{printf("%s||8\tGeneCount||20\t%i|2|%i\n",$1,$2,$2)}' tmp_count.txt | sed 's/Cluster/clusterset/g' >> tmp_input.txt
printf "|- myCSP2|Spacer\n" >> tmp_input.txt
#append species expression data columns, each separated by a vertical spacer between the heatmap columns
grep -F "homo_sapiens_Exp" ${F} | sed 's/|2/|3/g' >> tmp_input.txt
printf "|- myCSP3|Spacer\n" >> tmp_input.txt
grep -F "rattus_norvegicus_Exp" ${F} | sed 's/|2/|3/g' >> tmp_input.txt
printf "|- myCSP4|Spacer\n" >> tmp_input.txt
grep -F "bos_taurus_Exp" ${F} | sed 's/|2/|3/g' >> tmp_input.txt
printf "|- myCSP5|Spacer\n" >> tmp_input.txt
#add enrichment term summaries
for C in `cut -f1 tmp_ordering.txt | sed 's/clusterset/Cluster/g'`
do
	grep -F -w ${C} ${SGO} | awk '{printf("%s||8\t%s||8\t%f|4|%0.1f\n",$1,$2,$3,$3)}' | sed 's/Cluster/clusterset/g' >> tmp_input.txt
done
printf "|- myCSPY|Spacer\n"  >> tmp_input.txt
for C in `cut -f1 tmp_ordering.txt | sed 's/clusterset/Cluster/g'`
do
	grep -F -w ${C} ${SM} | awk '{printf("%s||8\t%s||8\t%f|5|%0.1f\n",$1,$2,$3,$3)}' | sed 's/Cluster/clusterset/g' >> tmp_input.txt
done
#define heatmap color maps
export CACMAP="0:(255,255,255);1:(0,205,0);3:(255,255,0);5:(205,0,0);6:(105,0,0)"
export ECMAP="-2:(0,255,0);0:(0,0,0);2:(255,0,0)"
export CNTMAP="10:(255,255,255);100:(0,0,255)"
export GOMAP="0:(255,255,255);5:(255,0,0)"
export MMAP="0:(255,255,255);5:(0,155,0,0)"
cat tmp_input.txt | ${HMAWK} -vStrokeC="-" -vStrokeSC="black" -vL="ClusterAssignment GeneCount Activity GOEnrichment MotifEnrichment" -vD=0 -vLegendNum="7" -vC="${CACMAP} ${CNTMAP} ${ECMAP} ${GOMAP} ${MMAP}" -vFontSize=8 > example_plots/Summary_Ordered_By_CID.svg

