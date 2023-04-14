
INPUTDIR=/ahg/regev/users/sroy/compfuncgen/data/anajay/OS_KCl-DE/species_specific
INPUTDIR=/ahg/regev/users/sroy/compfuncgen/data/stressdata/ilandata/ALL/species_specific/
#OUTPUTDIR=/ahg/regev/users/sroy/compfuncgen/aurianproj/fgfiles/
OUTPUTDIR=/ahg/regev/users/sroy/compfuncgen/analysis/fgfiles/stressdata/ilandata
CLUSTERDIR=/seq/compbio/sroy/networklearning/fgraph/results/gmm_crosspecies/stressdata/ilan_heatstress_8spec_pp_postem/crossspecies_k5/randinit1/lf0.8_nlf0.8
mkdir -p $OUTPUTDIR
FGCONVTOOL=./genData
#for cutoff in 0.24 0.248
#SUBDIR=$OUTPUTDIR/c${cutoff}
SUBDIR=$OUTPUTDIR/speciesmodelfiles/
mkdir $SUBDIR
for SPECIES in Scer Cgla Scas Calb Kwal Klac Sjap Spom
do
	export SPECIES
	SPECIESLOWER=`echo $SPECIES | awk '{printf("%s\n",tolower(ENVIRON["SPECIES"]))}'`
	#EXPRFILE=$INPUTDIR/$SPECIES/high2_std2_merged.genexp
	if [ $SPECIES == "Cgla" -o $SPECIES == "Calb" ]
	then
		TEMPFILE=$INPUTDIR/$SPECIESLOWER/nointerpolate_heat42_0.50miss.geneexp
	else
		TEMPFILE=$INPUTDIR/$SPECIESLOWER/nointerpolate_heat_0.50miss.geneexp
	fi
	EXPRFILE=temp.txt
	awk '{if (NR>=1) print $0}' $TEMPFILE> $EXPRFILE
	MODELSUFF=$SUBDIR/$SPECIES
	TIMEPNT=`head -n1 $EXPRFILE| awk '{print NF-1}'` 
	echo "$FGCONVTOOL $EXPRFILE 1 $MODELSUFF $TIMEPNT fgraph 1"
	$FGCONVTOOL $EXPRFILE 1 $MODELSUFF $TIMEPNT fgraph 1
done
