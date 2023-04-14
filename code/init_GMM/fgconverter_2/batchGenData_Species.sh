
#INPUTDIR=/ahg/regev/users/sroy/compfuncgen/data/anajay/OS_KCl-DE/species_specific
#OUTPUTDIR=/ahg/regev/users/sroy/compfuncgen/aurianproj/fgfiles/
OUTPUTDIR=~/data/compfuncgen/heat/species_specific
INPUTDIR=/work/sroy/regev_lab/regev/sroy/compfuncgen/data/stressdata/ilandata/ALL/species_specific/
#OUTPUTDIR=/ahg/regev/users/sroy/compfuncgen/analysis/fgfiles/stressdata/HeOsOx_rps/
#OUTPUTDIR=/ahg/regev/users/sroy/compfuncgen/analysis/fgfiles/stressdata/OSstress_scersbay_highlydiv/
mkdir -p $OUTPUTDIR
FGCONVTOOL=./genData
#for cutoff in 0.24 0.248
	#SUBDIR=$OUTPUTDIR/c${cutoff}
	#SUBDIR=$OUTPUTDIR/
	#mkdir -p $SUBDIR
#for SPECIES in scer sbay cgla calb spom ylip klac scas
#for SPECIES in scer sbay cgla calb spom ylip klac scas
for SPECIES in Scer Scas Cgla Calb Klac Kwal Spom Sjap
do
	SPECIESLOWER=`echo $SPECIES | awk '{printf("%s",tolower($1))}'`
	#EXPRFILE=$INPUTDIR/$SPECIES/high2_std2_merged.genexp
	if [ $SPECIES ==  Cgla ]
	then
	EXPRFILE=$INPUTDIR/$SPECIESLOWER/nointerpolate_heat37_0.50miss.geneexp
	else if [ $SPECIES == Calb ]
	then
	EXPRFILE=$INPUTDIR/$SPECIESLOWER/nointerpolate_heat42_0.50miss.geneexp
	else
	EXPRFILE=$INPUTDIR/$SPECIESLOWER/nointerpolate_heat_0.50miss.geneexp
	fi
	fi

	MEASUREMENTS=`awk '{print NF-1}' $EXPRFILE| sort -u`
	MODELSUFF=$OUTPUTDIR/$SPECIES
	echo "$FGCONVTOOL $EXPRFILE 1 $MODELSUFF $MEASUREMENTS fgraph 1 none"
	$FGCONVTOOL $EXPRFILE 1 $MODELSUFF $MEASUREMENTS fgraph 1 none
done
