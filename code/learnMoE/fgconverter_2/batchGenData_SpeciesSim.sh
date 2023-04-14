
INPUTDIR=~/data/genenetweaver/multispec_learning/pg0.4_pm0.6
FGCONVTOOL=./genData
for SPECIES in Scer Sbay Cgla Calb Klac Spom Sjap
do
	EXPRFILE=$INPUTDIR/$SPECIES/processed/geneexp.txt
	MEASUREMENTS=`awk '{print NF-1}' $EXPRFILE| sort -u`
	MODELSUFF=$INPUTDIR/$SPECIES/processed/$SPECIES
	echo "$FGCONVTOOL $EXPRFILE 1 $MODELSUFF $MEASUREMENTS fgraph 1 none"
	$FGCONVTOOL $EXPRFILE 1 $MODELSUFF $MEASUREMENTS fgraph 1 none
done
