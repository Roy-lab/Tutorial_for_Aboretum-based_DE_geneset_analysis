
INPUTDIR=/ahg/regev/users/sroy/compfuncgen/data/anajay/PLEASE_CLUSTER_K3_K4_K5/
OUTPUTDIR=/ahg/regev/users/sroy/compfuncgen/analysis/fgfiles/anajay/PLEASE_CLUSTER_K3_K4_K5
mkdir -p $OUTPUTDIR
FGCONVTOOL=./genData
for SPECIES in KCl-Spom-diffexp KCl-Spom KCl-scer_hog1_2.0 KCl-spom_hog1_2.0
do
	TEMPFILE=$INPUTDIR/${SPECIES}_interpolate_0.50miss.geneexp
	EXPRFILE=temp.txt
	awk '{if (NR>=1) print $0}' $TEMPFILE> $EXPRFILE
	MODELSUFF=$OUTPUTDIR/$SPECIES
	echo "$FGCONVTOOL $EXPRFILE 1 $MODELSUFF 12"
done
