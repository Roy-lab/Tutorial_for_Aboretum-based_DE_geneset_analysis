DIR=~/data/compfuncgen/dawn_ilan_jay/
DIR=~/data/genenetweaver/multispec_learning/pg0.4_pm0.6
#for SPECIES in scer cgla scas calb klac spom
for SPECIES in Scer Sbay Klac Cgla Calb Sjap Spom
do
	#echo "./batchGenData.sh $DIR/$SPECIES/randpartitions/ $DIR/$SPECIES/randpartitions${SPECIES}"
	echo "./batchGenData.sh $DIR/$SPECIES/processed/randpartitions_n30/ $DIR/$SPECIES/processed/randpartitions_n30 ${SPECIES}"
	./batchGenData.sh $DIR/$SPECIES/processed/randpartitions_n30/ $DIR/$SPECIES/processed/randpartitions_n30 ${SPECIES}
done
