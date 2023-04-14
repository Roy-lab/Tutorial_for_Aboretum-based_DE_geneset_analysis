INDIR=../../../data/genenetweaver/p0.4/
#for N in 100 200 300 400 500
for N in 1000
do
for sparsity in 0.1 0.15
do
SUBDIR=$INDIR/net${N}/sparsity_p${sparsity}/processed
EXPCNT=`head -n1 $SUBDIR/net${N}_geneexp.txt | awk '{print NF}'`
EXPCNT=`expr $EXPCNT - 1`
echo "./genData $SUBDIR/net${N}_geneexp.txt 1 $SUBDIR/net${N} $EXPCNT fgraph 0 none"
./genData $SUBDIR/net${N}_geneexp.txt 1 $SUBDIR/net${N} $EXPCNT fgraph 0 none
done
done
