if [ $# -ne 2 ]
then
	echo "Usage: batchGen inputdir suffix"
	exit
fi

INPUTDIR=$1
SUFFIX=$2
OLDSETTINGS=dgen_dream.conf
NEWSETTINGS=dgen_dream_new.conf
for DATASETID in 1 2 3 4 5
do
	EXPRFILE=$INPUTDIR/net${DATASETID}/$SUFFIX.txt
	MODELFILE=$INPUTDIR/net${DATASETID}/$SUFFIX.model
	DATAFILE=$INPUTDIR/net${DATASETID}/$SUFFIX.data
	export EXPRFILE
	export MODELFILE
	export DATAFILE
	awk '{ if($1~/geneexpfile/) printf("geneexpfile\t%s\n",ENVIRON["EXPRFILE"])
		else if($1~/outputfile/) printf("outputfile\t%s\n",ENVIRON["DATAFILE"])
		else if($1~/pnlformatfile/) printf("pnlformatfile\t%s\n",ENVIRON["MODELFILE"])
		else printf("%s\n",$0)}' $OLDSETTINGS > $NEWSETTINGS
	cat $NEWSETTINGS
	echo "./genData $NEWSETTINGS"
	./genData $NEWSETTINGS
done
