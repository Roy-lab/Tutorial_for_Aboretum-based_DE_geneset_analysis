#!/bin/sh
set -u

GK=$1		# arboretum k number
ORDER=$2	# order file location
OGID=$3		# OGID file location
TREE=$4		# tree file location
CONFIG=$5	# config file location
BEST=$6		# name of representative sample
OUTDIR=$7	# output directory name, e.g. arb_output/

P=0.8		# transitioning probability (default is 0.8)
CODE=code/Arboretum/arboretum

$CODE -s ${ORDER} -e ${OGID} -k ${GK} -t ${TREE} -c ${CONFIG} -o ${OUTDIR} -b ${BEST} \
	-r yes -m learn -i uniform -p ${P} -u true -w true

echo ""
echo " - FINISHED ARBORETUM CLUSTERING"
echo ""

