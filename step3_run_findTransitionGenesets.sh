#!/bin/sh
set -u

#####################
## USER PARAMETERS ##
#####################
ARBOUT=$1	# arboretum output dir
ORDER=$2	# arboretum input order file location
OGID=$3		# arboretum input OGID file location
BEST=$4		# name of representative sample
OUTDIR=$5	# output dir name
SCRIPT_DIR=script/
CODE_DIR=code/
#####################

## STEP1: reorder arboretum cluster ID / generate input format
${SCRIPT_DIR}/reorder_reformat_arboretum_result.pl ${ARBOUT}/


## STEP2: run find_transition_genesets
${CODE_DIR}/scFindTransitioning/findTransitionGenesets_miss $ARBOUT $ORDER $OGID $BEST 0.05 $OUTDIR 5 $ORDER 0

echo ""
echo " - FINISHED FIND TRANSIONING GENESET CLUSTERING"
echo ""

