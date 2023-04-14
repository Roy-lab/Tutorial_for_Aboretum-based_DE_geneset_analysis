#!/bin/sh
set -u

#####################
## USER PARAMETERS ##
#####################
NUMSET=$1				# number of total samples
ORDER_FILE=$2			# order file for arboretum input
K=$3					# arboretum k number
SRCDIR=input_files/		# PLACE WHERE YOUR INPUT DATA VALUE FILES ARE IN
ARBSRCDIR=arb_input/	# PLACE WHERE THE ARBORETUM INPUT FILES (MAYBE YOU PREPARED) ARE IN
CODEDIR=code/init_GMM/
#####################

OUTNAME=merged_gmm_k${K}	# name of directory for init clusters

## STEP1: make merged matrix file
cat ${SRCDIR}/allgenenames.txt > merged_mat${K}.txt
for SET in  `cat $ORDER_FILE`
do
	cut -f2- ${ARBSRCDIR}/${SET}_meanval.txt > temp${K}
	paste merged_mat${K}.txt temp${K} > temp${K}.mat
	mv temp${K}.mat merged_mat${K}.txt
	rm temp${K}
done


## STEP2: prep GMM input format
COLNUM=$NUMSET
echo " - generating input data (col=${COLNUM})"
#./genData [input] [do_log_transform?] [output_prefix] [#_column] [output_format] [starting_line_index] [do_zero_mean?]
${CODEDIR}/genData merged_mat${K}.txt 1 ${OUTNAME} ${COLNUM} fgraph 0 none


## STEP3: run GMM
echo " - clustering GMM"
OUTDIR=${OUTNAME}
mkdir ${OUTDIR}
${CODEDIR}/learnMoE -m ${OUTNAME} -o ${OUTDIR} -l ${K} -e random -v 1
echo ""
echo " - ${OUTDIR} done"
rm merged_mat${K}.txt ${OUTNAME}.*


## STEP4: re-format cluster results and move GMM results to arboretum input dir
for SAMPLE in `cat $ORDER_FILE`
do
        cat ${OUTDIR}/clusterassign.txt | awk -v setname="${SAMPLE}" '{ printf setname "_" $0 "\n" }' > ${OUTDIR}/${SAMPLE}_initclust.txt
done
mv ${OUTDIR} ${ARBSRCDIR}


## STEP5: generate config file
CONFIG=${ARBSRCDIR}/config.k${K}.txt
if [ -f ${CONFIG} ]
then
	rm ${CONFIG}
fi
for SETNAME in `cat $ORDER_FILE`
do
        CLUSTER=${ARBSRCDIR}/${OUTDIR}/${SETNAME}_initclust.txt
        MATRIX=${ARBSRCDIR}/${SETNAME}_meanval.txt
        printf '%s\t%s\t%s\n' "${SETNAME}" "${CLUSTER}" "${MATRIX}" >> ${CONFIG}
done
echo " - config_file: ${CONFIG}"
echo ""
echo " - FINISHED PREP FOR ARBORETUM"
echo ""

