#!/usr/bin/env bash
set -u

ORDER=$1
GENELIST=$2

join_arr() {
	local IFS="$1"
	shift
	echo "$*"
}

i=0
for GENE in `cat ${GENELIST}`
do
	((i=i+1))
	printf '%s\t' "OG${i}_0"
	array=()
	for SETNAME in `cat ${ORDER}`
	do
		array+=("${SETNAME}_${GENE}")
	done
	join_arr , "${array[@]}"
done


