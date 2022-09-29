set -u

INPUTFILE=$1	# input filename
GK=$2		# k number of arboretum
MAX_EXPR=$3	# max value for expression
SCRIPT_DIR=script/

INPUTFN=$INPUTFILE
OUTFN="${INPUTFN%.*}"

MID_GK=$((($GK + 1) / 2))
MID_EXPR=$(($MAX_EXPR / 2))
cat ${INPUTFN} | ${SCRIPT_DIR}/Heatmap.awk -vStrokeC="-" -vStrokeSC="black" -vL="ClusterAssignment Expression Missing" -vD="space" -vLegendNum="${GK}" -vC="1:(0,205,0);${MID_GK}:(255,255,0);${GK}:(205,0,0) 0:(0,255,0);${MID_EXPR}:(0,0,0);${MAX_EXPR}:(255,0,0) 0:(105,105,105)" -vFontSize=8 > ${OUTFN}.svg


