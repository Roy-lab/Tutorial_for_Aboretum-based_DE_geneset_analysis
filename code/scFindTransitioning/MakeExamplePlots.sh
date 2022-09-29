export HMAWK=/mnt/ws/sysbio/roygroup/shared/scripts/figscripts/Heatmap.awk

cat example_out/ordered_clusterset_means.txt | ${HMAWK} -vStrokeC="-" -vStrokeSC="black" -vL="ClusterAssignment Activity" -vD="space" -vLegendNum="6" -vC="1:(0,205,0);3:(255,255,0);5:(205,0,0) -2.0:(0,255,0);0:(0,0,0);2.0:(255,0,0)" -vFontSize=8 > example_plots/ordered_clusterset_means.svg

cat example_out/clusterset100.txt | ${HMAWK} -vStrokeC="-" -vStrokeSC="black" -vL="ClusterAssignment Activity" -vD="space" -vLegendNum="6" -vC="1:(0,205,0);3:(255,255,0);5:(205,0,0) -2.0:(0,255,0);0:(0,0,0);2.0:(255,0,0)" -vFontSize=8 > example_plots/clusterset100.svg

cat example_out/all_assign.txt | ${HMAWK} -vStrokeC="-" -vStrokeSC="black" -vL="ClusterAssignment Activity" -vD="space" -vLegendNum="6" -vC="1:(0,205,0);3:(255,255,0);5:(205,0,0) -2.0:(0,255,0);0:(0,0,0);2.0:(255,0,0)" -vFontSize=8 > example_plots/all_assign.svg
