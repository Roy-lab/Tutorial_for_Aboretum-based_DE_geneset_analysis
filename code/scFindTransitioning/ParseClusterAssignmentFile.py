
import sys
import re
import math
if(len(sys.argv)<5):
    print 'usage ParseClusterAssignmentFile.py [assignment_data] [map_list] [species_list] [outputfile]'
    sys.exit(1)

print 'Reading in cluster assignemnt file\n'

CADATA=open(sys.argv[1]).readlines()
CADATA_D = {}

for l in range(0,len(CADATA)):
        line=CADATA[l].strip().split('\t')
	gene=line[0]
	assignments=line[1:]
        CADATA_D[gene]=assignments

print len(CADATA_D)
	
print 'Reading in assignemnt map\n'

MAP=open(sys.argv[2]).readlines()
MAP_L=[]
MAP_D={}
for l in range(0,len(MAP)):
        mappedid=int(MAP[l])
	MAP_L.append(mappedid-1)
	MAP_D[mappedid-1]=l
print MAP_L

SPC=open(sys.argv[3]).readlines()
SPC_L=[]
for l in range(0,len(SPC)):
        line=SPC[l].strip().split('\t')
        species=line[0]
        SPC_L.append(species)
print SPC_L

oFile=open(sys.argv[4],'w')
oFile.write('Loci')
for species in SPC_L:
	oFile.write('\t'+species)
oFile.write('\n')
for gene in CADATA_D:
	if "Dumm" in gene:
		continue
	oFile.write(gene)
	for l in range(0,len(CADATA_D[gene])):	
		ca=int(CADATA_D[gene][l])
		if ca >= 0:
			ca_o=int(MAP_D[ca])
			oFile.write('\t'+str(ca_o))
		if ca < 0:
			oFile.write('\t'+str(ca))
	oFile.write('\n')
