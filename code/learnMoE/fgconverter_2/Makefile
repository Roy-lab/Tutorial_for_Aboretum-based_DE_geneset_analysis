LFLAG = -lgsl -lgslcblas 
SRC = BioNetwork.C  BioNode.C  FrameWork.C  Gene.C  GeneManager.C  Graph.C  Interaction.C  InteractionManager.C  Node.C  Path.C  Protein.C  ProteinManager.C	Randomizer.C

INCLPATH = ../include
LIBPATH = ../lib
CC=g++
CFLAGS = -g

genData: $(SRC)
	$(CC) $(SRC) -I $(INCLPATH)  -L $(LIBPATH) $(LFLAG) $(CFLAGS) -o genData

