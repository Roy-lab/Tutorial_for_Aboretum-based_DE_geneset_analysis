#include <fstream>
#include <iostream>
#include <math.h>
#include <string.h>
#include "gsl/gsl_randist.h"
#include "Randomizer.H"

#include "Gene.H"
#include "GeneManager.H"
#include "Protein.H"
#include "ProteinManager.H"
#include "Interaction.H"
#include "InteractionManager.H"
#include "Node.H"
#include "Graph.H"
#include "BioNetwork.H"
#include "FrameWork.H"


FrameWork::FrameWork()
{
	vType=CONT;
}

FrameWork::~FrameWork()
{
}


int 
FrameWork::initDataSrc(const char* inFName)
{
	ifstream inFile(inFName);
	char buffer[1024];
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		if(strstr(buffer,"#")!=NULL)
		{
			continue;
		}
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		char keyName[256];
		char keyVal[256];
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				strcpy(keyName,tok);
			}
			else
			{
				strcpy(keyVal,tok);
			}
			tokCnt++;
			tok=strtok(NULL,"\t");
		}

		if(strcmp(keyName,"geneexpfile")==0)
		{
			geneMgr.readGeneData(keyVal);
			//geneMgr.standardizeGeneExpr();
			//geneMgr.scaleGeneExpr();
			//geneMgr.discretizeGeneExpr();
			//geneMgr.rankOrderGeneExpr();
		}
		else if(strcmp(keyName,"proteinexpfile")==0)
		{
			protMgr.setGeneMgr(&geneMgr);
			if(protMgr.readProteinData(keyVal)==-1)
			{
				return -1;
			}
			//protMgr.standardizeProteinExpr();
			//protMgr.scaleProteinExpr();
			//protMgr.discretizeProteinExpr();
			//protMgr.rankOrderProteinExpr();
		}
		else if(strcmp(keyName,"ppintfile")==0)
		{
			//P0-P1 is the same as P1-P0
			ppMgr.setDirectionality(true);
			ppMgr.setDelimiter("\t");
			ppMgr.readInteractions(keyVal);
			ppMgr.dumpInteractions("pptemp.txt");
		}
		else if(strcmp(keyName,"pdintfile")==0)
		{
			//P0-G1 is not the same as G1-P0
			pdMgr.setDirectionality(false);
			pdMgr.setDelimiter("\t");
			pdMgr.readInteractions(keyVal);
			pdMgr.dumpInteractions("pdtemp.txt");
		}
		else if(strcmp(keyName,"outputfile")==0)
		{
			strcpy(outFName,keyVal);
		}
		else if(strcmp(keyName,"pnlformatfile")==0)
		{
			strcpy(pnlGraphFName,keyVal);
		}
		else if(strcmp(keyName,"multdist")==0)
		{
			char* pos=strchr(keyVal,',');
			if(pos==NULL)
			{
				cout << "Bad format for multinomial dist" << endl;
				return -1;
			}
			*pos='\0';
			double flVal=atof(keyVal);
			mult_dist.push_back(flVal);
			flVal=atof(pos+1);
			mult_dist.push_back(flVal);
		}
		else if(strcmp(keyName,"expcntpergene")==0)
		{
			//Number of experiments per gene
			expCntPerGene=atoi(keyVal);
		}
		else if(strcmp(keyName,"genecnt")==0)
		{
			geneCnt=atoi(keyVal);
		}
		else if(strcmp(keyName,"hubgene")==0)
		{
			geneID=atoi(keyVal);
		}
		else if(strcmp(keyName,"maxpathlen")==0)
		{
			pathLen=atoi(keyVal);
		}
		else if(strcmp(keyName,"vartype")==0)
		{
			if(strcmp(keyVal,"continuous")==0)
			{
				vType=CONT;
			}
			else if(strcmp(keyVal,"discrete")==0)
			{
				vType=DISC;
			}
		}
	/*	else if(strcmp(keyName,"action")==0)
		{
			strcpy(actionStr,keyVal);	
		}*/
	}
	inFile.close();
	return 0;
}


int 
FrameWork::initDataSrc(const char* geneexpFName, int logtrans,const char* outsuffix, int aCnt,int dataStart, const char* datatransform)
{
	geneMgr.readGeneData(geneexpFName,logtrans,dataStart);
	if(strcmp(datatransform,"stnd")==0)
	{
		geneMgr.standardizeLogSpaceGeneExpr();
	}
	else if(strcmp(datatransform,"zeromean_logspace")==0)
	{
		geneMgr.zeroMeanLogSpaceGeneExpr();
	}
	else if(strcmp(datatransform,"zeromean")==0)
	{
		geneMgr.zeroMeanGeneExpr();
	}
	strcpy(outFName,outsuffix);
	sprintf(pnlGraphFName,"%s.model",outsuffix);
	expCntPerGene=aCnt;
	vType=CONT;
	return 0;
}

//Assume that in the conf file we specify a gene from where we want to build a network

//Assume that in the conf file we specify a gene from where we want to build a network
//From this gene, we find all the reachable nodes using protein-protein interactions and protein-dna interactions
int
FrameWork::generateData()
{
	BioNetwork bnw;
	bnw.setProteinManager(&protMgr);
	bnw.setGeneManager(&geneMgr);
	bnw.setPPInteractionManager(&ppMgr);
	bnw.setPDInteractionManager(&pdMgr);
	bnw.createNetwork();
	bnw.getReachability();
	bnw.showNetworkStat();
	map<int,int> measuredGeneIDs;
	map<string,int> pdIntr;
	map<string,int> ppIntr;
	Gene* gene=geneMgr.getGeneNode(geneID);
	if(bnw.getGraphInfo(gene->getName(),pathLen,measuredGeneIDs,pdIntr,ppIntr)==-1)
	{
		return -1;
	}
	//Want to generate expCnt samples each with geneCnt genes and proteins and their associated physical interactions
	// We need to obtain the experiment ids in which the genes we have selected above are deleted
/*	for(map<int,int>::iterator aIter=measuredGeneIDs.begin();aIter!=measuredGeneIDs.end();aIter++)
	{
		int gId=aIter->first;
		for(int i=0;i<expCntPerGene;i++)
		{
			int eId=(gId*expCntPerGene) + i;
			expIds.push_back(eId);
		}
			
	}*/
	//Get data for only those experiments where the hub genes is perturbed
	for(int i=0;i<expCntPerGene;i++)
	{
		int eId=(geneID*expCntPerGene) + i;
		expIds.push_back(eId);
		//expIds.push_back(i);
	}
			
	//bnw.getDiffExpGenes(expIds,geneID);		
	createGraph(measuredGeneIDs, pdIntr, ppIntr);
	phyNw.doTopologicalSort();
	vector<Node*>& topOrder=phyNw.getTopologicalSort();
	writeToFile(topOrder);
	phyNw.genPNLInputFormat(pnlGraphFName);
	return 0;
}

int
FrameWork::generateAllData()
{
	BioNetwork bnw;
	bnw.setProteinManager(&protMgr);
	bnw.setGeneManager(&geneMgr);
	bnw.setPPInteractionManager(&ppMgr);
	bnw.setPDInteractionManager(&pdMgr);
	bnw.createNetwork();
	map<int,int> measuredGeneIDs;
	map<string,int> pdIntr;
	map<string,int> ppIntr;
	//Get data for all experiments
	for(int i=0;i<expCntPerGene;i++)
	{
		expIds.push_back(i);
	}
	bnw.getAllNodeIDs(measuredGeneIDs);
	createGraph(measuredGeneIDs, pdIntr, ppIntr);
	phyNw.doTopologicalSort();
	vector<Node*>& topOrder=phyNw.getTopologicalSort();
	writeToFile(topOrder);
	writeToFile_MeanStd(topOrder);
	phyNw.genPNLInputFormat(pnlGraphFName);
	return 0;
}


int
FrameWork::generateAllData_TAB()
{
	BioNetwork bnw;
	bnw.setProteinManager(&protMgr);
	bnw.setGeneManager(&geneMgr);
	bnw.setPPInteractionManager(&ppMgr);
	bnw.setPDInteractionManager(&pdMgr);
	bnw.createNetwork();
	map<int,int> measuredGeneIDs;
	map<string,int> pdIntr;
	map<string,int> ppIntr;
	//Get data for all experiments
	for(int i=0;i<expCntPerGene;i++)
	{
		expIds.push_back(i);
	}
	bnw.getAllNodeIDs(measuredGeneIDs);
	createGraph(measuredGeneIDs, pdIntr, ppIntr);
	phyNw.doTopologicalSort();
	vector<Node*>& topOrder=phyNw.getTopologicalSort();
	writeToFile_TAB(topOrder);
	return 0;
}


/*
int
Framework::generate()
{
	if(strcmp(actionStr,"gendata")==0)
	{
		generateAllData();
	}
	else if(strcmp(actionStr,"genmodel")==0)
	{
		generateModel();
	}
	else
	{
		cout <<"Unknown actionStr " << actionStr << endl;
	}
	return 0;
}*/

//Node Ids represent protein as well as gene ids
//In this function we are going to create the graph, i.e add nodes as well as edges
//A protein node is a child of the gene node with the same id
//A protein-protein interaction node has a child for the dynamic attribute. This child
//also has protein nodes as parents
//A protein-dna interaction node has a child for the dynamic attribute. This dynamic
//node has a parent which is a protein node. This dynamic node has a gene node as child whose
//transcription level is controlled by the protein node.
int
FrameWork::createGraph(map<int,int>& nodeIds,map<string,int>& pdIntIds, map<string,int>& ppIntIds)
{
	vector<int> expLevels;
	expLevels.push_back(BioNode::LOW);
	expLevels.push_back(BioNode::MEDIUM);
	expLevels.push_back(BioNode::HIGH);

	vector<int> edgeVals;
	edgeVals.push_back(0);
	edgeVals.push_back(1);
	//For every node in the nodeId add a protein and gene node and make the protein a child of the gene node
	for(map<int,int>::iterator aIter=nodeIds.begin();aIter!=nodeIds.end();aIter++)
	{
		Gene* aGene=geneMgr.getGeneNode(aIter->first);
		const char* geneName=aGene->getName();
		phyNw.addNode(geneName,Node::GENE,aIter->first,expLevels);
		Protein* aProt=protMgr.getProteinNode(aIter->first);
		if(aProt!=NULL)
		{
			const char* proteinName=aProt->getName();
			phyNw.addNode(proteinName,Node::PROTEIN,aIter->first,expLevels);
			//Edge is added from first node to second node. This is a parent-child reln,
			//where geneName is the parent and proteinName is the child
			phyNw.addEdge(geneName,proteinName);
		}
	}
	//For every pdna interaction add a node for the static attribute, another node for dynamic 
	//attribute which is the child of the first. Then make this node and the protein node
	//as the parent of the child node
	for(map<string,int>::iterator aIter=pdIntIds.begin();aIter!=pdIntIds.end();aIter++)
	{
		char sName[256];
		sprintf(sName,"%spg.I",aIter->first.c_str());
		//Here the ids in the source data don't make sense
		phyNw.addNode(sName,Node::STAT_PG,-1,edgeVals);
		char dName[256];
		sprintf(dName,"%spg.N",aIter->first.c_str());
		phyNw.addNode(dName,Node::DYN_PG,-1,edgeVals);
		phyNw.addEdge(sName,dName);

		//Now get the protein node and the gene node and connect them
		int protId=0;
		int geneId=0;
		int srcPos=0;
		int destPos=0;
		char tempStr[256];
		const char* idPtr=aIter->first.c_str();
		while(idPtr[srcPos]!='\0')
		{
			if(idPtr[srcPos]=='-')
			{
				tempStr[destPos]='\0';
				protId=atoi(tempStr);
				destPos=0;
			}
			else
			{
				tempStr[destPos]=idPtr[srcPos];
				destPos++;
			}
			srcPos++;
		}
		tempStr[destPos]='\0';
		geneId=atoi(tempStr);

		Gene* gNode=geneMgr.getGeneNode(geneId);
		Protein* pNode=protMgr.getProteinNode(protId);
		//Make dName as a parent of gNode
		phyNw.addEdge(dName,gNode->getName());
		//Make pNode as a parent of dName 
		phyNw.addEdge(pNode->getName(),dName);
	}
	//For every protein-protein interaction do something similar
	for(map<string,int>::iterator aIter=ppIntIds.begin();aIter!=ppIntIds.end();aIter++)
	{
		char sName[256];
		//Here the ids in the source data don't make sense
		sprintf(sName,"%spp.I",aIter->first.c_str());
		phyNw.addNode(sName,Node::STAT_PP,-1,edgeVals);
		char dName[256];
		sprintf(dName,"%spp.N",aIter->first.c_str());
		phyNw.addNode(dName,Node::DYN_PP,-1,edgeVals);
		phyNw.addEdge(sName,dName);

		//Now get the protein nodes
		int protId1=0;
		int protId2=0;
		int srcPos=0;
		int destPos=0;
		char tempStr[256];
		const char* idPtr=aIter->first.c_str();
		while(idPtr[srcPos]!='\0')
		{
			if(idPtr[srcPos]=='-')
			{
				tempStr[destPos]='\0';
				protId1=atoi(tempStr);
				destPos=0;
			}
			else
			{
				tempStr[destPos]=idPtr[srcPos];
				destPos++;
			}
			srcPos++;
		}
		tempStr[destPos]='\0';
		protId2=atoi(tempStr);

		Protein* pNode1=protMgr.getProteinNode(protId1);
		Protein* pNode2=protMgr.getProteinNode(protId2);
		//Make pNode1 and pNode2 as a parents of dName
		phyNw.addEdge(pNode1->getName(),dName);
		phyNw.addEdge(pNode2->getName(),dName);
	}
	
	return 0;
}

int
FrameWork::writeToFile_PNL(vector<Node*>& topOrder)
{
	gsl_rng* g_rng=gsl_rng_alloc(gsl_rng_default);
	char exprFName[1024];
	sprintf(exprFName,"%s.data",outFName);
	ofstream oFile(exprFName);
	for(int eId=0;eId<expIds.size();eId++)
	{
		for(int i=0;i<topOrder.size();i++)
		{
			//Get the node here
			Node* anode=topOrder[i];
			//Depending upon the node type take appropriate action
			Node::NodeType nType=anode->getNodeType();
			switch(nType)
			{
				case Node::GENE:
				{
					int idInData=anode->getIdInData();
					Gene* aGene=geneMgr.getGeneNode(idInData);
					double val=aGene->getExpLevelAt(expIds[eId]);
					oFile <<" " << val;
					break;
				}
				case Node::PROTEIN:
				{
					int idInData=anode->getIdInData();
					Protein* aProtein=protMgr.getProteinNode(idInData);
					double val=aProtein->getExpLevelAt(expIds[eId]);
					oFile <<" " << val;
					break;
				}
				case Node::STAT_PG:
				case Node::STAT_PP:
				{
					//Use the multinomial distribution to set a value here
					double rNo=gsl_ran_flat(g_rng,0,1);
					if(rNo<=mult_dist[0])
					{
						oFile << " " << 1;
					}
					else
					{
						oFile << " " << 0;
					}
					break;
				}
				case Node::DYN_PP:
				case Node::DYN_PG:
				{
					//Set to 1 with 50-50 chance. It really does not matter
					double rNo=gsl_ran_flat(g_rng,0,1);
					if(rNo<0.5)
					{
						oFile <<" " << 1;
					}
					else
					{
						oFile <<" " << 1;
					}
					break;
				}
			}
		}
		oFile << endl;
	}
	oFile.close();
	return 0;
}


int
FrameWork::writeToFile(vector<Node*>& topOrder)
{
	gsl_rng* g_rng=gsl_rng_alloc(gsl_rng_default);
	char exprFName[1024];
	sprintf(exprFName,"%s.data",outFName);
	ofstream oFile(exprFName);
	for(int eId=0;eId<expIds.size();eId++)
	{
		if(eId==167)
		{
			cout <<"Stop here " << endl;
		}
		int shown=0;
		for(int i=0;i<topOrder.size();i++)
		{
			//Get the node here
			Node* anode=topOrder[i];
			//Depending upon the node type take appropriate action
			Node::NodeType nType=anode->getNodeType();
			switch(nType)
			{
				case Node::GENE:
				{
					int idInData=anode->getIdInData();
					Gene* aGene=geneMgr.getGeneNode(idInData);
				//	int dval=aGene->getDiscreteExpLevelAt(expIds[eId]);
					//double cval=aGene->getRankedExpLevelAt(expIds[eId]);
					double cval=aGene->getExpLevelAt(expIds[eId]);
					if(isnan(cval))
					{
						cout <<"Skipping " << aGene->getName() << endl;
					}
					else
					{
						if(shown>0)
						{
							oFile <<"\t";
						}
						if(vType==CONT)
						{
							oFile << i<<"=["<< cval << "]" ;
						}
						shown++;
					}
					break;
				}
				case Node::PROTEIN:
				{
					int idInData=anode->getIdInData();
					Protein* aProtein=protMgr.getProteinNode(idInData);
					//int dval=aProtein->getDiscreteExpLevelAt(expIds[eId]);
					//double cval=aProtein->getRankedExpLevelAt(expIds[eId]);
					double cval=aProtein->getExpLevelAt(expIds[eId]);
					if(i>0)
					{
						oFile <<"\t";
					}
					if(vType==CONT)
					{
						oFile << i<<"=["<< cval << "]" ;
					}
					else
					{
				//		oFile << i<<"=["<< dval << "]" ;
					}
			/*		if((aProtein->getID()==geneID) || (aProtein->getID()==87))
					{
						cout << " "<< aProtein->getName() << ": " << aProtein->getExpLevelAt(expIds[eId]);
					}*/
					break;
				}
				case Node::STAT_PG:
				case Node::STAT_PP:
				{
					//Use the multinomial distribution to set a value here
					if(i>0)
					{
						oFile <<"\t";
					}
					oFile << i<<"=[0|"<<mult_dist[1]
							 <<",1|" <<mult_dist[0]<<"]";
					break;
				}
				case Node::DYN_PP:
				case Node::DYN_PG:
				{
					//Set to 1 with 50-50 chance. It really does not matter
				/*	if(rNo<0.5)
					{
						oFile <<" " << 1;
					}
					else
					{
						oFile <<" " << 1;
					}*/
					break;
				}
			}

		}
		oFile << endl;
	}
	oFile.close();
	return 0;
}


int
FrameWork::writeToFile_TAB(vector<Node*>& topOrder)
{
	gsl_rng* g_rng=gsl_rng_alloc(gsl_rng_default);
	char exprFName[1024];
	sprintf(exprFName,"%s.geneexp",outFName);
	ofstream oFile(exprFName);
	for(int i=0;i<topOrder.size();i++)
	{
		Node* anode=topOrder[i];
		int idInData=anode->getIdInData();
		Gene* aGene=geneMgr.getGeneNode(idInData);
		oFile <<aGene->getName();
		for(int eId=0;eId<expIds.size();eId++)
		{
			double cval=aGene->getExpLevelAt(expIds[eId]);
			//oFile <<"\t" << log(cval);
			oFile <<"\t" << cval;
		}
		oFile << endl;
	}
	oFile.close();
	return 0;
}

int
FrameWork::writeToFile_MeanStd(vector<Node*>& topOrder)
{
	char exprFName[1024];
	sprintf(exprFName,"%s.meansd",outFName);
	ofstream mFile(exprFName);
	for(int i=0;i<topOrder.size();i++)
	{
		Node* anode=topOrder[i];
		int idInData=anode->getIdInData();
		Gene* aGene=geneMgr.getGeneNode(idInData);
		mFile <<aGene->getName() <<"\t" << aGene->getMean() <<"\t" << aGene->getStdev() << endl;
	}
	mFile.close();
	return 0;
}


int
main(int argc, const char** argv)
{
	if((argc!=2)&&(argc!=8))
	{
		cout << "Usage genData conffile " << endl;
		cout << "OR Usage genData expression_file logtrans[0|1] output_suffix measures_per_gene outpuformattype datastart datatransform[stnd|zermean|zeromean_logspace|none]"  << endl;
		return 0;
	}
	FrameWork fw;
	if(argc==2)
	{
		if(fw.initDataSrc(argv[1])==-1)
		{
			cout << "An error occured during initialization" << endl;
			return 0;
		}
	}
	else if(argc==8)
	{
		fw.initDataSrc(argv[1],atoi(argv[2]),argv[3],atoi(argv[4]),atoi(argv[6]),argv[7]);
	}
	if(strcmp(argv[5],"fgraph")==0)
	{
		fw.generateAllData();
	}
	else if(strcmp(argv[5],"tab")==0)
	{
		fw.generateAllData_TAB();
	}
	else
	{
		cout <<"Unable to understand format" << argv[5] << endl;
	}
	return 0;
}
