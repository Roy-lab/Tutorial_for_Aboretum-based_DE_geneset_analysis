#ifndef _FRAMEWORK
#define _FRAMEWORK

#include <vector>

using namespace std;

#define CONT 1
#define DISC 2

class InteractionManager;
class ProteinManager;
class GeneManager;
class Graph;

class FrameWork
{
	public:
		FrameWork();
		~FrameWork();
		int initDataSrc(const char*);
		int initDataSrc(const char* geneexpFName,int logtrans, const char* outsuffix, int aCnt,int,const char* transform);
		int generateData();
		int generateAllData();
		int generateAllData_TAB();
	private:
		int createGraph(map<int,int>&, map<string,int>&,map<string,int>&);
		int writeToFile_PNL(vector<Node*>&);
		int writeToFile(vector<Node*>&);
		int writeToFile_TAB(vector<Node*>& topOrder);
		int writeToFile_MeanStd(vector<Node*>& topOrder);
		char outFName[256];
		char pnlGraphFName[256];
		//This is for interactions that are present
		//The less uniform this distribution, the closer we are to what we know
		vector<double> mult_dist;
		InteractionManager pdMgr;
		InteractionManager ppMgr;
		ProteinManager protMgr;
		GeneManager geneMgr;
		int expCntPerGene;
		vector<int> expIds;
		//For now number of genes also specifies the number of proteins
		int geneCnt;
		
		//The ID of the gene from which we will trace paths
		int geneID;
		//The max path length we will use to find reachable nodes
		int pathLen;
		Graph phyNw;
		int vType;
};
#endif
