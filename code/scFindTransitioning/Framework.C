#include <iostream>
#include <fstream>
#include <vector>
#include <string.h>
#include <stdlib.h>
#include "HierarchicalClusterNode.H"
#include "HierarchicalCluster.H"
#include "OptimalLeafOrder.H"
#include "GeneExpManager.H"

#include "GeneMap.H"
#include "MappedOrthogroup.H"
#include "MappedOrthogroupReader.H"

#include "Framework.H"
using namespace std;
Framework::Framework()
{
}

Framework::~Framework()
{
}

/*This function will read in the name of the clusterassignment file and also the mark values that we will need. We need to store the mark values column wise as per
 * the columns of the data matrix */
int 
Framework::readDataMatrix(const char* aDirName, int maxmissing)	//JSmod
{
	char aFName[1024];
	sprintf(aFName,"%s/allcelltypes_clusterassign_brk.txt",aDirName);
	ifstream inFile(aFName);
	char* buffer=NULL;
	int bufflen=0;
	string strbuff;
	int gid=1;
	int totalLoci=0;
	while(inFile.good())
	{
		getline(inFile,strbuff);
		if(strbuff.length()<=0)
		{
			continue;
		}
		if(bufflen<=strbuff.length())
		{
			if(buffer!=NULL)
			{
				delete[] buffer;
			}
			bufflen=strbuff.length()+1;
			buffer=new char[bufflen];
		}
		strcpy(buffer,strbuff.c_str());
		if(strstr(buffer,"Loci")!=NULL)
		{
			readColumns(buffer);
			continue;
		}
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;	
		string geneName;
		vector <double> valueSet;
		map <int,int> frequency;
		int misses=0;
		int total=0;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				geneName.append(tok);
			}
			else
			{
				double v=atof(tok);
				if(v<0)
				{
					misses+=1;
				}
				total+=1;
				/*if(v>0)
				{
					node->attrib[tokCnt-1]=1;
				}
				else
				{
					node->attrib[tokCnt-1]=0;
				}
				node->attrib[tokCnt-1]=v+1;*/
				valueSet.push_back(v+1);
				int intVal=(int) (v+1);	
				if(frequency.find(intVal)==frequency.end())
				{
					frequency[intVal]=1;
				}
				else
				{
					frequency[intVal]=frequency[intVal]+1;
				}
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		totalLoci++;
		if(frequency.size()==1)
		{
			frequency.clear();
			valueSet.clear();
			continue;
		}
		if(misses>maxmissing)	//JSmod
		{
			continue;
		}	
		HierarchicalClusterNode* node=new HierarchicalClusterNode;
		for(int i=0;i<valueSet.size();i++)
		{
			node->expr.push_back(valueSet[i]);
			//node->attrib[i]=valueSet[i];
		}
		node->nodeName.append(geneName.c_str());
		nodeSet[geneName]=node;
		backup[geneName]=node;
		nameIDMap[geneName]=gid;
		node->size=1;
		gid++;
		frequency.clear();
		valueSet.clear();
	}
	cout <<"Found " << nodeSet.size() << " loci that change cluster assignment of total " << totalLoci << endl;
	if(attribNameIDMap.size()==0)
	{
		cout << "Did not find any cell types!" <<endl;
		exit(0);
	}
	inFile.close();
	//Now read the mark profiles. This should be read in the same order as the order of attribIDNameMap
	for(map<int,string>::iterator cIter=attribIDNameMap.begin();cIter!=attribIDNameMap.end();cIter++)
	{
		if(cIter->second.find("Anc")==0)
		{
			continue;
		}
		char geneFName[1024];
		sprintf(geneFName,"%s/%s_exprtab.txt",aDirName,cIter->second.c_str());
		GeneExpManager* geExp=new GeneExpManager;
		geExp->readExpression_Withheader(geneFName);
		markProfileSet[cIter->second]=geExp;
		//cout << cIter->second << endl;
	}
	cout <<"Read " << markProfileSet.size() << " different datasets " << endl;
	createMarkHeaders();
	return 0;
}

int
Framework::readOGIDs(const char* celltypeOrder,const char* aFName)
{
	//mor.setSpeciesMapping(attribIDNameMap);
	mor.readSpeciesMapping(celltypeOrder);
	mor.readFile(aFName);
	return 0;
}

int
Framework::setSrcCellType(const char* src)
{
	srcCelltype.append(src);
	return 0;
}

int 
Framework::readSpeciesOrder(const char* sFile )
{
	cout << "Reading SpeciesOrder from "<< sFile << endl;
        string buffstr;
        ifstream inFile(sFile);
        char* buffer = new char[1024];
        int bufflen=1024;
        while(inFile.good()){
                inFile.getline(buffer,1023);
                if(strlen(buffer)<=0)
                {
                        continue;
                }
                if(bufflen<=buffstr.length())
                {
                        delete[] buffer;
                        bufflen=buffstr.length()+1;
                        buffer=new char[bufflen];
                }
                int tokCnt=0;
                char* tok=strtok(buffer,"\t");
                //define variables to parse out
                string name;
                //parse string into the needed variables.
                while(tok!=NULL)
                {
                        if(tokCnt==0)
                        {
                                name.append(tok);
                        }
                        tok=strtok(NULL,"\t");
                        tokCnt++;
                }
                speciesOrder.push_back(name);
        }
        inFile.close();
}

//This function will generate one file which is the assignment of genes in all cell types and this will include everything. 
//The second set of files will be the clustersets that will include all the individual transitions. I think we should also generate the cluster mark profiles for the clusterset.
//It's better to output in heatmap.awk compatible format right here.
int 
Framework::generateTransitioningGeneSets(double threshold,const char* outdir,int mingenesetSize)
{
	cout << "Entering Framework::generateTransitioningGeneSets" << endl;
	map<int,map<string,int>*> modules;
	map<int,HierarchicalClusterNode*> attribs;
	cluster.setDistanceType(HierarchicalCluster::CITYBLOCK);
	cluster.cluster(modules,nodeSet,threshold,attribs);
	double** dist=cluster.getDist();
	char aFName[1024];
	sprintf(aFName,"%s/all_assign.txt",outdir);
	ofstream oFile(aFName);

	// DC add: create the geneset file
	char gFName[1024];
	sprintf(gFName, "%s/all_genesets.txt",outdir);
	ofstream gFile(gFName);
	//SK:  add file that is purely the cluster assignment matrix
	char bFName[1024];
        sprintf(bFName,"%s/all_genes_clusterassignment_matrix.txt",outdir);
        ofstream bFile(bFName);
	bFile << "Gene";
	for(int i=0;i<attribIDNameMap.size();i++)
	{
		bFile << "\t" << attribIDNameMap[i];
	}
	bFile << endl;
	int c=0;
	int atcnt=0;
	olo.setDist(dist);
	for(map<int,HierarchicalClusterNode*>::iterator aIter=attribs.begin();aIter!=attribs.end();aIter++)
	{
		//SR added this check of 100 elements in the set.
		map<string,int>* members=modules[c];
		vector <string>* ordering =new vector<string>;
		if(members->size()>100)
		{
			for(map<string,int>::iterator mIter=members->begin();mIter!=members->end();mIter++)
			{
				ordering->push_back(mIter->first);
			}
		}
		else
		{
			olo.setHierarchicalClusterNode(aIter->second);
			olo.reorder(*ordering);
		}
		for(int i=0;i<ordering->size();i++)
		{
			HierarchicalClusterNode* hc=backup[(*ordering)[i]];
			//for(map<int,double>::iterator aIter=hc->attrib.begin();aIter!=hc->attrib.end();aIter++)
			for(int j=0;j<hc->expr.size();j++)
			{
				string& colName=attribIDNameMap[j];
				if(hc->expr[j]>=0)oFile<<(*ordering)[i]<< "||6\t" <<colName <<"||8\t" << hc->expr[j]<<"|1|"<< hc->expr[j]<< endl;
				else oFile<<(*ordering)[i]<< "||6\t" <<colName <<"||8\t" << hc->expr[j]<<"|3|"<< hc->expr[j]<< endl;
			}
			//Now put a vertical spacer
			//we need this only once
			if(c==0)
			{
				oFile <<"|- MyVspacer1|Spacer"<< endl;
			}
			//Not print out out the mark profiles too
			//for(map<int,double>::iterator aIter=hc->attrib.begin();aIter!=hc->attrib.end();aIter++)
			////for(int j=0;j<speciesOrder.size();j++)
			//SR turning this around from printing in speciesOrder to the attribute order
			for(int j=0;j<hc->expr.size();j++)
			{
				//string& colName=speciesOrder[j];//attribIDNameMap[j];
				string& colName=attribIDNameMap[j];
				//cout << colName << endl;
				GeneExpManager* geMgr=markProfileSet[colName];
				vector<string>& markNames=geMgr->getColNames();
				//here we have to be careful because we don't know what the locus name in the other cell type is
				if(colName.find("Anc")!=std::string::npos)
                                {
                                        continue;
                                }
				const char* lociCelltype=NULL;
				if(colName!=srcCelltype)
				{
					lociCelltype=getNameinCelltype((*ordering)[i].c_str(),colName.c_str());
				}
				else
				{
					lociCelltype=(*ordering)[i].c_str();
				}
				if(lociCelltype==NULL || hc->expr[i]<0)
				{
					for(int m=0;m<markNames.size();m++)
                                        {
                                                oFile<< (*ordering)[i]<<"||6\t"<< colName << "_" << markNames[m] << "||6\t" << -100 <<"|3" << endl;
                                        }
					oFile <<"|- cellmarkerVspacer"<<(j+1)<<"|Spacer"<< endl;
					continue;
				}
				vector<double>* expvals=geMgr->getExp(lociCelltype);
				if(expvals!=NULL)
				{
					for(int m=0;m<markNames.size();m++)
					{
						oFile<< (*ordering)[i]<<"||6\t"<< colName << "_" << markNames[m] << "||6\t" << (*expvals)[m]<<"|2" << endl;
					}
				}
				else
				{
					for(int m=0;m<markNames.size();m++)
                                        {
						oFile<< (*ordering)[i]<<"||6\t"<< colName << "_" << markNames[m] << "||6\t" << -100 <<"|3" << endl;
					}
				}
				if(c==0 && i==0 && (j<hc->expr.size()-1))
				{
					oFile <<"|- cellmarkerVspacer"<<(j+1)<<"|Spacer"<< endl;
				}
			}
			atcnt=hc->expr.size();	
		}
		//Now put a horizontal spacer		
		oFile <<"|Spacer||"<<c << " |-"<< endl;
		if((*ordering).size()>=mingenesetSize)
		{
			// DC
			// Write out the gene set for this cluster
			// ClusterID	gene#gene#...gene
			gFile << "Cluster" << c << "\t" << (*ordering)[0];
			for(int i=1;i<ordering->size();i++)
			{
				gFile << "#" << (*ordering)[i];
			}
			gFile << endl;
			// end DC update
			clusterset[c]=ordering;
			char cFName[1024];
			sprintf(cFName,"%s/clusterset%d.txt",outdir,c);
                        ofstream coFile(cFName);
			//Show header
			for(int i=0;i<ordering->size();i++)
			{
				HierarchicalClusterNode* hc=backup[(*ordering)[i]];	
				//for(map<int,double>::iterator aIter=hc->attrib.begin();aIter!=hc->attrib.end();aIter++)
				bFile << (*ordering)[i];
				for(int j=0;j<hc->expr.size();j++)
				{
					bFile << "\t" << hc->expr[j];
					if(hc->expr[j]>0)coFile<<(*ordering)[i] <<"||6\t"<< attribIDNameMap[j] << "||8\t" << hc->expr[j]<<"|1|"<< hc->expr[j]<<endl;
					else coFile<<(*ordering)[i] <<"||6\t"<< attribIDNameMap[j] << "||8\t" << hc->expr[j]<<"|3|"<< hc->expr[j]<<endl;
				}
				bFile << endl;
				coFile <<"|- MyVspacer1|Spacer"<< endl;
				//Not print out out the mark profiles too
				//for(int j=0;j<speciesOrder.size();j++)
				//SR changing this to print the values in the order of cluster assginment.
				for(int j=0;j<hc->expr.size();j++)
				{
					//string& colName=speciesOrder[j]; //attribIDNameMap[j];
					string& colName=attribIDNameMap[j];
					GeneExpManager* geMgr=markProfileSet[colName];
					vector<string>& markNames=geMgr->getColNames();
					//vector<double>* expvals=geMgr->getExp((*ordering)[i]);
					if(colName.find("Anc")!=std::string::npos)
                                	{
                                        	continue;
                                	}
					const char* lociCelltype=NULL;
                                	if(colName!=srcCelltype)
                                	{
                                        	lociCelltype=getNameinCelltype((*ordering)[i].c_str(),colName.c_str());
                                	}
                                	else
                                	{
                                        	lociCelltype=(*ordering)[i].c_str();
                                	}
					//const char* lociCelltype=getNameinCelltype((*ordering)[i].c_str(),colName.c_str());
					if(lociCelltype==NULL)
					{
						for(int m=0;m<markNames.size();m++)
						{
							coFile<< (*ordering)[i]<<"||6\t"<< colName << "_" << markNames[m] << "||6\t" << -100 <<"|3" << endl;
							//cout << "line " << 203 << " " << (*ordering)[i]<<"||6\t"<< colName << "_" << markNames[m] << "||6\t" << -100 <<"|3" << endl;
						}
						coFile <<"|- cellmarkerVspacer"<<(j+1)<<"|Spacer"<< endl;
						continue;
					}
                                	vector<double>* expvals=geMgr->getExp(lociCelltype);
					if(expvals!=NULL)
					{
						for(int m=0;m<markNames.size();m++)
						{
							coFile<< (*ordering)[i]<<"||6\t"<< colName << "_" << markNames[m] << "||6\t" << (*expvals)[m]<<"|2" << endl;
							//cout << "line " << 203 << " " << (*ordering)[i]<<"||6\t"<< colName << "_" << markNames[m] << "||6\t" << (*expvals)[m]<<"|2" << endl;
						}
					}
					else
					{
						for(int m=0;m<markNames.size();m++)
						{
							coFile<< (*ordering)[i]<<"||6\t"<< colName << "_" << markNames[m] << "||6\t" << -100 <<"|3" << endl;
							//#cout << "line " << 203 << " " << (*ordering)[i]<<"||6\t"<< colName << "_" << markNames[m] << "||6\t" << -100 <<"|3" << endl;
						}
					}
					coFile <<"|- cellmarkerVspacer"<<(j+1)<<"|Spacer"<< endl;
				}
			}
			coFile.close();
			bFile << "DUMMY" << c;
                	for(int x=0;x<attribIDNameMap.size();x++)
                	{
                        	bFile << "\t-100";
                	}
                	bFile << endl;
		}
		c++;
	}
	oFile.close();
	bFile.close();
	gFile.close(); // Close geneset file
	cout << "Exiting Framework::generateTransitioningGeneSets" << endl;
	return 0;
}

int
Framework::generateOrderedClusterMeans(const char* outDirName)
{
	cout << "Entering Framework::generateOrderedClusterMeans" << endl;
	map<string,vector<double>*> meanC_Set;
	map<string,vector<int>*> meanD_Set;
	map<string,HierarchicalClusterNode*> meanNodeSet;
	int gid=0;
	//We need to compute the means, and reorder them for the clusters that are of reasonable size
	cout << "Size of Clusterset is " << clusterset.size() << endl;;
	bool first=true;
	for(map<int,vector<string>*>::iterator aIter=clusterset.begin();aIter!=clusterset.end();aIter++)
	{
		vector <string>* clusterMembers=aIter->second;
		vector <double>* meanC=getMeanContinuous(clusterMembers,first);
		vector <int>*  meanD=getMeanDiscrete(clusterMembers);
		HierarchicalClusterNode* node = new HierarchicalClusterNode;
		for(int i=0;i<meanC->size();i++)
		{
			node->expr.push_back((*meanC)[i]);
		}
		char cname[1024];
		sprintf(cname,"clusterset%d",aIter->first);
		node->nodeName.append(cname);
		/*cout << node->nodeName << endl;
		for(int i=0;i<meanC->size();i++)
		{
			cout << (*meanC)[i] << endl ;
		}
		for(int i=0;i<meanD->size();i++)
                {
                        cout << (*meanD)[i] << endl ;
                }*/
		meanNodeSet[cname]=node;
		meanC_Set[cname]=meanC;
		meanD_Set[cname]=meanD;
		nameIDMap[cname]=gid;
		node->size=1;
		gid++;
		first=false;
	}
	cout << "Obtained means" << endl;
	map<int,map<string,int>*> modules;
	HierarchicalCluster clusterMeans;
	clusterMeans.setDistanceType(HierarchicalCluster::PEARSON);
	clusterMeans.cluster(modules,meanNodeSet,1);
	olo.setHierarchicalClusterNode(clusterMeans.getRoot());
	vector <string> ordering;
	olo.reorder(ordering);
	cout << ordering.size() << endl;
	char aFName[1024];
	sprintf(aFName,"%s/ordered_clusterset_means.txt",outDirName);
	ofstream oFile(aFName);
	for(int i=0;i<ordering.size();i++)
	{
		
		cout << ordering[i] << endl;
		HierarchicalClusterNode* hc=meanNodeSet[ordering[i]];
		vector<int>* meanD=meanD_Set[ordering[i]];
		for(int j=0;j<meanD->size();j++)
		{
			if((*meanD)[j]>0)
			{
				//cout << (*ordering)[i] << "||8\t" << attribIDNameMap[j] <<"||8\t" << (*meanD)[j]<<"|1|"<<(*meanD)[j]<< endl;
				oFile << ordering[i] << "||8\t" << attribIDNameMap[j] <<"||8\t" << (*meanD)[j]<<"|1|"<<(*meanD)[j]<< endl;
			}
			else
			{
				//cout << (*ordering)[i] << "||8\t" << attribIDNameMap[j] <<"||8\t" << (*meanD)[j]<<"|3|"<<(*meanD)[j]<< endl;
				oFile << ordering[i] << "||8\t" << attribIDNameMap[j] <<"||8\t" << (*meanD)[j]<<"|3|"<<(*meanD)[j]<< endl;
			}
		}
		for(map<int,string>::iterator aIter=attribIDNameMap.begin();aIter!=attribIDNameMap.end();aIter++)
		{
			//Get the cellName 
			string& cellName=(string&)aIter->second; 
			if(cellName.find("Anc")!=std::string::npos)
                	{
                        	continue;
                	}
			if(i==0 && aIter->first<attribIDNameMap.size())
                        {       
                                //cout << "|- cellmarkerVspacer" << (aIter->first+1) << "|Spacer" << endl;      
                                oFile << "|- cellmarkerVspacer" << (aIter->first+1) << "|Spacer" << endl;
                        }
			//Get the marks
			GeneExpManager* geMgr=markProfileSet[cellName];
			vector<string>& colNames=geMgr->getColNames();
			char cellMark[1024];
			for(int j=0;j<colNames.size();j++)
			{
				sprintf(cellMark,"%s_%s",cellName.c_str(),colNames[j].c_str());
				string attribName(cellMark);
				int markID=attribNameIDMap_Profile[attribName];
				if(hc->expr[markID]!=-100) oFile << ordering[i] << "||8\t" << meanLabels[markID] <<"||8\t" << hc->expr[markID] << "|2" << endl;
				else oFile << ordering[i] << "||8\t" << meanLabels[markID] << "||8\t" << hc->expr[markID] << "|3" << endl;
				//if(hc->expr[markID]!=-100) cout << (*ordering)[i] << "||8\t" << attribName <<"||8\t" << hc->expr[markID] << "|2" << endl;
                                //else cout << (*ordering)[i] << "||8\t" << attribName << "||8\t" << hc->expr[markID] << "|3" << endl;
			}
		}
	}
	oFile.close();
	char bFName[1024];
        sprintf(bFName,"%s/ordered_mean_clusterassign_matrix.txt",outDirName);
        ofstream mFile(bFName);
	mFile << "Species";
	//prep just one initial example:
	cout << ordering[0] << endl;
        HierarchicalClusterNode* hc=meanNodeSet[ordering[0]];
        vector<int>* meanD=meanD_Set[ordering[0]];
	for(int j=0;j<meanD->size();j++)
	{
		mFile << "\t" << attribIDNameMap[j];
	}
	mFile << endl;
	for(int i=0;i<ordering.size();i++)
        {

                cout << ordering[i] << endl;
                HierarchicalClusterNode* hc=meanNodeSet[ordering[i]];
                vector<int>* meanD=meanD_Set[ordering[i]];
		mFile << ordering[i];
                for(int j=0;j<meanD->size();j++)
                {
                	mFile << "\t" << (*meanD)[j];
                }
		mFile << endl;
	}
	cout << "Exiting Framework::generateOrderedClusterMeans" << endl;
	return 0;
}

//This will create a pseudo header using the ordering of the cell types in the cmint assignment set
int
Framework::createMarkHeaders()
{
	int markID=0;
	for(map<int,string>::iterator aIter=attribIDNameMap.begin();aIter!=attribIDNameMap.end();aIter++)
	{
		if(aIter->second.find("Anc")==0)
		{
			continue;
		}
		//for(int x=0;x<speciesOrder.size();x++)
		//{		
		string& cellName=(string&)aIter->second; 
		//Get the marks
		GeneExpManager* geMgr=markProfileSet[cellName];
		vector<string>& colNames=geMgr->getColNames();
		char cellMark[1024];		
		for(int i=0;i<colNames.size();i++)
		{
			sprintf(cellMark,"%s_%s",cellName.c_str(),colNames[i].c_str());
			string attribName(cellMark);
			//cout << markID << "\t" << attribName;
			attribIDNameMap_Profile[markID]=attribName;
			attribNameIDMap_Profile[attribName]=markID;
			markID++;	
		}
	}
	return 0;
}

vector<double>*
Framework::getMeanContinuous(vector<string>* geneSet,bool fillLabels)
{
	vector<double>* theMean=new vector<double>;	
	//get the cell types in the order of the cellNameIDMap
	//for(map<int,string>::iterator aIter=attribIDNameMap.begin();aIter!=attribIDNameMap.end();aIter++)
	int count=0;
	if(fillLabels && meanLabels.size()>0)
	{
		meanLabels.clear();
	}	
	for(int x=0;x<speciesOrder.size();x++)
	{
		string& colName=speciesOrder[x]; //aIter->second;
		if(colName.find("Anc")!=std::string::npos)
                {
			continue;
                }
		GeneExpManager* geMgr=markProfileSet[colName];
		vector<string>& colNames=geMgr->getColNames();
		for(int j=0;j<colNames.size();j++)
		{
			if(fillLabels)
			{
				string label;
				label.append(speciesOrder[x]);
				label.append("_");
				label.append(colNames[j]);
				meanLabels.push_back(label);
				count++;
			}	
			double s=0;
			//SK: initialize for count of genes with data.
			double c=0;
			for(int i=0;i<geneSet->size();i++)
			{
				//SK: update for cases of missing genes
				const char* lociCelltype=NULL;
				if(colName!=srcCelltype)
				{
					//cout << "called for " << colName << endl;
					lociCelltype=getNameinCelltype((*geneSet)[i].c_str(),colName.c_str());
				}
				else
				{
					lociCelltype=(*geneSet)[i].c_str();
                                }
				if(lociCelltype==NULL)
				{
					continue;
				}
				vector<double>* vals=geMgr->getExp(lociCelltype);
				if(vals!=NULL)
				{
					s = s+(*vals)[j];	
					c += 1;
				}	
			}
			//SK: count average with respect to the number of genes with data and not the number of genes perse.
			if(c>=1)
			{
				s=s/c; //geneSet->size();
				theMean->push_back(s);
			}
			else{
				theMean->push_back(-100);
			}
		}
	}
	return theMean;
}

vector<int>*
Framework::getMeanDiscrete(vector<string>* geneSet)
{
	vector<int>* theMean=new vector<int>;
	for(map<int,string>::iterator aIter=attribIDNameMap.begin();aIter!=attribIDNameMap.end();aIter++)
	{
		map<int,int> frequency;
		for(int i=0;i<geneSet->size();i++)
		{
			string& gname=(*geneSet)[i];
			HierarchicalClusterNode* hc=backup[gname];
			int val=(int)hc->expr[aIter->first];
			if(frequency.find(val)==frequency.end())
			{
				frequency[val]=1;
			}
			else
			{
				frequency[val]=frequency[val]+1;
			}
		}
		//Now get the majority. In case of a tie we will use the max of the assignments (hopefully there won't be ties!)
		int maxCnt=-1;
		int maxVal=-1;
		for(map<int,int>::iterator cIter=frequency.begin();cIter!=frequency.end();cIter++)
		{
			if(cIter->second>maxCnt)
			{
				maxCnt=cIter->second;
				maxVal=cIter->first;
			}
		}
		//Now check for ties
		for(map<int,int>::iterator cIter=frequency.begin();cIter!=frequency.end();cIter++)
		{
			if(cIter->second==maxCnt && maxVal!=cIter->first)
			{
				//cout <<"Oops found a tie of " << maxVal<< " with " << cIter->first << endl;
				if(cIter->first>maxVal)
				{
					maxVal=cIter->first;
				}
			}
		}
		frequency.clear();
		theMean->push_back(maxVal);
	}
	return theMean;
}



int
Framework::readColumns(char* buffer)
{
	int tokCnt=0;
	char* tok=strtok(buffer,"\t");
	while(tok!=NULL)
	{
		if(tokCnt>0)
		{
			string colName(tok);
			attribNameIDMap[colName]=tokCnt-1;
			attribIDNameMap[tokCnt-1]=colName;
		}
		tok=strtok(NULL,"\t");
		tokCnt++;
	}
	return 0;
}

const char*
Framework::getNameinCelltype(const char* locus, const char* celltype)
{
	//cout << "Touch Framework::getNameinCelltype " << "\t" << locus << "\t" << celltype << endl;
	const char* name=NULL;
	//try first to get orthogroup based on source species gene name 
	MappedOrthogroup* mgrp;
	string in;
	in.append(locus);
	mgrp=mor.getMappedOrthogroup(locus,srcCelltype.c_str());
	if(mgrp!=NULL)
	{
		map<int,map<string,string>*> geneSets=mgrp->getGeneSets();
		for(map<int,map<string,string>*>::iterator lIter=geneSets.begin();lIter!=geneSets.end();lIter++)
		{
			map<string,string>* levelSet=lIter->second;
			//adding second condition. 
			if(levelSet->find(celltype)!=levelSet->end() && levelSet->find(srcCelltype)!=levelSet->end() && levelSet->find(srcCelltype)->second==in)
			{
				name=levelSet->find(celltype)->second.c_str();
				break;		
			}
		}
	}
	//a check for the chiclids case, not generally useful for all species. 
	//if(mgrp==NULL && locus[0]=='a' && locus[1]=='b')
	//{
	//	cout << "Asbu gene not matched to an OGID\n";
	//}
	//ogid group was not found with a gene name so check if it's a OG id string. 
	if(mgrp==NULL && locus[0]=='O' && locus[1]=='G' && isdigit(locus[2])) 
	{  
		char* pos=strchr((char*)locus,'*');
        	int ogid=-100;
		int dup=-100;
        	if(pos!=NULL)
        	{      
               		*pos='\0';
        	}
		pos=strchr((char*)locus,'_');
		ogid=atoi(locus+2);
		dup=atoi(pos+1)-1;
		cout << "Found OGID: " << ogid << " on level: " << dup << endl;
		mgrp = mor.getMappedOrthogroups().find(ogid)->second;
		map<int,map<string,string>*> geneSets=mgrp->getGeneSets();
		if(geneSets.find(dup)!=geneSets.end())
		{
                	map<string,string>* levelSet=geneSets.find(dup)->second;
                	if(levelSet!=NULL && levelSet->find(celltype)!=levelSet->end())
                	{
                		name=levelSet->find(celltype)->second.c_str();
				cout << "Found " << name << " for " << celltype << endl;
                	}
		}
		else
		{
			cout << "Gene set for this level " << dup << " in " << ogid << " not found; something wrong with OGIDS input. Recommend to check and restart!" << endl;
			//may add exit here if necessary.
		}
	}
	return name;
}

int
main(int argc, const char** argv)
{
	if(argc!=10)	//JSmod
	{
		//This is configured for cmint and expects specific files. So make sure to have those files in the datadir
		cout <<"Usage: ./findTransitionGenesets cmint_result_dir celltypeorder OGID_file srccelltype threshold outputdir maxgenesetsize plots_species_order maxmissing" << endl;	//JSmod
		return 0;
	}
	Framework fw;
	fw.readDataMatrix(argv[1],atoi(argv[9]));	//JSmod
	fw.readSpeciesOrder(argv[8]);
	fw.readOGIDs(argv[2],argv[3]);
	fw.setSrcCellType(argv[4]);
	char outputDirCmdF[1024];			//JSmod
	sprintf(outputDirCmdF,"mkdir -p  %s",argv[6]);	//JSmod
	system(outputDirCmdF);				//JSmod
	fw.generateTransitioningGeneSets(atof(argv[5]),argv[6],atoi(argv[7]));
	fw.generateOrderedClusterMeans(argv[6]);
	return 0;
}

