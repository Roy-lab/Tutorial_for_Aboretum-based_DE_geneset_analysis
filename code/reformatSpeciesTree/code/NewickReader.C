#include <string>
#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
using namespace std;
#include "GeneTree.H"
#include "NewickReader.H"

NewickReader::NewickReader()
{
	root=NULL;
}

NewickReader::~NewickReader()
{
}

GeneTree*
NewickReader::readTree(const char* aFName)
{
	ifstream inFile(aFName);
	string strbuff;
	getline(inFile,strbuff);
	int len=strbuff.length();
	char* buffer=new char[len+1];
	strcpy(buffer,strbuff.c_str());
	inFile.close();
	GeneTree* gtree=parseNewickFormat(buffer);
	delete[] buffer;
	return gtree;	
}

GeneTree*
NewickReader::parseNewickFormat(char* newickStr)
{
	int depthCnt=1;
	int index=0;
	char nodeName[256];
	char distance[256];
	int nodeindex=0;
	int distindex=0;
	GeneTree* currNode=root;
	bool popName=true;
	while(newickStr[index]!='\0')
	{
		char c=newickStr[index];
		switch (c)
		{
			case '(':
			{
				GeneTree* newNode=new GeneTree;
				char aname[25];
				sprintf(aname,"Anc%d",depthCnt);
				newNode->name.append(aname);
				newNode->species.append(aname);
				if(currNode!=NULL)
				{
					//cout<<"Found child for "<< currNode->getName() << endl;
					if((currNode->leftchild==NULL) && (strstr(currNode->name.c_str(),"Anc")!=NULL))
					{
						newNode->parent=currNode;
						currNode->leftchild=newNode;
					}
					else
					{
						GeneTree* parent=currNode->parent;
						if(parent==NULL)
						{
							cout <<"Bad format! No parent for " << currNode->name << endl;
							exit(0);
						}
						if(parent->rightchild!=NULL)
						{
							cout <<"Bad format! Right child of parent " << parent->name << " already set to " << parent->rightchild->name << endl;
							exit(0);
						}
						parent->rightchild=newNode;
						newNode->parent=parent;
					}
				}
				else
				{
					root=newNode;	
				}	
				currNode=newNode;
				depthCnt++;
				break;
			}
			case ')':
			{
				if(nodeindex>0)
				{
					nodeName[nodeindex]='\0';
					GeneTree* newChild=new GeneTree;
					newChild->name.append(nodeName);
					newChild->species.append(nodeName);
					newChild->parent=currNode->parent;
					GeneTree* parent=currNode->parent;
					parent->rightchild=newChild;
					nodeindex=0;
					distindex=0;
				}
				//cout <<"Moving up from " << currNode->getName() << " to " << currNode->getParent()->getName()<< endl;
				//Need to move up the tree
				currNode=currNode->parent;
				break;
			}
			case ',':
			{
				popName=true;
				if(nodeindex>0)
				{
					nodeName[nodeindex]='\0';
					GeneTree* newChild=new GeneTree;
					newChild->name.append(nodeName);
					newChild->species.append(nodeName);
					newChild->parent=currNode;
					currNode->leftchild=newChild;
					currNode=newChild;
				}
				
				distance[distindex]='\0';
				nodeindex=0;
				distindex=0;
				break;
			}
			case ';':
			{
				if(nodeindex>0)
				{
					nodeName[nodeindex]=0;
				}
				break;
			}
			case ':':
			{
				popName=false;	
				break;
			}
			default:
			{
				if(popName)
				{
					nodeName[nodeindex]=c;
					nodeindex++;
				}
				else
				{
					distance[distindex]=c;
					distindex++;
				}
			}
		}
		index++;
	}
	return root;
}
