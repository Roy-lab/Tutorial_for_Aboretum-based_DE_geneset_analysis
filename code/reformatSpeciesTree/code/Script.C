#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <stdlib.h>
#include <algorithm>
#include <cstdlib>
#include "GeneTree.H"
#include "NewickReader.H"

using namespace std;

GeneTree*
readTreeFromFile(const char* treeFName)
{
        NewickReader nreader;
        GeneTree* gtree=nreader.readTree(treeFName);
        return gtree;
}


int
showTree(GeneTree* anode)
{
        GeneTree* lchild=anode->leftchild;
        if(lchild!=NULL)
        {
                if(anode->nodeType==1)
                {
                        cout << lchild->species << "\tleft\t" << anode->species << endl;
                }
                else if(anode->nodeType==2)
                {
                        cout << "(Duplicate) " << lchild->species << "\tleft\t" << anode->species << endl;
                }
                else
                {
                        cout <<"Unrecognized node type " << anode->nodeType << endl;
                }
        }
        GeneTree* rchild=anode->rightchild;
        if(rchild!=NULL)
        {
                if(anode->nodeType==1)
                {
                        cout << rchild->species << "\tright\t" << anode->species << endl;
                }
                else if(anode->nodeType==2)
                {
                        cout << "(Duplicate) " << rchild->species << "\tright\t" << anode->species << endl;
                }
                else
                {
                        cout <<"Unrecognized node type " << anode->nodeType << endl;
                }
        }
        if(lchild!=NULL)
        {
                showTree(lchild);
        }
        if(rchild!=NULL)
        {
                showTree(rchild);
        }
        return 0;
}

int
main(int argc, const char** argv)
{
        if(argc!=2)
        {
                cout <<"Usage: reformatSpeciesTree speciestree" << endl;
                return 0;
        } 
        GeneTree* tree=readTreeFromFile(argv[1]);
        showTree(tree);
	return 0;
}
