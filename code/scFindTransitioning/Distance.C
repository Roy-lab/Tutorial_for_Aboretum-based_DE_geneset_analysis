#include <math.h>
#include <iostream>
#include "Distance.H"
#include <stdlib.h>
using namespace std;

Distance::Distance()
{
}

Distance::~Distance()
{
}

double
Distance::computeSymmKLDivergence(double m1, double v1, double m2, double v2)
{
	double kl1=computeKLDivergence(m1,v1,m2,v2);
	double kl2=computeKLDivergence(m2,v2,m1,v1);
	double klsim=(kl1+kl2)/2;
	return klsim;
}

double 
Distance::computeKLDivergence(double m1, double v1, double m2, double v2)
{
	double t1=log(v2/v1);
	double t2=((m1-m2)*(m1-m2))/(v2);
	double t3=v1/v2;
	double kld=(t1+t2+t3-1)/2;
	return kld;
}

double 
Distance::computeZstat(double m1,double v1,double m2, double v2,int sampleCnt)
{
	double temp1=m1-m2;
	double temp2=sqrt((v1+v2)/sampleCnt);
	double zstat=temp1/temp2;
	return zstat;
}

double
Distance::computeCC(vector<double>& v1, vector<double>& v2)
{
	double cc=0;
	double m1=0;
	double n1=0;
	for(int i=0;i<v1.size();i++)
	{
		if(v1[i]!=-100)
		{
			m1=m1+v1[i];
			n1+=1;
		}
	}
	m1=m1/n1;

	double m2=0;
	double n2=0;
	for(int i=0;i<v2.size();i++)
	{
		if(v2[i]!=-100)
		{
			m2=m2+v2[i];
			n2+=1;
		}
	}
	m2=m2/n2;	
	//cout << m1 << "\t" << n1 << "\t" << m2 << "\t" << n2 << endl;
	double xx=0;
	double yy=0;
	double xy=0;
	double oppRel=0;
	int n=0;
	for(int i=0;i<v1.size();i++)
	{
		if(v1[i]!=-100 && v2[i]!=-100)
		{
			n++;
			double diff1=v1[i]-m1;
			xx=xx+(diff1*diff1);
			double diff2=v2[i]-m2;
			yy=yy+(diff2*diff2);
			xy=xy+(diff1*diff2);
			//cout << diff1 << "\t" << diff2 << "\t" << v1[i]  << "\t" << v2[i] << endl;
			if( (diff1*diff2) < 0)
			{
				oppRel++;
			}
		}
	}
	//cout << xy << "\t" << xx << "\t" << yy << "\t" << (xy*xy)/(xx*yy) << endl;
	if(n==0)
	{
		return 0;
	}
	cc=sqrt((xy*xy)/(xx*yy));
	double threshold=n1/2.0;
	if(oppRel > threshold)
	{
		cc=cc*(-1);
	}
	return cc;
}
