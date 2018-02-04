#include <iostream>
#include"Graph.h"
#include<sys/time.h>
int main()
{
	ofstream outfile;
	outfile.open("data.txt", ios::app);
    dijkstor d2=dijkstor();
    parallelpush d1=parallelpush();
    d1.fuzhi(d2);
    double nflow=0;
    double nmxflow=0;
    int count=0;
    ERGraph graph(500,1,d2,d1);
    //12,1,s=42//
    graph.prepush(20,20,400,outfile);
   /* for(int n=100;n<5000;n*=2)
    	{	ERGraph graph(n,1,d2,d1);
    	   	double flow=0;
    	    double mxflow=0;
    		for(int i=0;i<40;i++)
    		{
    			int s=rand()%n;
    			int t=s;
    			while(t==s)t=rand()%n;
    			pair<int,int>dd=graph.prepush(s,t,100,outfile);
    			flow+=dd.first;
    			mxflow+=dd.second;
        		count++;
    		}
    		nflow+=flow;
    		nmxflow+=mxflow;
    		cout<<"****************************************"<<endl;
    		cout<<n<<"average flow is "<<(double)flow/(double)count<<endl;
    		cout<<n<<"average ratio is"<<(double)flow/(double)mxflow<<endl;
    	}
    	cout<<"tt average flow is "<<(double)nflow/(double)count<<endl;
        cout<<"tt average ratio is"<<(double)nflow/(double)nmxflow;*/
}
