#include <iostream>
#include"Graph.h"
#include<sys/time.h>
int main()
{
	ofstream outfile;
	outfile.open("data.txt", ios::app);
    parallelpush d1=parallelpush();
    dijkstor d2=dijkstor();
    double nflow=0;
    double nmxflow=0;
    int count=0;
    ERGraph graph(400,1,d2,d1);
    graph.prepush(0,398,400,outfile);
    for(int n=100;n<5000;n*=2)
    	{	ERGraph graph(n,1,d2,d1);
    	   	double flow=0;
    	    double mxflow=0;
    		for(int i=0;i<40;i++)
    		{
    			int s=rand()%n;
    			int t=s;
    			while(t==s)
    			 t=rand()%n;
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
    		/*cout<<"graph finied"<<endl;
    		double s=0,p=0;
    		int min=INT_MAX;
    		//for(int i=0;i<20;i++)
    		int bs=0,bt=0;
    		for(int j=0;j<10;j++)
    			for(int i=0;i<10;i++)
    				{
    				int out=1000000;
    				if(i!=j)
    					out=graph.prepush(i,j,100,outfile).first;
    				if(out<min)
    					{
    						bs=i;
    						bt=j;
    						min=out;
    					}
    				}
    		cout<<"all out "<<endl;
    		cout<<bs<<" "<<bt<<" "<<min<<endl;
    		/*for(int i=0;i<20;i++)
    			{
    				int s=rand()%n;
    				int t=s;
    				while(t==s)
    					t=rand()%n;
    				pair<int,int> a=graph.prepush(s,t,n,outfile);
    				s+=a.first;
    				p+=a.second;
    			}
    		outfile<<LY<<" "<<n<<" "<<s/20<<" "<<p/20<<endl;*/
    	}
    	cout<<"tt average flow is "<<(double)nflow/(double)count<<endl;
        cout<<"tt average ratio is"<<(double)nflow/(double)nmxflow;
}
