#include <iostream>
#include"Graph.h"
#include<sys/time.h>
int main()
{
	ofstream outfile;
	outfile.open("data.txt", ios::app);
    parallelpush d1=parallelpush();
    dijkstor d2=dijkstor();
    for(int n=100;n<101;n*=2)
    	{	ERGraph graph(100,1,d2,d1);
    		graph.prepush(18,48,100,outfile);
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
}
