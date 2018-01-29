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
    	{	ERGraph graph(50,1,d2,d1);
    		cout<<"graph finied"<<endl;
    		double s=0,p=0;
    		//for(int i=0;i<20;i++)
    		graph.prepush(30,18,100,outfile);
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
