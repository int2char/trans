//
// Created by root on 17-5-9.
//
#ifndef CSPALLOC_PATHALG_H
#include"limits.h"
#define CSPALLOC_PATHALG_H
#define INFCOST INT_MAX/2
#include<bits/stdc++.h>
#include <unistd.h>
#include"edge.h"
#include<sys/time.h>
#include<queue>
#define ML 50
#define BS 5
#define WD 8
#ifndef LY
	#define LY 1
#endif
#define inf INT_MAX/2
using namespace std;
class pairless {
    public:
        bool operator()(pair<int,int>&a,pair<int,int>&b){
            return a.second<b.second;
        }
};
class algbase {
    protected:
        vector<int> getrout(int &s, int &t, vector<edge> &edges, vector<int> &pre) {
            vector<int> rout;
            int pp = pre[t];
            while (pp >= 0) {
                rout.push_back(pp);
                pp = pre[edges[pp].s];
            }
            reverse(rout.begin(), rout.end());
            return rout;
        }
    public:
        algbase(){};
        virtual bool cutcake(int)=0;
        virtual void init(vector<edge>&extenedges,vector<vector<int>>&relate,ginfo)=0;
	 	virtual vector<int> routalg(int s,int t,int bw)=0;
	 	virtual pair<int,int> prepush(int s,int t,int bw)=0;
};
class dijkstor:public algbase{
    public:
		vector<vector<int>>neibour;
		vector<int>ancestor;
		int edgesize,nodenum,pesize,pnodesize,maxbw;
		vector<vector<vector<int>>>mask;
		vector<edge>edges;
		vector<int>dist;
		vector<int>pre;
		vector<int>leveln;
		vector<int>exn2n;
		vector<vector<int>>rela;
		vector<int>order;
		vector<vector<int>>neie;
		vector<vector<int>>nein;
		vector<int>height;
		vector<int>value;
		vector<int>esign;
		vector<int>ordernode;
		int W;
        dijkstor(){};
        void topsort()
        {
        	cout<<" in top sort "<<endl;
        	queue<int>zero;
        	vector<int>orz(nodenum*LY,-1);
        	order=orz;
        	for(int i=0;i<nodenum*LY;i++)
        		zero.push(i);
        	int biao=0;
			while(!zero.empty())
			{
				int node=zero.front();
				zero.pop();
				order[node]=biao++;
				for(int i=0;i<nein[node].size();i++)
				{
					if((--ancestor[nein[node][i]])==0)
						zero.push(nein[node][i]);
				}
			}
			ordernode=order;
        }
        virtual bool cutcake(int index){
        	cout<<"cut "<<index<<endl;
        	if(maxbw-(index+1)*10>=0)
        			maxbw-=(index+1)*10;
			else
				{
					cout<<"failure"<<endl;
					return false;
				}
        	vector<int>tmp;
        	for(int i=0;i<edgesize;i++)
        		tmp.push_back(i);
        	mask[index].push_back(tmp);
        	leveln[index]++;
        	return true;
        }
        virtual void init(vector<edge>&extenedges,vector<vector<int>>&relate,ginfo ginf){
        	maxbw=500;
        	rela=relate;
        	edgesize=extenedges.size(),nodenum=ginf.enodesize,pesize=ginf.pesize,pnodesize=ginf.pnodesize;
        	exn2n=ginf.exn2n;
            vector<edge> exedges;
        	edges=extenedges;
        	vector<vector<int>>nd(nodenum*LY,vector<int>());
        	neibour=nd;
        	vector<int>ad(nodenum*LY,0);
        	ancestor=ad;	
    		vector<int>bl(BS,0);
    		leveln=bl;
    		vector<vector<vector<int>>>nm(BS,vector<vector<int>>());
    		mask=nm;
    		vector<int>dd(nodenum*LY,inf);
    		dist=dd;
    		vector<int>pp(nodenum*LY,-1);
    		pre=pp;
    		//topsort();
    		
    		//prepush init
			pesize=edges.size();
			W=WD+1;
			vector<vector<int>>tneie(nodenum*LY,vector<int>());
			vector<vector<int>>tnein(nodenum*LY,vector<int>());
			vector<int>theight(nodenum*LY,0);
			vector<int>tvalue(nodenum*LY,0);
			vector<int>tesign(edges.size()*LY,1);
			neie=tneie;
			nein=tnein;
			height=theight;
			value=tvalue;
			esign=tesign;
			for(int k=0;k<LY;k++)
			{
				int startn=k*nodenum;
				int starte=k*pesize;
				for(int i=0;i<edges.size();i++)
					for(int j=0;j<W-1;j++)
					{
						int s=edges[i].s*W+j+startn;
						int t=edges[i].t*W+j+1+startn;
	        			ancestor[t]++;
						neie[s].push_back(i+1+starte);
						neie[t].push_back(-(i+1+starte));
						nein[s].push_back(t);
						nein[t].push_back(s);
					}
			}
			topsort();
        }
        virtual vector<int> routalg(int s,int t,int bw){
        		cout<<"begin to push "<<ordernode.size()<<endl;
        		int tnode=-1;
				vector<int>visited(nodenum*LY,0);
				vector<int>pre(nodenum*LY,-1);
				int vflag=1;
				priority_queue<pair<int,int>,vector<pair<int,int>>,pairless>que;
				for(int i=0;i<LY;i++)
					{	
						que.push(make_pair(s+nodenum*i,0));
	        			visited[s+nodenum*i]=1;
					}
				while(!que.empty()&&vflag)
				{
					int node=que.top().first;
					que.pop();
					cout<<node<<endl;
					for(int i=0;i<nein[node].size();i++)
					{
						if(neie[node][i]!=0)
						{	
							int to=nein[node][i];
							if(visited[to]==0){
								pre[to]=node;que.push(make_pair(to,ordernode[to]));visited[to]=1;
							}
							else
								continue;
							if((to%nodenum)/W==t){tnode=to;vflag=0;break;}
						}
					}
				}
				cout<<"finished"<<endl;
				/*vector<int>rout;
				int pp=pre[tnode];
				cout<<"tnode is "<<tnode<<endl;
				while(pp>=0)
				{
					rout.push_back(pp);
					pp=pre[pp];
				}
				cout<<endl;
				if(rout.size()>0)
				{
					for(int i=0;i<rout.size();i++)
						cout<<rout[i]<<" ";
					cout<<endl;
					return rout;
				}
				else
				{
					if(!cutcake(index))return vector<int>();
				}*/
	 	}
        static bool compare(pair<int,int>&a,pair<int,int>&b)
        {
        	if(a.second<b.second)
        		return true;
        	return false;
        }
        virtual pair<int,int> prepush(int s,int t,int bw)
        {
        	time_t start,end;
        	start=clock();
        	for(int j=0;j<neie[440].size();j++)
			{
				//int eid=abs(neie[pre][j])-1;
				//if(neie[pre][j]<0&&esign[eid]<0)
				cout<<neie[440][j]<<" ";
			}
        	cout<<endl;
        	cout<<"**********************************"<<endl;
        	cout<<"serial: "<<LY<<" "<<pnodesize<<" "<<s<<" "<<t<<endl;
			for(int i=0;i<LY*nodenum;i++)
			{
				height[i]=0;
				value[i]=0;
			}
			for(int i=0;i<edges.size()*LY;i++)
				esign[i]=1;
        	for(int k=0;k<LY;k++)
        		{
        		int startn=k*nodenum;
        		for(int i=0;i<W;i++)
        			height[startn+s*W+i]=W+1;
        		}
        	for(int k=0;k<LY;k++)
        		{
        			for(int i=0;i<edges.size();i++)
        				if(edges[i].s==s)
        					{
        					value[k*nodenum+W*edges[i].t+1]=1;
        					cout<<k*nodenum+W*edges[i].t+1<<endl;
        					esign[k*pesize+i]*=-1;
        					}
        		}
        	cout<<"gggg"<<endl;
        	int mark=1;
        	int cc=0;
        	while(mark==1)
        	{
        		mark=0;
        		for(int i=0;i<nodenum*LY;i++)
					{
        				int bi=i%nodenum;
        				if(value[i]>0&&bi/W!=s&&bi/W!=t)
						{
        					int flag=0;
							int minheight=INT_MAX;
							for(int j=0;j<nein[i].size();j++)
							{
								if(value[i]>0&&(esign[abs(neie[i][j])-1]*neie[i][j])>0)
								{
									int to=nein[i][j];
									if(height[i]==height[to]+1)
									{
										value[i]--;
										value[to]++;
										esign[abs(neie[i][j])-1]*=-1;
										mark=1;
										flag=1;
									}
									else
										minheight=min(minheight,height[to]);
								}
							}
							if(value[i]>0&&minheight<INT_MAX&&flag==0)
									height[i]=minheight+1,mark=1;
						}
    				}
        		cc++;
        	}
        	int flow=0;
        	for(int i=0;i<LY*W*pnodesize;i++)
        		if(value[i]!=0)
        			{
        				int bi=i%nodenum;
        				if(bi/W==t)flow+=value[i];
        			}
        	cout<<"flow is: "<<flow<<endl;
        	int count=0;
        	for(int i=0;i<edges.size()*LY;i++) 
        		if(esign[i]<0)
        		{	
        			cout<<i<<" "<<edges[i].s<<" "<<edges[i].t<<endl;
        			for(int k=0;k<W;k++)
        				cout<<height[edges[k].s*W+i]<<" ";
        			cout<<endl;
        			for(int k=0;k<W;k++)
        			    cout<<height[edges[k].t*W+i]<<" ";
        			cout<<endl;
        			
        		}
        		else
        			esign[i]=0;
        	while(check(s,t));
        	end=clock();
        	cout<<"CPU time is: "<<end-start<<endl;
        	//cout<<"count is: "<<count<<endl;
        	//cout<<"die is: "<<cc<<endl;
        	//return make_pair(flow,end-start);
        	return make_pair(0,0);
        }
        int check(int s,int t)
        {
			int tnode=-1;
			vector<int>visited(nodenum,0);
			vector<int>pre(nodenum,-1);
			vector<int>pree(nodenum,-1);
			int vflag=1;
			queue<int>que;
			for(int i=0;i<W;i++)
				que.push(s*W+i);
			visited[s*W]=1;
			while(!que.empty()&&vflag)
			{
				int node=que.front();
				que.pop();
				for(int i=0;i<nein[node].size();i++)
				{
					int ne=neie[node][i];
					int eid=abs(ne)-1;
					if(esign[eid]*ne<0)
					{	
						int to=nein[node][i];
						if(visited[to]==0){
							pre[to]=node;
							pree[to]=eid;
							que.push(to);
							visited[to]=1;
						}
						if(to/W==t){tnode=to;vflag=0;break;}
					}
				}
			}
			int prn=tnode;
			if(tnode>=0)
			{
				int prn=tnode;
				while(prn/W!=s)
				{
					cout<<prn/W<<" ";
					esign[pree[prn]]*=-1;
					prn=pre[prn];
				}
				cout<<prn<<" ";
			}
			cout<<endl;
			return (tnode>=0)?1:0;
        };
};
class parallelor:public algbase
{
	private:
		edge *dev_edges,*aedges;
		int*dev_m,*m,*dev_pre,*pre,*pred,*dev_pred,*dev_d,*d,*dev_mask,*mask,*dev_leveln,*leveln;
		int*dev_rela,*rela;
		int presize,dsize,masksize,levelnsize;
		int edgesize,nodenum,pesize,pnodesize;
		int neisize,duansize;
		int *choosel,*dev_choosel;
		int *rout,*dev_rout;
		int *routn,*dev_routn,*order,*dev_order;
		vector<int>hleveln,ancestor;
		int maxbw;
		int *dev_qian,*qian,*dev_qbeg,*qbeg,*dev_qsize,*qsize;
		epair *dev_nei,*nei;
		int *dev_duan,*duan;
		int *dev_beg,*beg;
		int *dev_value,*value;
		int *dev_height,*height;
		vector<vector<int>>neibour;
		vector<edge>edges;
		vector<int>ordernode;
		void allocate(int maxn,int maxedges);
		void copydata(int s,vector<edge>&edges,int nodenum);
		void dellocate();
		int W;
		int *st,*te,*dev_st,*dev_te;
		int *chan,*dev_chan;
		vector<vector<int>>neibn;
		int *mark,*dev_mark;
	public:
	 	 parallelor();
	 	 void topsort();
	 	 virtual bool cutcake(int index);
	     virtual pair<int,int> prepush(int s,int t,int bw){};
	 	 virtual void init(vector<edge>&extenedges,vector<vector<int>>&relate,ginfo);
	 	 void initprepush(vector<edge>&extenedges,vector<vector<int>>&relate,ginfo ginf);
	 	 virtual vector<int> routalg(int s,int t,int bw);
	 	 virtual ~parallelor(){dellocate();}
};
class parallelpush:public algbase
{
	private:
		vector<edge>edges;
		int nodenum;
		int pnodesize;
		int*dev_h,*h,*dev_v,*v,*dev_neie,*neie,*dev_nein,*nein;
		int*dev_esign,*esign;
		int *dev_emark,*emark,*mark,*dev_mark;
		int *minarray;
		int max;
		int *st,*te,*dev_st,*dev_te;
		int W;
		void dellocate();
	public:
		 parallelpush();
	 	 void topsort(){};
	 	 virtual bool cutcake(int index){};
	     virtual pair<int,int> prepush(int s,int t,int bw);
	 	 virtual void init(vector<edge>&extenedges,vector<vector<int>>&relate,ginfo);
	 	 virtual vector<int> routalg(int s,int t,int bw){};
	 	 virtual ~parallelpush();
};
#endif //CSPALLOC_PATHALG_H
