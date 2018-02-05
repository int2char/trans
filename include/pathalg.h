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
#define WD 5
#ifndef LY 
	#define LY 10
#endif
#define inf INT_MAX/2
using namespace std;
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
		vector<vector<int>>ynein;
		vector<int>height;
		vector<int>value;
		vector<int>esign;
		vector<int>ordernode;
		vector<int>source;
		vector<int>ends;
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
    		
			pesize=edges.size();
			W=WD+1;
			vector<vector<int>>tneie(nodenum*LY,vector<int>());
			vector<vector<int>>tnein(nodenum*LY,vector<int>());
			vector<vector<int>>tynein(pnodesize*LY,vector<int>());
			vector<int>theight(nodenum*LY,0);
			vector<int>tvalue(nodenum*LY,0);
			vector<int>tesign(edges.size()*LY,1);
			neie=tneie;
			nein=tnein;
			ynein=tynein;
			height=theight;
			value=tvalue;
			esign=tesign;
			vector<int>tsource(pnodesize,0);
			vector<int>tends(pnodesize,0);
			source=tsource;
			ends=tends;
			//cout<<"pnodesize is :"<<pnodesize<<endl;
			for(int k=0;k<LY;k++)
			{
				int startn=k*nodenum;
				int starte=k*pesize;
				for(int i=0;i<edges.size();i++)
					{
					//cout<<edges[i].s<<" "<<edges[i].t<<endl;
					ynein[edges[i].s+k*pnodesize].push_back((i+k*edges.size()+1));
					ynein[edges[i].t+k*pnodesize].push_back(-(i+k*edges.size()+1));
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
			}
			//topsort();
        }
        virtual vector<int> routalg(int s,int t,int bw){
        	
	 	}
        static bool compare(pair<int,int>&a,pair<int,int>&b)
        {
        	if(a.second<b.second)
        		return true;
        	return false;
        }
        
        virtual pair<int,int> prepush(int slen,int tlen,int bw)
        {
        	time_t start,end;
        	cout<<endl;
        	//cout<<"**********************************"<<endl;
        	//cout<<"serial: "<<LY<<" "<<pnodesize<<" "<<s<<" "<<t<<endl;
			for(int i=0;i<LY*nodenum;i++)
			{
				height[i]=0;
				value[i]=0;
			}
			for(int i=0;i<edges.size()*LY;i++)
				esign[i]=1;
        	vector<int>ifhas(nodenum,0);
        	
        	srand(1);
        	int ccc=0;
        	for(int i=0;i<pnodesize;i++)
        	{
        		while(slen>0)
        			{
        			int j=rand()%pnodesize;
        			if(source[j]==0)
        				{
        					slen--;
        					source[j]=1;
        					for(int k=0;k<LY;k++)
        						{
        						value[k*nodenum+j*W]=1;
        						}
        				}
        			}
        	}
        	//cout<<" ccc is "<<ccc<<endl;
        	for(int i=0;i<pnodesize;i++)
        	{
        		while(tlen>0)
        			{
        			//cout<<"???"<<endl;
        			int j=rand()%pnodesize;
        			//cout<<j<<endl;
        			if(source[j]==0&&ends[j]==0)
        				{
        					ends[j]=1;
        					tlen--;
        				}
        			}
        	}
        	//cout<<"out it "<<endl;
        	for(int i=0;i<LY*edges.size();i++)
        	{
        		int ran=rand()%100;
        		if(ran<50)
        			esign[i]=0;
        	}
        	start=clock();
        	int mark=1;
        	int cc=0;
        	int flow=0;
        	int time=0;
        	int fff=INT_MAX;
        	int iasd=50;
        	while(mark==1)
        	{
        		mark=0;
        		time++;
        		for(int i=0;i<nodenum*LY;i++)
					{
        				int bi=i%nodenum;
        				//cout<<value[i]<<endl;
        				if(value[i]>0&&ends[bi/W]==0)
						{
        					//cout<<"in"<<endl;
        					int flag=0;
							int minheight=INT_MAX;
							for(int j=0;j<nein[i].size();j++)
							{
								int eid=abs(neie[i][j])-1;
								bool btest=(esign[eid]*neie[i][j]>0)?true:false;
								bool b1=(esign[eid]>0&&neie[i][j]>0)?true:false;
								bool b2=(esign[eid]<0&&neie[i][j]<0&&i%W==abs(esign[eid]))?true:false;
								int to=nein[i][j];
								if(value[i]>0&&(b1||b2)&&esign[eid]!=0)
								{
									if(height[i]>W+1&&i%W==0&&source[(i%nodenum)/W]==1){
										value[i]=0;
										continue;
									}
									{
										if(height[i]==height[to]+1)
										{
											value[i]--;
											value[to]++;
											//cout<<"push "<<i<<"to "<<to<<endl;
											if(ends[(to%nodenum)/W]==1)
												flow++;
											int off=to%W;
											esign[eid]=(esign[eid]>0)?-off:1;
											mark=1;
											flag=1;
										}
										else
											{
											minheight=min(minheight,height[to]);
											}
									}
								}
							}
							if(value[i]>0&&minheight<INT_MAX&&flag==0)
									{
										height[i]=minheight+1,mark=1;
									}
						}
    				}
        		cc++;
        	}
        	end=clock();
        	flow=0;
        	for(int i=0;i<nodenum*LY;i++)
        		if(value[i]!=0)
        			{
        				flow+=value[i];
        			}
        	cout<<"CPU flow is:"<<flow<<endl;
        	int count=0;
        	for(int i=edges.size();i<edges.size()*2;i++) 
        		if(esign[i]<0)
        		{	count++;
        			int id=i%edges.size();
        		}
        		else
        			esign[i]=0;
        	cout<<"CPU time is:"<<end-start<<endl;
        	//int maxf=checkR(s,t);
        	return make_pair(flow,end-start);
        };
    	public:
        int checkhop(int s,int t,vector<int>&esign,vector<int>&value,vector<int>&ends)
        {
        	cout<<"in check"<<endl;
        	cout<<nodenum<<endl;
        	int max=0;
        	for(int h=0;h<LY;h++)
        	{
        		cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<endl;
        		int nodeoff=h*nodenum;
        		cout<<"in"<<endl;
				for(int j=0;j<pnodesize;j++)
				{
					//cout<<"in "<<ends[j]<<endl;
					if(ends[j]>0)
					{
						//cout<<"**********************route of "<<j<<":"<<endl;
						for(int i=1;i<W;i++)
						{
							int tnode=nodeoff+j*W+i;
							while(value[tnode]>0)
							{
								int node=tnode;
								int off=i;
								max=0;
								cout<<node<<" ";
								while(off>0)
								{
									max++;
									if(max>WD+1)
										{
											//cout<<"loop occurs,erro happened!"<<endl;
											//break;
										}
									//cout<<" "<<endl;
									for(int k=0;k<nein[node].size();k++)
										{
										int eid=abs(neie[node][k])-1;
										//cout<<edges[eid].s<<endl;
										//cout<<eid<<" "<<esign[eid]<<" "<<neie[node][k]<<" "<<off<<endl;
										//cout<<node<<endl;
										if(neie[node][k]<0&&(-esign[eid]==off))
											{
											esign[eid]*=-1;
											off=off-1;
											node=edges[eid%edges.size()].s*W+off+nodeoff;
											cout<<node<<" ";
											//cout<<s<<endl;
											break;
											}
										}
								}
								cout<<endl;
								value[tnode]--;
								//cout<<endl;
							}
						}
					}
				}
        	}
        }
        int checkback(int s,int t)
        {
        	int tnode=-1;
			vector<int>visited(nodenum,0);
			vector<int>pre(nodenum,-1);
			vector<int>pree(nodenum,-1);
			int vflag=1;
			queue<int>que;
			que.push(s*W);
			visited[s*W]=1;
			while(!que.empty()&&vflag)
			{
				int node=que.front();
				que.pop();
				for(int i=0;i<nein[node].size();i++)
				{
					int eid=abs(neie[node][i])-1;
					if(esign[eid]<0)
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
			cout<<tnode%W<<endl;
			if(tnode>=0)
			{
				int prn=tnode;
				while(prn!=s*W)
				{
					cout<<prn/W<<" ";
					esign[pree[prn]]*=-1;
					prn=pre[prn];
				}
				cout<<prn/W<<" ";
			}
			cout<<endl;
			return (tnode>=0)?1:0;
        };
        int trance(int s,int t)
        {
        	int maxf=0;
        	int tnode=1;
        	while(tnode>=0)
			{
				tnode=-1;
				vector<int>visited(pnodesize,0);
				vector<int>pre(pnodesize,-1);
				vector<int>pree(pnodesize,-1);
				int vflag=1;
				queue<int>que;
				for(int i=0;i<W;i++)
					que.push(s);
				visited[s]=1;
				while(!que.empty()&&vflag)
				{
					int node=que.front();
					que.pop();
					for(int i=0;i<ynein[node].size();i++)
					{
						int ne=ynein[node][i];
						int eid=abs(ne)-1;
						if(esign[eid]<0&&ne>0)
						{	
							int to=edges[eid].t;
							if(visited[to]==0){
								pre[to]=node;
								pree[to]=eid;
								que.push(to);
								visited[to]=1;
							}
							if(to==t){tnode=to;vflag=0;break;}
						}
					}
				}
				int prn=tnode;
				int len=0;
				//cout<<"tnode is "<<tnode<<endl;
				if(tnode>=0)
				{
					int prn=tnode;
					while(prn!=s)
					{
						len++;
						esign[pree[prn]]*=-1;
						//cout<<prn<<endl;
						prn=pre[prn];
					}
					//cout<<prn<<" ";
				}
				if(len<=WD&&tnode>=0)maxf++;
			}
        	return maxf;
        }
        int checkR(int s,int t)
                {
        			cout<<"in mzx flow"<<endl;
        			for(int i=0;i<edges.size();i++)
        				esign[i]=1;
        			int tnode=10;
        			int maxf=0;
        			while(tnode>=0)
        			{
						tnode=-1;
						vector<int>visited(pnodesize,0);
						vector<int>pre(pnodesize,-1);
						vector<int>pree(pnodesize,-1);
						int vflag=1;
						queue<int>que;
						for(int i=0;i<W;i++)
							que.push(s);
						visited[s]=1;
						while(!que.empty()&&vflag)
						{
							int node=que.front();
							que.pop();
							for(int i=0;i<ynein[node].size();i++)
							{
								int ne=ynein[node][i];
								int eid=abs(ne)-1;
								if(esign[eid]*ne>0)
								{	
									int to=edges[eid].t;
									if(visited[to]==0){
										pre[to]=node;
										pree[to]=eid;
										que.push(to);
										visited[to]=1;
									}
									if(to==t){tnode=to;vflag=0;break;}
								}
							}
						}
						int prn=tnode;
						int len=0;
						if(tnode>=0)
						{
							int prn=tnode;
							while(prn!=s)
							{
								len++;
								esign[pree[prn]]*=-1;
								prn=pre[prn];
							}
							//cout<<prn<<" ";
						}
						if(len<=WD&&tnode>=0)maxf++;
        			}
        			//cout<<"max f is "<<maxf<<endl;
        			maxf=trance(s,t);
        			cout<<"max ff is "<<maxf<<endl;
        			return maxf;
         };
        
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
					cout<<prn<<" ";
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
		int*source;
		int*ends;
		int*dev_source,dev_ends;
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
		int *dev_ends;
		int *dev_source;
		int W;
		int *source;
		int *ends;
		void dellocate();
		vector<vector<int>>ynein;
		dijkstor* dilor;
	public:
		 parallelpush();
	 	 void topsort(){};
	 	 bool checkR(int s,int t);
	 	 virtual bool cutcake(int index){};
	     virtual pair<int,int> prepush(int s,int t,int bw);
	 	 virtual void init(vector<edge>&extenedges,vector<vector<int>>&relate,ginfo);
	 	 void fuzhi(dijkstor&d){dilor=&d;}
	 	 virtual vector<int> routalg(int s,int t,int bw){};
	 	 virtual ~parallelpush();
};
#endif //CSPALLOC_PATHALG_H
