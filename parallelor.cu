#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include"pathalg.h"
static const int WORK_SIZE =258;
__global__ void BFShigh(int t,int *m,int *st,int *te,int *d,int *chan,int round,int edgesize,int nodenum)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	if(i>=edgesize)return;
	int from=st[i];
	if (chan[from]<0)return;
	chan[from]=-1;
	int to=te[i];
	d[to]=round;
	if((to%nodenum)/(WD+1)==t)*m=1;
}
__global__ void initchan(int s,int *chan,int *d,int *pred,int nodenum)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	if(i>=nodenum*LY)return;
	int bi=i%nodenum;
	int W=WD+1;
	chan[i]=(bi/W==s)?1:-1;
	d[i]=(bi/W==s)?0:inf;
	pred[i]=d[i];
}
__global__ void chanchan(int *m,int *pred,int *d,int *chan,int nodenum)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	if(i>=nodenum*LY)return;
	chan[i]=-1;
	if(d[i]<pred[i])
	{
		chan[i]=1;
		pred[i]=d[i];
	}
}
void parallelor::copydata(int s,vector<edge>&edges,int nodenum){
	memset(pre,-1,sizeof(int)*nodenum);
	*m=0;
	for(int i=0;i<nodenum;i++)
		d[i]=INT_MAX/2;
	d[s]=0;
	for(int i=0;i<edges.size();i++)
		aedges[i]=edges[i];
	cudaMemcpy(dev_edges,aedges,edges.size()* sizeof(edge),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_m,m,sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_d,d,sizeof(int)*nodenum,cudaMemcpyHostToDevice);
	cudaMemcpy(dev_pre,pre,sizeof(int)*nodenum,cudaMemcpyHostToDevice);
};
void parallelor::dellocate(){
	/*delete[]d;
	delete[]pre;
	delete[]aedges;
	delete m;
	cudaFree(dev_edges);
	cudaFree(dev_m);
	cudaFree(dev_d);
	cudaFree(dev_pre);*/
};
void parallelor::allocate(int maxn,int maxedge){
	m=new int;
	d=new int[maxn],pre=new int[maxn];
	aedges=new edge[maxedge];
	cudaMalloc(&dev_edges, sizeof(edge)*maxedge);
	cudaMalloc((void**)&dev_d,maxn*sizeof(int));
	cudaMalloc((void**)&dev_pre,maxn*sizeof(int));
	cudaMemcpy(duan,dev_duan,duansize*sizeof(int),cudaMemcpyDeviceToHost);
	cudaMalloc((void**)&dev_m,sizeof(int));
}
bool parallelor::cutcake(int index){

	cout<<"cut "<<index<<endl;
	if(maxbw-(index+1)*10>=0)
		maxbw-=(index+1)*10;
	else
		{
			cout<<"failure"<<endl;
			return false;
		}
	hleveln[index]++;
	return true;
};
void parallelor::topsort()
{
	cout<<" in top sort "<<endl;
	queue<int>zero;
	vector<int>order(nodenum*LY,-1);
	for(int i=0;i<nodenum*LY;i++)
		zero.push(i);
	int biao=0;
	while(!zero.empty())
	{
		int node=zero.front();
		zero.pop();
		order[node]=biao++;
		for(int i=0;i<neibn[node].size();i++)
		{
			if((--ancestor[neibn[node][i]])==0)
				zero.push(neibn[node][i]);
		}
	}
	vector<pair<int,int>>tmp;
	for(int i=0;i<order.size();i++)
		tmp.push_back(make_pair(i,order[i]));
	sort(tmp.begin(),tmp.end(),pairless());
	for(int i=0;i<order.size();i++)
		ordernode.push_back(tmp[i].first);
};
void parallelor::init(vector<edge>&extenedges,vector<vector<int>>&relate,ginfo ginf){
	cout<<"in cuda init"<<endl;
	nodenum=ginf.enodesize;
	edges=extenedges;
	mark=new int;
	*mark=0;
	W=WD+1;
	int *d,*dev_d,*pred,*dev_pred;
	st=new int[2*WD*edges.size()*LY];
	te=new int[2*WD*edges.size()*LY];
	chan=new int[nodenum*LY];
	d=new int[nodenum*LY];
	pred=new int[nodenum*LY];
	vector<vector<int>>nein(nodenum*LY,vector<int>());
	vector<int>as(nodenum*LY,0);
	ancestor=as;
	neibn=nein;
	cout<<"gsdfs"<<endl;
	for(int k=0;k<LY;k++)
	{
		int startn=k*nodenum;
		for(int i=0;i<edges.size();i++)
			for(int j=0;j<W-1;j++)
			{
				int s=edges[i].s*W+j+startn;
				int t=edges[i].t*W+j+1+startn;
				ancestor[t]++;
				neibn[s].push_back(t);
				neibn[t].push_back(s);
			}
	}
	cout<<"before sort "<<endl;
	topsort();
	int count=0;
	cout<<"sort out "<<endl;
	for(int i=0;i<nodenum*LY;i++)
		for(int j=0;j<neibn[ordernode[i]].size();j++)
		{
			st[count]=ordernode[i];
			te[count]=neibn[ordernode[i]][j];
			count++;
		}
	cout<<"asdasd"<<endl;
	for(int i=0;i<nodenum*LY;i++)
	{
		chan[i]=-1;
		d[i]=INT_MAX/2;
		pred[i]=d[i];
	}
	cout<<"hrerr"<<endl;
	cudaMalloc((void**)&dev_chan,nodenum*LY*sizeof(int));
	cudaMalloc((void**)&dev_st,LY*WD*edges.size()*sizeof(int));
	cudaMalloc((void**)&dev_te,LY*WD*edges.size()*sizeof(int));
	cudaMalloc((void**)&dev_d,LY*nodenum*sizeof(int));
	cudaMalloc((void**)&dev_pred,LY*nodenum*sizeof(int));
	cudaMalloc((void**)&dev_mark,sizeof(int));
	cudaMemcpy(dev_chan,chan,nodenum*LY*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_te,te,LY*WD*edges.size()*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_st,st,LY*WD*edges.size()*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_d,d,LY*nodenum*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_pred,pred,LY*nodenum*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_mark,mark,sizeof(int),cudaMemcpyHostToDevice);
	cout<<"get out"<<endl;
};

void parallelor::initprepush(vector<edge>&extenedges,vector<vector<int>>&relate,ginfo ginf){
	cout<<"in cuda init"<<endl;
	maxbw=500;
	//allocate in cuda
	nodenum=ginf.enodesize;
	edges=extenedges;
	cout<<"out cuda init"<<endl;
}
parallelor::parallelor()
{

};
vector<int> parallelor:: routalg(int s,int t,int bw)
{
	cout<<"blasting "<<endl;
	int E=2*edges.size()*WD*LY;
	int kk=1;
	for(int i=0;i<1;i++)
	{
		*mark=0;
		initchan<< <(nodenum*LY/WORK_SIZE)+1, WORK_SIZE >> >(s,dev_chan,dev_d,dev_pred,nodenum);
		cudaMemcpy(dev_m,&mark, sizeof(int), cudaMemcpyHostToDevice);
		do{
			cudaMemcpy(chan,dev_chan,nodenum*sizeof(int), cudaMemcpyDeviceToHost);
			int cc=0;
			BFShigh << <(E/WORK_SIZE)+1, WORK_SIZE >> >(t,dev_m,dev_st,dev_te,dev_d,dev_chan,kk,E,nodenum);
			chanchan<< <(nodenum*LY/WORK_SIZE)+1, WORK_SIZE >> >(dev_m,dev_pred,dev_d,dev_chan,nodenum);
			cudaMemcpy(mark, dev_m, sizeof(int), cudaMemcpyDeviceToHost);
			kk++;
		}
		while(*mark==0);
		cout<<"out here is !"<<endl;
		//cout<<"kk is: "<<kk<<endl;
	}
	cout<<"out routalg"<<endl;
	return vector<int>();
};
int fls(int x)
{
	int position;
	int i;
	if(x!=0)
		for(i=(x>>1),position=0;i!=0;++position)
			i>>=1;
	else
		position=-1;
	return pow(2,position+1);
}