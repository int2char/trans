#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include"pathalg.h"
static const int WORK_SIZE =258;
parallelpush::parallelpush()
{
};
void parallelpush::init(vector<edge>&extenedges,vector<vector<int>>&relate,ginfo ginf){
	nodenum=ginf.enodesize;
	pnodesize=ginf.pnodesize;
	edges=extenedges;
	W=WD+1;
	h=new int[W*pnodesize*LY];
	v=new int[W*pnodesize*LY];
	mark=new int;
	vector<vector<int>>rawneie(pnodesize,vector<int>());
	for(int i=0;i<edges.size();i++)
		{
			int s=edges[i].s;
			int t=edges[i].t;
			rawneie[s].push_back(i+1);
			rawneie[t].push_back(-(i+1));
		}
	max=0;
	for(int i=0;i<rawneie.size();i++)
		if(rawneie[i].size()>max)max=rawneie[i].size();
	max++;
	neie=new int[pnodesize*max];
	for(int i=0;i<pnodesize;i++)
		{
			for(int j=0;j<max;j++)
			{
				if(j<rawneie[i].size())
					neie[i*max+j]=rawneie[i][j];
				else
					neie[i*max+j]=INT_MAX;
			}
		}
	emark=new int[LY*edges.size()];
	esign=new int[LY*edges.size()];
	st=new int[edges.size()];
	te=new int[edges.size()];
	for(int i=0;i<edges.size();i++)
		{
			st[i]=edges[i].s*W;
			te[i]=edges[i].t*W;
		}
	source=new int[pnodesize];
	ends=new int[pnodesize];
	cudaMalloc((void**)&dev_h,LY*W*pnodesize*sizeof(int));
	cudaMalloc((void**)&dev_mark,sizeof(int));
	cudaMalloc((void**)&dev_v,LY*W*pnodesize*sizeof(int));
	cudaMalloc((void**)&dev_neie,max*pnodesize*sizeof(int));
	cudaMalloc((void**)&dev_esign,LY*edges.size()*sizeof(int));
	cudaMalloc((void**)&dev_emark,LY*edges.size()*sizeof(int));
	cudaMalloc((void**)&dev_st,edges.size()*sizeof(int));
	cudaMalloc((void**)&dev_ends,pnodesize*sizeof(int));
	cudaMalloc((void**)&dev_source,pnodesize*sizeof(int));
	cudaMalloc((void**)&dev_te,edges.size()*sizeof(int));
	cudaMemcpy(dev_mark,mark,sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_neie,neie,max*pnodesize*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_st,st,edges.size()*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_te,te,edges.size()*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_ends,ends,pnodesize*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_source,source,pnodesize*sizeof(int),cudaMemcpyHostToDevice);
};
__global__ void push2(int*dev_h,int*dev_v,int* dev_esign,int* dev_emark,int*st,int*te,int*neie,int N,int W,int E,int*mark,int max,int*dev_source,int*dev_ends)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int bi=i%N;
	int value=dev_v[i];
	int node=bi/W;
	if(i>=N*LY||value==0||dev_ends[node]==1||dev_h[i]>2*W+1)return;
	if(bi%W==0&&dev_h[i]>W+1&&dev_source[node]==1){dev_v[i]=0;return;}
	int ly=i/N;
	int off=i%W;
	int h=dev_h[i];
	int b=node*max;
	int minheight=INT_MAX;
	int ebj,nbj,hnbj,eid,seid;
	int flag=0;
	for(int j=0;j<max;j++)
	{
		ebj=neie[b+j];
		if(ebj<INT_MAX&&value>0){
			seid=abs(ebj)-1;
			eid=ly*E+seid;
			nbj=-1;
			bool b1=ebj>0&&dev_esign[eid]>0;
			bool b2=ebj<0&&dev_esign[eid]<0&&(abs(dev_esign[eid])==off);
			if(b1||b2)
			{
				if(ebj>0&&off<W-1)
					nbj=te[seid]+off+1;
				if(ebj<0&&off>0)
					nbj=st[seid]+off-1;
				if(nbj<0)continue;
				nbj+=ly*N;
				hnbj=dev_h[nbj];
				if(dev_emark[eid]==0&&h==hnbj+1)
				{
					dev_emark[eid]++;
					value--;
					*mark=1;
				}
				minheight=min(minheight,hnbj);
			}
		}
		else
			break;
	}
	if(value>0&&minheight<INT_MAX){dev_h[i]=minheight+1;*mark=1;}
};
__global__ void aggregate3(int* dev_esign,int* dev_v,int* dev_emark,int* dev_st,int* dev_te,int*dev_h,int W,int E,int N)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	if(i>=E*LY)return;
		int s,t;
	int bottom=(i/E)*N;
	int bi=i%E;
	if(dev_emark[i]>0)
	{
		if(dev_esign[i]>0)
		{
			s=dev_st[bi];
			t=dev_te[bi]+1;
		}
		if(dev_esign[i]<0)
		{
			t=dev_st[bi];
			s=dev_te[bi]+1;
		}
		s+=bottom;
		t+=bottom;
		for(int k=0;k<W-1;k++)
			{
				int h1=dev_h[s+k];
				int h2=dev_h[t+k];
				if(dev_v[s+k]>0&&h1==h2+1)
				{
					atomicSub(&dev_v[s+k],1);
					atomicAdd(&dev_v[t+k],1);
					dev_esign[i]=(dev_esign[i]>0)?-(k+t)%W:1;
					break;
				}
			}
			
	}
	dev_emark[i]=0;
};
pair<int,int> parallelpush::prepush(int slen,int tlen,int bw)
{
	for(int i=0;i<LY*edges.size();i++)
		emark[i]=0;
	for(int k=0;k<LY;k++)
		for(int i=0;i<edges.size();i++)
			esign[i+k*edges.size()]=1;
	for(int i=0;i<W*LY*pnodesize;i++)
		{
			h[i]=0;
			v[i]=0;
		}
	for(int i=0;i<pnodesize;i++)
	{
		source[i]=0;
		ends[i]=0;
	}
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
						v[k*nodenum+j*W]=1;
						}
				}
			}
	}
	for(int i=0;i<pnodesize;i++)
	{
		while(tlen>0)
			{
			int j=rand()%pnodesize;
			if(source[j]==0&&ends[j]==0)
				{
					ends[j]=1;
					tlen--;
				}
			}
	}
	for(int i=0;i<LY*edges.size();i++)
	{
		int ran=rand()%100;
		if(ran<50)
			esign[i]=0;
	}
	cudaMemcpy(dev_h,h,LY*W*pnodesize*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v,v,LY*W*pnodesize*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_esign,esign,LY*edges.size()*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_emark,emark,LY*edges.size()*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_ends,ends,pnodesize*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_source,source,pnodesize*sizeof(int),cudaMemcpyHostToDevice);
	*mark=1;
	int time=0;
	time_t start,end;
	start=clock();
	int flag=0;
	int fl2=1;
	while(*mark!=0)
		{
			if(time%20==0)
				{	*mark=0;
					cudaMemcpy(dev_mark,mark,sizeof(int),cudaMemcpyHostToDevice);
				}
			
			push2<<<LY*nodenum/WORK_SIZE+1,WORK_SIZE>>>(dev_h,dev_v,dev_esign,dev_emark,dev_st,dev_te,dev_neie,nodenum,W,edges.size(),dev_mark,max,dev_source,dev_ends);
			cudaMemcpy(emark,dev_emark,LY*edges.size()*sizeof(int),cudaMemcpyDeviceToHost);
			aggregate3<<<LY*edges.size()/WORK_SIZE+1,WORK_SIZE>>>(dev_esign,dev_v,dev_emark,dev_st,dev_te,dev_h,W,edges.size(),W*pnodesize);
			
			if(time%20==0)
				cudaMemcpy(mark,dev_mark,sizeof(int),cudaMemcpyDeviceToHost);
			time++;
	    }	
	int flow=0;
	cudaMemcpy(v,dev_v,LY*nodenum*sizeof(int),cudaMemcpyDeviceToHost);
	end=clock();
	for(int i=0;i<W*pnodesize*LY;i++)
			if(v[i]!=0)
				{
					int bi=i%(W*pnodesize);
					if(ends[bi/W]==1)flow+=v[i];
	}
	cout<<"GPU flow is "<<flow<<endl;
	cudaMemcpy(esign,dev_esign,LY*edges.size()*sizeof(int),cudaMemcpyDeviceToHost);
	vector<int>vesign,vvalue,vends;
	for(int i=0;i>edges.size()*LY;i++)
		vesign.push_back(esign[i]);
	for(int i=0;i<nodenum*LY;i++)
		vvalue.push_back(v[i]);
	for(int i=0;i<pnodesize;i++)
		vends.push_back(ends[i]);
	//dilor->checkhop(0,0,vesign,vvalue,vends);
    cout<<"GPU time is: "<<end-start<<endl;
	return make_pair(flow,end-start);
};

void parallelpush:: dellocate()
{
	/*delete[] h;
	delete[] minarray;
	delete[] v;
	delete[] mark;
	delete[] neie;
	delete[] nein;
	delete[]emark;
	delete[]esign;*/
	/*cudaFree(dev_h);
	cudaFree(dev_mark);
	cudaFree(dev_v);
	cudaFree(dev_neie);
	cudaFree(dev_nein);
	cudaFree(dev_esign);
	cudaFree(dev_emark);*/
}
parallelpush:: ~parallelpush(){};



