#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include"pathalg.h"
static const int WORK_SIZE =258;
parallelpush::parallelpush()
{
	cout<<"fuck c++,rubish!!!!"<<endl;
};
void parallelpush::init(vector<edge>&extenedges,vector<vector<int>>&relate,ginfo ginf){
	cout<<"in cuda init"<<endl;
	nodenum=ginf.enodesize;
	pnodesize=ginf.pnodesize;
	edges=extenedges;
	cout<<"out cuda init"<<endl;
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
	cout<<"max is: "<<max<<endl;
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

	cudaMalloc((void**)&dev_h,LY*W*pnodesize*sizeof(int));
	cudaMalloc((void**)&dev_mark,sizeof(int));
	cudaMalloc((void**)&dev_v,LY*W*pnodesize*sizeof(int));
	cudaMalloc((void**)&dev_neie,max*pnodesize*sizeof(int));
	cudaMalloc((void**)&dev_esign,LY*edges.size()*sizeof(int));
	cudaMalloc((void**)&dev_emark,LY*edges.size()*sizeof(int));
	cudaMalloc((void**)&dev_st,edges.size()*sizeof(int));
	cudaMalloc((void**)&dev_te,edges.size()*sizeof(int));
	cudaMemcpy(dev_mark,mark,sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_neie,neie,max*pnodesize*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_st,st,edges.size()*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_te,te,edges.size()*sizeof(int),cudaMemcpyHostToDevice);
};
__global__ void push(int*dev_h,int*dev_v,int* dev_esign,int* dev_emark,int*dev_neie,int*dev_nein,int N,int max,int W,int s,int t,int*mark)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int bi=i%N;
	int value=dev_v[i];
	if(i>=N*LY||value==0||bi/W==s||bi/W==t)return;
	int h=dev_h[i];
	int b=i*max;
	int minheight=INT_MAX;
	for(int j=0;j<max;j++)
	{
		int nbj=dev_nein[b+j];
		if(value>0&&nbj<INT_MAX)
		{
			int ebj=dev_neie[b+j];
			int hnbj=dev_h[nbj];
			int eid=abs(ebj)-1;
			if((ebj^dev_esign[eid])>0)
			{
				if(dev_emark[eid]>INT_MAX/2&&h==hnbj+1)
				{
					dev_emark[eid]=(ebj>0)?nbj:i;
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
__global__ void push1(int*dev_h,int*dev_v,int* dev_esign,int* dev_emark,int*dev_neie,int*dev_nein,int N,int max,int W,int s,int t,int*mark)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int bi=i%N;
	int value=dev_v[i];
	*mark=1;
	if(i>=N*LY||value==0||bi/W==s||bi/W==t)return;
	int h=dev_h[i];
	int b=i*max;
	int minheight=INT_MAX;
	for(int j=0;j<max;j++)
	{
		int nbj=dev_nein[b+j];
		if(value>0&&nbj<INT_MAX)
		{
			int ebj=dev_neie[b+j];
			int hnbj=dev_h[nbj];
			int eid=abs(ebj)-1;
			if((ebj^dev_esign[eid])>0)
			{
				if(dev_emark[eid]==0&&h==hnbj+1)
				{
					dev_emark[eid]++;
					value--;
					*mark=1;
				}
				minheight=(minheight<hnbj)?minheight:hnbj;
			}
		}
		else
			break;
	}
	if(value>0&&minheight<INT_MAX){dev_h[i]=minheight+1;*mark=1;}
};
__global__ void push2(int*dev_h,int*dev_v,int* dev_esign,int* dev_emark,int*st,int*te,int*neie,int N,int W,int E,int s,int t,int*mark,int max)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int bi=i%N;
	int value=dev_v[i];
	int node=bi/W;
	if(i>=N*LY||value==0||node==s||node==t)return;
	int ly=i/N;
	int off=i%W;
	int h=dev_h[i];
	int b=node*max;
	int minheight=INT_MAX;
	int ebj,nbj,hnbj,eid,seid;
	for(int j=0;j<max;j++)
	{
		ebj=neie[b+j];
		if(value>0&&ebj<INT_MAX){
			seid=abs(ebj)-1;
			eid=ly*E+seid;
			nbj=-1;
			bool btest=(ebj^dev_esign[eid])>0;
			//bool b1=ebj>0&&dev_esign[eid]>0;
			//bool b2=ebj<0&&dev_esign[eid]<0&&abs(dev_esign[eid])==off;
			if(btest)
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
		else
		{
			t=dev_st[bi];
			s=dev_te[bi]+1;
		}
		s+=bottom;
		t+=bottom;
		for(int k=0;k<W;k++)
			{
				int h1=dev_h[s+k];
				int h2=dev_h[t+k];
				if(dev_v[s+k]>0&&h1==h2+1)
				{
					atomicSub(&dev_v[s+k],1);
					atomicAdd(&dev_v[t+k],1);
					//dev_esign[i]=(dev_esign[i]>0)?-(k+t)%W:1;
					dev_esign[i]*=-1;
					break;
				}
			}
	}

	dev_emark[i]=0;
};
__global__ void pushrelable(int*dev_h,int*dev_v,int* dev_esign,int* dev_emark,int*dev_neie,int*dev_nein,int N,int max,int W,int s,int t,int*mark)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int bi=i%N;
	int value=dev_v[i];
	if(i>=N*LY||value==0||bi/W==s||bi/W==t)return;
	int h=dev_h[i];
	int b=i*max;
	int minheight=INT_MAX;
	for(int j=0;j<max;j++)
	{
		int nbj=dev_nein[b+j];
		if(value>0&&nbj<INT_MAX)
		{
			int ebj=dev_neie[b+j];
			int hnbj=dev_h[nbj];
			int eid=abs(ebj)-1;
			if((ebj^dev_esign[eid])>0)
			{
				if(dev_emark[eid]==i)
				{
					dev_emark[eid]++;
					atomicSub(&dev_v[i],1);
					atomicAdd(&dev_v[nbj],1);
					value--;
					dev_esign[eid]*=-1;
					*mark=1;
				}
				minheight=min(minheight,hnbj);
			}
		}
		else
			break;
	}
	if(value>0&&minheight<INT_MAX/2){dev_h[i]=minheight+1;*mark=1;}
};

__global__ void aggregate4(int* dev_esign,int* dev_v,int* dev_emark,int* dev_st,int* dev_te,int*dev_h,int W,int E)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	if(i>=E*LY)return;
	int s,t;
	dev_emark[i]=INT_MAX;
	if(dev_esign[i]>0)
	{
		s=dev_st[i];
		t=dev_te[i]+1;
	}
	else
	{
		t=dev_st[i];
		s=dev_te[i]+1;
	}
	for(int k=0;k<W;k++)
		{
			int h1=dev_h[s+k];
			int h2=dev_h[t+k];
			if(dev_v[s+k]>0&&h1==h2+1)
			{
				dev_emark[i]=s+k;
				break;
			}
		}
};
__global__ void aggregate2(int* dev_esign,int*dev_v,int* dev_emark,int W,int E)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	if(i>=E*LY)return;
	int emid=dev_emark[i];
	if(emid<INT_MAX)
	{
		int s=abs(dev_esign[i])-2+emid%W;
		if(dev_esign[i]>0)
		{	atomicSub(&dev_v[s],1);
			atomicAdd(&dev_v[emid],1);
		}
		else
		{	atomicAdd(&dev_v[s],1);
			atomicSub(&dev_v[emid],1);
		}
		dev_esign[i]*=-1;
	}
	dev_emark[i]=INT_MAX;
};
__global__ void aggregate5(int* dev_esign,int* dev_v,int* dev_emark,int* dev_st,int* dev_te,int*dev_h,int W,int E)
{
        int i = threadIdx.x + blockIdx.x*blockDim.x;
        if(i>=E*LY*W)return;
        int s,t;
        int eid=i/W;
        int k=i%W;
        if(dev_esign[eid]>0)
        {
                s=dev_st[eid];
                t=dev_te[eid]+1;
        }
        else
        {
                t=dev_st[eid];
                s=dev_te[eid]+1;
        }
        int h1=dev_h[s+k];
        int h2=dev_h[t+k];
        if(dev_v[s+k]>0&&h1==h2+1)
                dev_emark[eid]=s+k;
};

__global__ void relable(int*dev_h,int*dev_v,int N,int*mark,int*dev_nein,int*dev_neie,int *dev_esign,int max,int W,int s,int t)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int bi=i%N;
	if(i>=N*LY||dev_v[i]==0||bi/W==s||bi/W==t)return;
	int b=i*max;
	int mini=INT_MAX;
	for(int j=0;j<max;j++)
	{
		if(dev_nein[b+j]<INT_MAX)
		{
			if((dev_neie[b+j]^dev_esign[abs(dev_neie[b+j])-1])>0)
				mini=min(mini,dev_h[dev_nein[b+j]]);
		}
		else
			break;
	}
	if(mini!=INT_MAX)
		dev_h[i]=mini+1,*mark=1;
};

pair<int,int> parallelpush::prepush(int s,int t,int bw)
{
	cout<<"**********************************"<<endl;
	cout<<"parral: "<<LY<<" "<<pnodesize<<" "<<s<<" "<<t<<endl;
	time_t start,end;
	start=clock();
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
	for(int k=0;k<LY;k++)
		{
		for(int i=0;i<edges.size();i++)
			if(edges[i].s==s)
				{
				v[k*W*pnodesize+W*edges[i].t+1]=1;
				esign[k*edges.size()+i]*=-1;
				}
		}
	for(int k=0;k<LY;k++)
		{
		int start=k*W*pnodesize;
		for(int i=s*W;i<s*W+W;i++)
			h[i+start]=WD;
		}
	cudaMemcpy(dev_h,h,LY*W*pnodesize*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_v,v,LY*W*pnodesize*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_esign,esign,LY*edges.size()*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_emark,emark,LY*edges.size()*sizeof(int),cudaMemcpyHostToDevice);
	*mark=1;
	int time=0;
	cout<<"max is "<<max<<endl;
	while(*mark!=0)
	{
		if(time%100==0)
			{*mark=0;
			cudaMemcpy(dev_mark,mark,sizeof(int),cudaMemcpyHostToDevice);}
		push2<<<LY*W*pnodesize/WORK_SIZE+1,WORK_SIZE>>>(dev_h,dev_v,dev_esign,dev_emark,dev_st,dev_te,dev_neie,W*pnodesize,W,edges.size(),s,t,dev_mark,max);
		//push1<<<LY*W*pnodesize/WORK_SIZE+1,WORK_SIZE>>>(dev_h,dev_v,dev_esign,dev_emark,dev_neie,dev_nein,W*pnodesize,max,W,s,t,dev_mark);
		//aggregate2<<<LY*edges.size()/WORK_SIZE+1,WORK_SIZE>>>(dev_esign,dev_v,dev_emark,W,edges.size());
		aggregate3<<<LY*edges.size()/WORK_SIZE+1,WORK_SIZE>>>(dev_esign,dev_v,dev_emark,dev_st,dev_te,dev_h,W,edges.size(),W*pnodesize);
		/*cudaMemcpy(emark,dev_emark,LY*edges.size()*sizeof(int),cudaMemcpyDeviceToHost);
		for(int i=0;i<LY*edges.size();i++)
			if(emark[i]>0)
				cout<<"gota... "<<i<<"s:"<<st[i]<<" "<<te[i]<<" "<<emark[i]<<endl;*/
		//relable<<<LY*W*pnodesize/WORK_SIZE+1,WORK_SIZE>>>(dev_h,dev_v,W*pnodesize,dev_mark,dev_nein,dev_neie,dev_esign,max,W,s,t);
		//aggregate2<<<LY*edges.size()/WORK_SIZE+1,WORK_SIZE>>>(dev_esign,dev_v,dev_emark,W,edges.size(),W*pnodesize,dev_mark);
		if(time%100==0)
			cudaMemcpy(mark,dev_mark,sizeof(int),cudaMemcpyDeviceToHost);
		/*cudaMemcpy(v,dev_v,LY*W*pnodesize*sizeof(int),cudaMemcpyDeviceToHost);
		cudaMemcpy(h,dev_h,LY*W*pnodesize*sizeof(int),cudaMemcpyDeviceToHost);
		cudaMemcpy(esign,dev_esign,LY*edges.size()*sizeof(int),cudaMemcpyDeviceToHost);
		int flow=0;
		cout<<"************* "<<time<<endl;
		for(int i=0;i<LY*W*pnodesize;i++)
			if(v[i]!=0)
				{
					int bi=i%(W*pnodesize);
					if(bi/W==t)flow+=v[i];
					cout<<i/(W*pnodesize)<<"\t"<<bi<<"\t"<<bi/W<<"\t"<<bi%W<<"\t"<<h[i]<<"\t"<<v[i]<<endl;
					/*if(i==319)
					{
						for(int j=0;j<max;j++)
							if(nein[i*max+j]<INT_MAX)
								cout<<neie[i*max+j]<<" "<<esign[abs(neie[i*max+j])-1]<<" "<<h[nein[i*max+j]]<<endl;
					}*/
				//}
		//cout<<"mark "<<*mark<<endl;
		time++;
	}
	cudaMemcpy(mark,dev_mark,sizeof(int),cudaMemcpyDeviceToHost);
	end=clock();
	cout<<"GPU time is: "<<end-start<<endl;
	cudaMemcpy(v,dev_v,LY*W*pnodesize*sizeof(int),cudaMemcpyDeviceToHost);
	cudaMemcpy(h,dev_h,LY*W*pnodesize*sizeof(int),cudaMemcpyDeviceToHost);
	int flow=0;
	for(int i=0;i<LY*W*pnodesize;i++)
		if(v[i]!=0)
			{
				int bi=i%(W*pnodesize);
				if(bi/W==t)flow+=v[i];
				//cout<<i/(W*pnodesize)<<" "<<bi<<" "<<bi/W<<" "<<bi%W<<" "<<h[i]<<" "<<v[i]<<endl;
			}
	cudaMemcpy(esign,dev_esign,LY*edges.size()*sizeof(int),cudaMemcpyDeviceToHost);
	int count=0;
	for(int i=0;i<edges.size()*LY;i++)
		if(esign[i]<0)
			count++;
	cout<<"resort"<<endl;
	/*for(int i=0;i<edges.size();i++)
		{
			if(esign[i]<0)
			{
				int sorce=edges[i].t*W;
				if(sorce/W==t)
				{
					int pre=edges[i].s*W;
					cout<<pre<<" ";
					while((pre/W)!=s)
					{
						int flag=0;
						for(int h=0;h<W;h++)
						{
							pre++;
							for(int k=0;k<max;k++)
								{
									if(nein[pre*max+k]<INT_MAX)
										if(esign[abs(neie[pre*max+k])-1]<0&&neie[pre*max+k]<0)
										{
											esign[abs(neie[pre*max+k])-1]*=-1;
											pre=edges[abs(neie[pre*max+k])-1].s*W;
											cout<<pre<<" ";
											flag=1;
										}
										if(flag==1)break;
								}
							if(flag==1)break;
						}
					}
					cout<<endl;
				}
			}
		}*/
	cout<<"flow is"<<flow<<endl;
	cout<<"count is "<<count<<endl;
	cout<<"die is "<<time<<endl;
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
	cudaFree(dev_h);
	cudaFree(dev_mark);
	cudaFree(dev_v);
	cudaFree(dev_neie);
	cudaFree(dev_nein);
	cudaFree(dev_esign);
	cudaFree(dev_emark);
}
parallelpush:: ~parallelpush(){dellocate();};



