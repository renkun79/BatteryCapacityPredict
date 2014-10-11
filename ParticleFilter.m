function [s,q,L]=ParticleFilter(s_Forward,z,F,t,model)
%%参数说明
%      s_Forward--前一时刻的粒子
%      z--测量值
%      F--状态转移方程
%      N--粒子的个数
%      GNoise--高斯白噪声
%      s_Pre--粒子
%      z_Pre--每个粒子的测量
%      o--方差向量
%      s--该时刻最终的粒子
%      q--粒子权值
%      L--该模型的Likelihood
%      s_average--未加权的平均
%%
N=length(s_Forward);
s_Pre=zeros(1,N);
z_Pre=zeros(1,N);
s=zeros(1,N);
q=zeros(1,N);
z_temp=zeros(1,N);
sigma=sqrt(0.0001);
%%
%预测
for i=1:N
    GNoise=randn(1);
    s_Pre(i)=F(s_Forward(i),t,model)+sigma*GNoise;
end;
%%
%测量
for i=1:N
    GNoise=randn(1);
    z_Pre(1,i)=s_Pre(1,i)+sigma*GNoise;
end;
%%
%重要性函数
for i=1:N
    d=z_Pre(i)-z;
    q(1,i)=normpdf(d,0,sigma);
end;
sumq=sum(q(:));
%归一化
for i=1:N
    q(1,i)=q(1,i)/sumq;
end;
%%
%计算状态估计的似然值
s_average=0;
for i=1:N
    s_average=s_average+q(1,i)*s_Pre(i);
end;
for i=1:N
    d=s_Pre(i)-s_average;
    q(1,i)=q(1,i)*normpdf(d,0,sigma);
end;
sumq=sum(q(:));
%归一化
for i=1:N
    q(1,i)=q(1,i)/sumq;
end;
s=s_Pre;
%%
%重采样
for i=1:N
    u=rand;
    qtempsum=0;
    for j=1:N
        qtempsum=qtempsum+q(j);
        if (qtempsum>u)
            s(1,i)=s_Pre(1,j);
            z_temp(1,i)=z_Pre(1,j);
            break;
        end;
    end;
end; 
z_Pre=z_temp;
%%
%求模型的likelihood   
z_average=0;
for i=1:N
    z_average=z_average+z_Pre(1,i).*q(i);
end;
%残差
rj=0;
for i=1:N
    rj=rj+abs(z_Pre(1,i)-z);
end;
rj=rj./N;
%L=gaussmf(rj,[sigma 0]);
L=normpdf(rj,0,sigma);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





























