function [s,q,L]=ParticleFilter(s_Forward,z,F,t,model)
%%����˵��
%      s_Forward--ǰһʱ�̵�����
%      z--����ֵ
%      F--״̬ת�Ʒ���
%      N--���ӵĸ���
%      GNoise--��˹������
%      s_Pre--����
%      z_Pre--ÿ�����ӵĲ���
%      o--��������
%      s--��ʱ�����յ�����
%      q--����Ȩֵ
%      L--��ģ�͵�Likelihood
%      s_average--δ��Ȩ��ƽ��
%%
N=length(s_Forward);
s_Pre=zeros(1,N);
z_Pre=zeros(1,N);
s=zeros(1,N);
q=zeros(1,N);
z_temp=zeros(1,N);
sigma=sqrt(0.0001);
%%
%Ԥ��
for i=1:N
    GNoise=randn(1);
    s_Pre(i)=F(s_Forward(i),t,model)+sigma*GNoise;
end;
%%
%����
for i=1:N
    GNoise=randn(1);
    z_Pre(1,i)=s_Pre(1,i)+sigma*GNoise;
end;
%%
%��Ҫ�Ժ���
for i=1:N
    d=z_Pre(i)-z;
    q(1,i)=normpdf(d,0,sigma);
end;
sumq=sum(q(:));
%��һ��
for i=1:N
    q(1,i)=q(1,i)/sumq;
end;
%%
%����״̬���Ƶ���Ȼֵ
s_average=0;
for i=1:N
    s_average=s_average+q(1,i)*s_Pre(i);
end;
for i=1:N
    d=s_Pre(i)-s_average;
    q(1,i)=q(1,i)*normpdf(d,0,sigma);
end;
sumq=sum(q(:));
%��һ��
for i=1:N
    q(1,i)=q(1,i)/sumq;
end;
s=s_Pre;
%%
%�ز���
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
%��ģ�͵�likelihood   
z_average=0;
for i=1:N
    z_average=z_average+z_Pre(1,i).*q(i);
end;
%�в�
rj=0;
for i=1:N
    rj=rj+abs(z_Pre(1,i)-z);
end;
rj=rj./N;
%L=gaussmf(rj,[sigma 0]);
L=normpdf(rj,0,sigma);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





























