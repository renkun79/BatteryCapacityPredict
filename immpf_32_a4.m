clear all;
clc;
echo off;
%%
%Interactive Multiple Model Particle Filter
%Model:Battery
%Author:Ye Zhe
%CreateDate:2014.7.9
%School:Harbin Institute of Technology,China
%%
%%%初始化%%%%%%%
N=500;       %粒子数目
f1=@F1;
f2=@F2;
f3=@F3;
predict_Start=27;
% 载入量测数据
load('a4shiyan.mat');
cycleTimes=length(a4);
x0=a4(1);  % 初始状态
%%
%绘制量测数据
figure(1)
plot(a4,'r','linewidth',2);
grid on;hold on
hold off
title('测量值');
xlabel('cycleTimes'); ylabel('capacity/%');
%%
%%%%%%%%%%%%%%%%%%%%%%%%IMMPF%%%%%%%%%%%%%%%%%%%%%
%初始化
s1=zeros(cycleTimes,N);%三个模型的粒子和对应权值
s2=zeros(cycleTimes,N);
s3=zeros(cycleTimes,N);
s=zeros(cycleTimes,N);%综合之后的粒子
q1=zeros(cycleTimes,N);
q2=zeros(cycleTimes,N);
q3=zeros(cycleTimes,N);
for i=1:N
    s1(1,i)=x0;
    s2(1,i)=x0;
    s3(1,i)=x0;
    s(1,i)=x0;
end;

x1_IMM = zeros(1,cycleTimes);      %模型1IMM算法状态估计值
x2_IMM = zeros(1,cycleTimes);      %模型2IMM算法状态估计值
x3_IMM = zeros(1,cycleTimes);      %模型3IMM算法状态估计值
x_pro_IMM = zeros(1,cycleTimes);   %IMM算法模型综合状态估计值
x_result=zeros(1,cycleTimes);
x_result(1)=x0;
P_IMM=zeros(1,cycleTimes);       %IMM算法模型综合状态协方差
P1_IMM=0;
P2_IMM=0;
P3_IMM=0;              %IMM算法各模型协方差
r1_IMM=0;
r2_IMM=0;
r3_IMM=0;
S1_IMM=0;
S2_IMM=0;
S3_IMM=0;
x_pro_IMM(:,1)=x0;

%模型转移概率矩阵
pij=[0.9,0.05,0.05;
    0.1,0.8,0.1;
    0.05,0.15,0.8];

u_IMM=zeros(3,cycleTimes);
u_IMM(:,1)=[0.5,0.3,0.2]';  %IMM算法模型概率

x1_IMM(1)=x0;x2_IMM(1)=x0;x3_IMM(1)=x0;  %IMM算法各模型初始状态

%三个模型的Likelihoods
L1=zeros(1,cycleTimes);
L2=zeros(1,cycleTimes);
L3=zeros(1,cycleTimes);
figure(2)
%%
for t=2:predict_Start
    x1_IMM(t) = 0;      %模型1IMM算法状态估计值
    x2_IMM(t) = 0;      %模型2IMM算法状态估计值
    x3_IMM(t) = 0;      %模型3IMM算法状态估计值
     %第一步Interacting（只针对IMM算法）
    c_j=pij'*u_IMM(:,t-1);
    
    %计算模型混合概率
    ui1=(1/c_j(1))*pij(:,1).*u_IMM(:,t-1);
    ui2=(1/c_j(2))*pij(:,2).*u_IMM(:,t-1);
    ui3=(1/c_j(3))*pij(:,3).*u_IMM(:,t-1);    

    %粒子滤波过程
    [s1(t,:),q1(t,:),L1(t)]=ParticleFilter(s1(t-1,:),a4(t),f1,t);
    [s2(t,:),q2(t,:),L2(t)]=ParticleFilter(s2(t-1,:),a4(t),f2,t);
    [s3(t,:),q3(t,:),L3(t)]=ParticleFilter(s3(t-1,:),a4(t),f3,t);
    
    %模型概率更新
    c=L1(t)*c_j(1)+L2(t)*c_j(2)+L3(t)*c_j(3);
    u_IMM(1,t)=(1/c)*L1(t)*c_j(1);
    u_IMM(2,t)=(1/c)*L2(t)*c_j(2);
    u_IMM(3,t)=(1/c)*L3(t)*c_j(3);
    %模型综合
    for i=1:N
        x1_IMM(t)=x1_IMM(t)+s1(t,i)*q1(t,i);
        x2_IMM(t)=x2_IMM(t)+s2(t,i)*q2(t,i);
        x3_IMM(t)=x3_IMM(t)+s3(t,i)*q3(t,i);
    end;
    x_result(t)=x1_IMM(t)*u_IMM(1,t)+x2_IMM(t)*u_IMM(2,t)+x3_IMM(t)*u_IMM(3,t);
    for i=1:N
        %s(t,i)=x_result(t);
        %%%%%%此处三个模型粒子的融合方式有问题%%%%%%%%%%%%%
        s(t,i)=s1(t,i)*q1(t,i)*u_IMM(1,t)+s2(t,i)*q2(t,i)*u_IMM(2,t)+s3(t,i)*q3(t,i)*u_IMM(3,t);
    end;
    %plot(t,x_result(t),'linewidth',2);
    plot(t,x1_IMM(t),'r*');
    plot(t,x2_IMM(t),'y*');
    plot(t,x3_IMM(t),'b*');
    grid on ;hold on
    title('模型一估计');
    xlabel('cycleTimes'); ylabel('Capacity/%');
end;
%%
for t=predict_Start+1:cycleTimes
    x1_IMM(t)=0;
    x2_IMM(t)=0;
    x3_IMM(t)=0;
    for i=1:N
        s1(t,i)=f1(s1(t-1,i),t);
        s2(t,i)=f2(s2(t-1,i),t);
        s3(t,i)=f3(s3(t-1,i),t);
        x1_IMM(t)=x1_IMM(t)+s1(t,i)*q1(predict_Start,i);
        x2_IMM(t)=x2_IMM(t)+s2(t,i)*q2(predict_Start,i);
        x3_IMM(t)=x3_IMM(t)+s3(t,i)*q3(predict_Start,i);
    end;
    %x1_IMM(t)=x1_IMM(t)/N;
    %x2_IMM(t)=x1_IMM(t)/N;
    %x3_IMM(t)=x1_IMM(t)/N;
    x_result(t)=x1_IMM(t)*u_IMM(1,predict_Start)+x2_IMM(t)*u_IMM(2,predict_Start)+x3_IMM(t)*u_IMM(3,predict_Start);
    plot(t,x1_IMM(t),'r*');
    plot(t,x2_IMM(t),'y*');
    plot(t,x3_IMM(t),'b*');
    grid on ;hold on
end;
%%
plot(2:cycleTimes,x_result(2:cycleTimes),'k','linewidth',2);
plot(2:cycleTimes,x1_IMM(2:cycleTimes),'y','linewidth',2);
plot(2:cycleTimes,x2_IMM(2:cycleTimes),'g','linewidth',2);
plot(2:cycleTimes,x3_IMM(2:cycleTimes),'b','linewidth',2);
plot(a4,'r','linewidth',2);
hold off
%%
figure
plot(a4,'r','linewidth',2);grid on;hold on
plot(2:cycleTimes,x_result(2:cycleTimes),'b','linewidth',2);
hold off
title('滤波结果与真实值的对比');
xlabel('cycleTimes'); ylabel('Capacity/%');
%%
%误差分析
%               d ------交互式的误差
%               d1------模型1的误差
%               d2------模型2的误差
%               d3------模型3的误差
d=zeros(1,cycleTimes);
d1=zeros(1,cycleTimes);
d2=zeros(1,cycleTimes);
d3=zeros(1,cycleTimes);
for t=1:cycleTimes
    d(t)=abs(x_result(t)-a4(t));
    d1(t)=abs(x1_IMM(t)-a4(t));
    d2(t)=abs(x2_IMM(t)-a4(t));
    d3(t)=abs(x3_IMM(t)-a4(t));
end;
figure
plot(d,'k','linewidth',2);grid on;
hold on
plot(d1,'y');
plot(d2,'g');
plot(d3,'b');

hold off
title('误差分析');

