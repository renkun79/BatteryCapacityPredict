clear all;
clc;
echo off;
%%
%Interactive Multiple Model Particle Filter
%Model:Battery
%Author:Ye Zhe
%CreateDate:2014.7.9
%School:Harbin Institute of Technology,China
%说明:有4组数据a1shiyan.mat,a2shiyan.mat,a3shiyan.mat,a4shiyan.mat四组电池数据分别对应
%model的值为1,2,3,4
%%
%%%初始化%%%%%%%
N=500;       %粒子数目
f1=@F1;
f2=@F2;
f3=@F3;
% 载入量测数据
%%
%%%更换模型需要修改的参数%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('a4shiyan.mat');
initial_data=a4;%原始量测数据
model=4;
predict_Start=32;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
cycleTimes=length(initial_data);
x0=initial_data(1);  % 初始状态
%%
%绘制量测数据
figure
plot(initial_data,'r','linewidth',1);grid on;hold on
plot([0 cycleTimes],[0.72 0.72],'c');

hold off
title('测量值');
xlabel('cycleTimes'); ylabel('capacity/%');
legend('真实测量');
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

x1_IMM = zeros(cycleTimes,1);      %模型1IMM算法状态估计值
x2_IMM = zeros(cycleTimes,1);      %模型2IMM算法状态估计值
x3_IMM = zeros(cycleTimes,1);      %模型3IMM算法状态估计值
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
pij=[0.7,0.25,0.05;
    0.25,0.7,0.05;
    0.05,0.05,0.9];

u_IMM=zeros(3,cycleTimes);
u_IMM(:,1)=[0.2,0.3,0.5]';  %IMM算法模型概率

x1_IMM=x0;x2_IMM=x0;x3_IMM=x0;  %IMM算法各模型初始状态

%三个模型的Likelihoods
L1=zeros(1,cycleTimes);
L2=zeros(1,cycleTimes);
L3=zeros(1,cycleTimes);
figure
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
    [s1(t,:),q1(t,:),L1(t)]=ParticleFilter(s1(t-1,:),initial_data(t),f1,t,model);
    [s2(t,:),q2(t,:),L2(t)]=ParticleFilter(s2(t-1,:),initial_data(t),f2,t,model);
    [s3(t,:),q3(t,:),L3(t)]=ParticleFilter(s3(t-1,:),initial_data(t),f3,t,model);
    
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
        %%%%%%此处三个模型粒子的融合方式可以更改%%%%%%%%%%%%
        s(t,i)=s1(t,i)*q1(t,i)*u_IMM(1,t)+s2(t,i)*q2(t,i)*u_IMM(2,t)+s3(t,i)*q3(t,i)*u_IMM(3,t);
    end;
    %plot(t,x_result(t),'linewidth',4);
    plot(t,x1_IMM(t),'y*')%legend('模型1');
    grid on ;hold on;
    plot(t,x2_IMM(t),'g*');%legend('模型2');
    plot(t,x3_IMM(t),'b*');%legend('模型3');
    
    title('模型一估计');
    xlabel('cycleTimes'); ylabel('Capacity/%');
end;
%%
for t=predict_Start+1:cycleTimes
    x1_IMM(t)=0;
    x2_IMM(t)=0;
    x3_IMM(t)=0;
    for i=1:N
        s1(t,i)=f1(s1(t-1,i),t,model);
        s2(t,i)=f2(s2(t-1,i),t,model);
        s3(t,i)=f3(s3(t-1,i),t,model);
        x1_IMM(t)=x1_IMM(t)+s1(t,i)*q1(predict_Start,i);
        x2_IMM(t)=x2_IMM(t)+s2(t,i)*q2(predict_Start,i);
        x3_IMM(t)=x3_IMM(t)+s3(t,i)*q3(predict_Start,i);
    end;
    %x1_IMM(t)=x1_IMM(t)/N;
    %x2_IMM(t)=x1_IMM(t)/N;
    %x3_IMM(t)=x1_IMM(t)/N;
    x_result(t)=x1_IMM(t)*u_IMM(1,predict_Start)+x2_IMM(t)*u_IMM(2,predict_Start)+x3_IMM(t)*u_IMM(3,predict_Start);
    plot(t,x1_IMM(t),'y*');
    plot(t,x2_IMM(t),'g*');
    plot(t,x3_IMM(t),'b*');
    grid on ;hold on
end;
h1=plot(2:cycleTimes,x_result(2:cycleTimes),'k','linewidth',1);
h2=plot(2:cycleTimes,x1_IMM(2:cycleTimes),'y','linewidth',1);
h3=plot(2:cycleTimes,x2_IMM(2:cycleTimes),'g','linewidth',1);
h4=plot(2:cycleTimes,x3_IMM(2:cycleTimes),'b','linewidth',1);
h5=plot(initial_data,'r','linewidth',1);
plot([0 cycleTimes],[0.72 0.72],'c');
legend([h1 h2 h3 h4 h5],'综合结果','模型1','模型2','模型3','真实数据');
hold off
%%
figure
plot(initial_data,'r','linewidth',1);grid on;hold on
plot(2:cycleTimes,x_result(2:cycleTimes),'b','linewidth',1);
plot([0 cycleTimes],[0.72 0.72],'c');
hold off
title('滤波结果与真实值的对比');
xlabel('cycleTimes'); ylabel('Capacity/%');
legend('真实数据','综合结果');
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
fc=zeros(1,cycleTimes);
fc1=zeros(1,cycleTimes);
fc2=zeros(1,cycleTimes);
fc3=zeros(1,cycleTimes);
for t=1:cycleTimes
    d(t)=abs(x_result(t)-initial_data(t));
    d1(t)=abs(x1_IMM(t)-initial_data(t));
    d2(t)=abs(x2_IMM(t)-initial_data(t));
    d3(t)=abs(x3_IMM(t)-initial_data(t));
    fc(t)=d(t)*d(t);
    fc1(t)=d1(t)*d1(t);
    fc2(t)=d2(t)*d2(t);
    fc3(t)=d3(t)*d3(t);
end;
figure
plot(d,'k','linewidth',1);grid on;hold on
plot(d1,'y','linewidth',1);
plot(d2,'g','linewidth',1);
plot(d3,'b','linewidth',1);
legend('综合结果','模型1','模型2','模型3');
hold off
title('误差分析');
%%
%平均误差
figure
plot([0 cycleTimes],[sum(d)/cycleTimes sum(d)/cycleTimes],'k');hold on;
plot([0 cycleTimes],[sum(d1)/cycleTimes sum(d1)/cycleTimes],'y');
plot([0 cycleTimes],[sum(d2)/cycleTimes sum(d2)/cycleTimes],'g');
plot([0 cycleTimes],[sum(d3)/cycleTimes sum(d3)/cycleTimes],'b');
legend('综合结果','模型1','模型2','模型3');
hold off
title('平均误差');
%%
figure
plot([0 cycleTimes],[sum(fc)/cycleTimes sum(fc)/cycleTimes],'k');hold on;
plot([0 cycleTimes],[sum(fc1)/cycleTimes sum(fc1)/cycleTimes],'y');
plot([0 cycleTimes],[sum(fc2)/cycleTimes sum(fc2)/cycleTimes],'g');
plot([0 cycleTimes],[sum(fc3)/cycleTimes sum(fc3)/cycleTimes],'b');
legend('综合结果','模型1','模型2','模型3');
hold off
title('方差');
