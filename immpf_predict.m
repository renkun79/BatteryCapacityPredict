clear all;
clc;
echo off;
%%
%Interactive Multiple Model Particle Filter
%Model:Battery
%Author:Ye Zhe
%CreateDate:2014.7.9
%EditDate
%2014/10/11 by Ye Zhe
%�����Ԥ��ʧЧʱ�̵�ģ��
%School:Harbin Institute of Technology,China
%%
%�޸ļ�¼:
%ȥ��F1 F2 F3Ԥ��Ĳ�����
%%
%˵��:��4������a1shiyan.mat,a2shiyan.mat,a3shiyan.mat,a4shiyan.mat���������ݷֱ��Ӧ
%model��ֵΪ1,2,3,4
%%
%%%��ʼ��%%%%%%%
N=500;       %������Ŀ
f1=@F1;
f2=@F2;
f3=@F3;
% ������������
%%
%%%����ģ����Ҫ�޸ĵĲ���%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('a4shiyan.mat');
initial_data=a4;%ԭʼ��������
model=4;
predict_Start=32;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
cycleTimes=length(initial_data);
x0=initial_data(1);  % ��ʼ״̬
%%
%������������
figure
plot(initial_data,'r','linewidth',1);grid on;hold on
plot([0 cycleTimes],[0.72 0.72],'c');

hold off
title('����ֵ');
xlabel('cycleTimes'); ylabel('capacity/%');
legend('��ʵ����');
%%
%%%%%%%%%%%%%%%%%%%%%%%%IMMPF%%%%%%%%%%%%%%%%%%%%%
%��ʼ��
s1=zeros(cycleTimes,N);%����ģ�͵����ӺͶ�ӦȨֵ
s2=zeros(cycleTimes,N);
s3=zeros(cycleTimes,N);
s=zeros(cycleTimes,N);%�ۺ�֮�������
q1=zeros(cycleTimes,N);
q2=zeros(cycleTimes,N);
q3=zeros(cycleTimes,N);
for i=1:N
    s1(1,i)=x0;
    s2(1,i)=x0;
    s3(1,i)=x0;
    s(1,i)=x0;
end;

x1_IMM = zeros(cycleTimes,1);      %ģ��1IMM�㷨״̬����ֵ
x2_IMM = zeros(cycleTimes,1);      %ģ��2IMM�㷨״̬����ֵ
x3_IMM = zeros(cycleTimes,1);      %ģ��3IMM�㷨״̬����ֵ
x_pro_IMM = zeros(1,cycleTimes);   %IMM�㷨ģ���ۺ�״̬����ֵ
x_result=zeros(1,cycleTimes);
x_result(1)=x0;
P_IMM=zeros(1,cycleTimes);       %IMM�㷨ģ���ۺ�״̬Э����
P1_IMM=0;
P2_IMM=0;
P3_IMM=0;              %IMM�㷨��ģ��Э����
r1_IMM=0;
r2_IMM=0;
r3_IMM=0;
S1_IMM=0;
S2_IMM=0;
S3_IMM=0;
x_pro_IMM(:,1)=x0;

%ģ��ת�Ƹ��ʾ���
pij=[0.7,0.25,0.05;
    0.25,0.7,0.05;
    0.05,0.05,0.9];

u_IMM=zeros(3,cycleTimes);
u_IMM(:,1)=[0.2,0.3,0.5]';  %IMM�㷨ģ�͸���

x1_IMM=x0;x2_IMM=x0;x3_IMM=x0;  %IMM�㷨��ģ�ͳ�ʼ״̬

%����ģ�͵�Likelihoods
L1=zeros(1,cycleTimes);
L2=zeros(1,cycleTimes);
L3=zeros(1,cycleTimes);
figure
%%
for t=2:predict_Start
    x1_IMM(t) = 0;      %ģ��1IMM�㷨״̬����ֵ
    x2_IMM(t) = 0;      %ģ��2IMM�㷨״̬����ֵ
    x3_IMM(t) = 0;      %ģ��3IMM�㷨״̬����ֵ
     %��һ��Interacting��ֻ���IMM�㷨��
    c_j=pij'*u_IMM(:,t-1);
    
    %����ģ�ͻ�ϸ���
    ui1=(1/c_j(1))*pij(:,1).*u_IMM(:,t-1);
    ui2=(1/c_j(2))*pij(:,2).*u_IMM(:,t-1);
    ui3=(1/c_j(3))*pij(:,3).*u_IMM(:,t-1);    

    %�����˲�����
    [s1(t,:),q1(t,:),L1(t)]=ParticleFilter(s1(t-1,:),initial_data(t),f1,t,model);
    [s2(t,:),q2(t,:),L2(t)]=ParticleFilter(s2(t-1,:),initial_data(t),f2,t,model);
    [s3(t,:),q3(t,:),L3(t)]=ParticleFilter(s3(t-1,:),initial_data(t),f3,t,model);
    
    %ģ�͸��ʸ���
    c=L1(t)*c_j(1)+L2(t)*c_j(2)+L3(t)*c_j(3);
    u_IMM(1,t)=(1/c)*L1(t)*c_j(1);
    u_IMM(2,t)=(1/c)*L2(t)*c_j(2);
    u_IMM(3,t)=(1/c)*L3(t)*c_j(3);
    %ģ���ۺ�
    for i=1:N
        x1_IMM(t)=x1_IMM(t)+s1(t,i)*q1(t,i);
        x2_IMM(t)=x2_IMM(t)+s2(t,i)*q2(t,i);
        x3_IMM(t)=x3_IMM(t)+s3(t,i)*q3(t,i);
    end;
    x_result(t)=x1_IMM(t)*u_IMM(1,t)+x2_IMM(t)*u_IMM(2,t)+x3_IMM(t)*u_IMM(3,t);
    for i=1:N
        %s(t,i)=x_result(t);
        %%%%%%�˴�����ģ�����ӵ��ںϷ�ʽ���Ը���%%%%%%%%%%%%
        s(t,i)=s1(t,i)*q1(t,i)*u_IMM(1,t)+s2(t,i)*q2(t,i)*u_IMM(2,t)+s3(t,i)*q3(t,i)*u_IMM(3,t);
    end;
    %plot(t,x_result(t),'linewidth',4);
    plot(t,x1_IMM(t),'y*')%legend('ģ��1');
    grid on ;hold on;
    plot(t,x2_IMM(t),'g*');%legend('ģ��2');
    plot(t,x3_IMM(t),'b*');%legend('ģ��3');
    
    title('ģ��һ����');
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
legend([h1 h2 h3 h4 h5],'�ۺϽ��','ģ��1','ģ��2','ģ��3','��ʵ����');
hold off
%%
figure
plot(initial_data,'r','linewidth',1);grid on;hold on
plot(2:cycleTimes,x_result(2:cycleTimes),'b','linewidth',1);
plot([0 cycleTimes],[0.72 0.72],'c');
hold off
title('�˲��������ʵֵ�ĶԱ�');
xlabel('cycleTimes'); ylabel('Capacity/%');
legend('��ʵ����','�ۺϽ��');
%%
%������
%               d ------����ʽ�����
%               d1------ģ��1�����
%               d2------ģ��2�����
%               d3------ģ��3�����
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
legend('�ۺϽ��','ģ��1','ģ��2','ģ��3');
hold off
title('������');
%%
%ƽ�����
figure
plot([0 cycleTimes],[sum(d)/cycleTimes sum(d)/cycleTimes],'k');hold on;
plot([0 cycleTimes],[sum(d1)/cycleTimes sum(d1)/cycleTimes],'y');
plot([0 cycleTimes],[sum(d2)/cycleTimes sum(d2)/cycleTimes],'g');
plot([0 cycleTimes],[sum(d3)/cycleTimes sum(d3)/cycleTimes],'b');
legend('�ۺϽ��','ģ��1','ģ��2','ģ��3');
hold off
title('ƽ�����');
%%
figure
plot([0 cycleTimes],[sum(fc)/cycleTimes sum(fc)/cycleTimes],'k');hold on;
plot([0 cycleTimes],[sum(fc1)/cycleTimes sum(fc1)/cycleTimes],'y');
plot([0 cycleTimes],[sum(fc2)/cycleTimes sum(fc2)/cycleTimes],'g');
plot([0 cycleTimes],[sum(fc3)/cycleTimes sum(fc3)/cycleTimes],'b');
legend('�ۺϽ��','ģ��1','ģ��2','ģ��3');
hold off
title('����');
%%
%����RUL��ģ��
%�����ڵ�ʱ�䲽������ÿ�����ӣ���ÿ�����ӵ���������һ���۲⣬��¼������ʧЧʱ�䣬�����ͳ��
%ģ��1��RUL����
RUL1=zeros(1,N);
for i=1:N
    temp=find(s1(:,i)<0.72);
    RUL1(i)=temp(1);
end;
figure;
subplot(2,1,1);
hist(RUL1,max(RUL1)-min(RUL1));
xlim([min(RUL1) max(RUL1)]);
xlabel('RUL');
title('ģ��1RUL');

%ģ��2��RUL����
RUL2=zeros(1,N);
for i=1:N
    temp=find(s2(:,i)<0.72);
    if (length(temp)>0)
      RUL2(i)=temp(1);
    end;
end;
%ģ��3��RUL����
RUL3=zeros(1,N);
for i=1:N
    temp=find(s3(:,i)<0.72);
    RUL3(i)=temp(1);
end;
subplot(2,1,2);
hist(RUL3,max(RUL3)-min(RUL3));
xlim([min(RUL3) max(RUL3)]);
xlabel('RUL');
title('ģ��3RUL');















