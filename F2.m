%%
%model=1--A1 
%model=2--A2
%model=3--A3
%model=4--A4
%model=5--A4参数估计得到参数
%%
function [x]=F2(x_Pre,k,model)
%%
%多项式模型
%参数说明
%x_Pre                 --前时刻状态
%  k                   --该时刻为第k个状态
%%
if model==1
    a=-3.3611e-005;
    b=-2.3363e-004;
end;
if model==2
    a=-3.4696e-005;
    b=-102500e-005;
      
end;
if model==3
    a=-3.5144e-005;
    b=2.0913e-004;
end;
if model==4
    a=-4.6603e-006;
    b=-1.2930e-004;   
end;
if model ==5
    a=-9.49109136977370e-05;
   b=-0.000210675223360429;
end;
% %%
% %参数估计测试模式
% if model==8
%    %27
%    a=-9.49109136977370e-05;
%    b=-0.000210675223360429;
% end;
x=x_Pre+2*a*k+a+b;
