%%
%model=1--A1 
%model=2--A2
%model=3--A3
%model=4--A4
%model=5--A4�������Ƶõ�����
%%
function [x]=F3(x_Pre,k,model)
%vģ��
if (model==1)
    a=0.0291;
    b=0.0316;
end;
if (model==2)
    a=0.0284;
    b=0.0308;
end;
if (model==3)
    a=0.0268;
    b=0.0291;
end;
if (model==4)
      a=0.0103087;
      b=0.0108064; 
end;
if (model==5)
    %27
    a=0.00533011222533502;
    b=0.00759227883810740;
end;
x=1./((b/a+(1./x_Pre-b/a))*exp(a));