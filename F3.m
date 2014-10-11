function [x]=F3(x_Pre,k,model)
%vÄ£ÐÍ
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
% if (model==8)
%     %27
%     a=0.00533011222533502;
%     b=0.00759227883810740;
% end;
x=1./((b/a+(1./x_Pre-b/a))*exp(a));