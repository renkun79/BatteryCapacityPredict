function [x]=F1(x_Pre,k,model)
%%
%ָ��ģ��
if (model==1)
   b=0.08162;
   c=0.92022;
   d=-0.001767;
end;
if (model==2)
   b=0.08507;
   c=0.9115;
   d=-0.001833;
end;
if (model==3)
    b=0.05205;
    c=0.9155;
    d=-0.001667;
end;
if(model==4)
    b=0.0736;
    c=0.90267;
    d=-0.0009538;
end;
if(model==5)
    b=0.0734692178454160;
    c=0.895787596682329;
    d=-0.00203272773460600;
end;
if(model==8)
    %32
    b=0.0734692178454160;
    c=0.895787596682329;
    d=-0.00203272773460600;
    %27   
    b=0.0746896556149248;
    c=0.899756534898285;
    d=-0.00380786118912128;

end;
x=exp(b)*x_Pre+c*exp(d*k)*(exp(d)-exp(b));