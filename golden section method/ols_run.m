% function f=object_func(x)
% f=exp(-x)+x^2;
% func=@object_func;
a=0;
b=1;
gold_ols(inline('exp(-x)+x.^2'),a,b);
global result
result;
x=0:1/1000:1;
y=exp(-x)+x.^2;
plot(x,y);
min=fminbnd('exp(-x)+x.^2',a,b);
