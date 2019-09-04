% 'zui su xia jiang fa'
% this algorithm only suitable for quadratic programming problem
% method parameter can be range 1-3 
% 1='using the best lambda_k'; 2=Armijo-Goldstern rule; 3=Wolfe-Powell rule
% x1&x2 are the initial points of iteration
% epsilon is the prescision of this algorithm
clear all
clc
format short
global result
global n
syms x y
func=x.^2+x*y+(1+y).^2;
x1=0;
x2=0;
epsilon=10.^(-6);
method=1;
steepest_descent_method(func,x1,x2,epsilon,method)
result;
n;
x0=result(1,1);
y0=result(2,1);

% the result
x0=vpa(x0,4)
y0=vpa(y0,4)
vpa(eval(func),10)

% find the min between[a,b]
fun=@(x)x(1).^2+x(1)*x(2)+(1+x(2)).^2;
min=fminsearch(fun,[100,100]);
% draw the pic of the func to find the min
% [x,y] = meshgrid(0:0.02:1);
% func=x.^2+x*y+(1+y).^2;
% func(x>10) = nan;
% func(y>10) = nan;
% func(x<-10) = nan;
% func(y<-10) = nan;
% surf(x,y,func);
% axis tight


% n = 50;
% x = linspace(0,1,n+1);
% y = cumsum([zeros(size(x));repmat((1-x)/n,n,1)]);
% x = repmat(x,size(y,1),1);
% z = 3*(x.^2)+3*(y.^2)+3*x.*y+1-3*x-3*y;
% mesh(x,y,z);
% axis tight