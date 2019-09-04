% 'zui su xia jiang fa'
% this algorithm only suitable for quadratic programming problem
% method parameter can be range 1-3 
% 1='using the best lambda_k'; 2=Armijo-Goldstern rule; 3=Wolfe-Powell rule
% x1&x2 are the initial points of iteration
% epsilon is the prescision of this algorithm
clear all
clc
global result
global n
syms x y
func=x.^3+x*y+(1+y).^2;
x1=0;
x2=0;
epsilon=10.^(-6);
method=1;
% max holds the number of max iteration
max=2;
steepest_descent_method_with_newton_lambda_k(func,x1,x2,epsilon,method,max)
result;
n;

% get the result
x0=result(1,1);
y0=result(2,1);
% x0=vpa(x0,4)
% y0=vpa(y0,4)
x=x0;
y=y0;
value_func=eval(func);

% find the min between[a,b]
fun=@(x)x(1).^2+x(1)*x(2)+(1+x(2)).^2;
min=fminsearch(fun,[100,100]);

% draw the pic of the func to find the min
[x,y] = meshgrid(-10:0.02:10);
func=x.^3+x*y+(1+y).^2;
% func(x>10) = nan;
% func(y>10) = nan;
% func(x<-10) = nan;
% func(y<-10) = nan;
surf(x,y,func);
axis tight
