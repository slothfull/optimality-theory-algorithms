clear
clc
global result
syms x y
% func=(6+x+y).^2+(2-3*x-3*y-x*y).^2;
func=x.^3+x*y+(1+y).^2;
% x1 x2 iteration initial value
x1=0;
x2=0;
epsilon=10^(-6); % iteration accuracy
lambda=1; % iteration step
n=2; % max number of iteration steps
newton(func,x1,x2,epsilon,lambda,n);
format rat
result;
vpa(result,10)
x=result(1,1);
y=result(2,1);
vpa(eval(func),10)
