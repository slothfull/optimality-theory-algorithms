% this algorithm can be used for searching 'best-lambda_k'='zui you bu chang'
clear all
clc
global result
% lambda get the initial value of x lambda=x0
lambda=6;
epsilon=10.^(-3);
syms x
% the variable of func should be x but not others!!!
func=x.^4-4*x.^3-6*x.^2-16*x+4;
ols_newton(func,lambda,epsilon);
% result is the x(min-value) find by algorithm newton
result;
% result=4.0000
x=4;
min_value_by_newton=eval(func);
%% plot figures
x=-1:0.01:5;
y=x.^4-4*x.^3-6*x.^2-16*x+4;
figure(1)
hold on
plot(x,y,'r');
% use matlab function to find the min
y1=@(x)x.^4-4*x.^3-6*x.^2-16*x+4;
[x,min_value]=fminbnd(y1,-1,5);
plot(x,min_value,'o');
xlabel('location')
ylabel('temperature')
title('The curve of temperature changing with the depth of the probe')
hold off