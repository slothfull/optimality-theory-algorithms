% this method is to solve "fei xian xing objective func,xian xing yue shu"
% this code is special for f(x1,x2) problem
% func=objective function;[x1;x2]=initial point;epsilon=algorithm accuracy
% n=max iteration steps; A&b=limit funcs: Ax<=b;Aeq&Beq=limit funcs:Aeqx=Beq
% ib is the lower bound of x1x2
% ub is the upper bound of x1x2
% if there no need to transmit certain parameter, get it=inf;
clear;
clc;
syms x y
func=x.^2+4*y.^2-34*x-32*y;
x1=1;x2=2;
epsilon=10.^(-3);
n=2;
A=[2,1;0,1];
b=[6;2];
Aeq=inf;
Beq=inf;
ib=[0;0];
ub=inf;
frank_wolfe(func,x1,x2,epsilon,n,A,b,Aeq,Beq,ib,ub)
global result;
result;
global temporal_result
temporal_result;
global optimality_discrimination_criteria
if optimality_discrimination_criteria==0
    disp('the result derived from frank_wolfe method is K-T point')
end
x1=result(1,1);
y1=result(2,1);
f=x1.^2+4*y1.^2-34*x1-32*y1;
% usage of linprog() can be found in 
% https://blog.csdn.net/ten_sory/article/details/54571525
