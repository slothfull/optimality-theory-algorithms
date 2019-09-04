clear all
clc
syms x y
func=x.^3+x*y+(1+y).^2;
H_k=[1,0;0,1];
epsilon=10.^(-6);
% 0: using best_lambda_k 
% 1: using wolfe_powell
wolfe_powell=0;
x1=0;
x2=0;
lambda0=inf;
n=2;
DFP_Quasi_newton_with_newton_for_lambda_k(func,x1,x2,H_k,epsilon,lambda0,wolfe_powell,n);
%% 
global result
result;
% x=temporal_result(1);
% y=temporal_result(2);
% vpa(eval(func),10)