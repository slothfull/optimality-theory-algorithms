clear all
clc
syms x y
func=x.^2+3*y.^2;
H_k=[2,1;1,1];
epsilon=10.^(-6);
wolfe_powell=1;
x1=1;
x2=1;
lambda0=inf;
n=2;
DFP_Quasi_newton(func,x1,x2,H_k,epsilon,lambda0,wolfe_powell,n);
global result
result;