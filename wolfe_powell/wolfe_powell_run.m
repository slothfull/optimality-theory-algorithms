clear all;
format short
global result
% inital 'bu chang' lambda_k
lambda=-1.2;
% acceptable parameters
sigma1=0.2;
sigma3=0.8;
% 'jia bu can shu' alpha
alpha=1.5;
% label intend to choose 1 from 3 step-limit conditions
label=1;

syms x
% this function should be the phai(lambda) but not the function f(x)
% while phai(lambda)=f(x_k+lambda*d_k)  x_k & d_k in one iteration in 
% one-dimension search canbe a constant number
func=3*x.^4-4*x.^3-12*x.^2;

wolfe_powell(func,lambda,sigma1,sigma3,alpha,label);
x=result;