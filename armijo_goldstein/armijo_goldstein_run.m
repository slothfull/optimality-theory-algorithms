clear all;
global result
lambda=-1.2;
sigma1=0.2;
sigma2=0.8;
alpha=1.5;

syms x
% this function should be the phai(lambda) but not the function f(x)
% while phai(lambda)=f(x_k+lambda*d_k)  x_k & d_k in one iteration in 
% one-dimension search canbe a constant number
func=3*x.^4-4*x.^3-12*x.^2;

armijo_goldstein(func,lambda,sigma1,sigma2,alpha);
% result is the acceptable 'bu chang' lambda_k
x=result;
