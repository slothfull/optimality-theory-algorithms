% UQP conjugate gradient method
% func takes two variables
% method parameter takes 1 2 3 which stand for Fletcher-Reeves & Dixon &
% Polak-Ribiere-Polyak method
% n is defined as the epochs of iterations
% x1 & x2 is the intial points
% H is derived from func
clear
clc
global result
syms x y
func=x.^3+x*y+(1+y).^2;
% func=x.^2+2*y.^2;
% attention the matlab function collect()
rectified_func=collect(func,x);

% x1 x2 iteration initial value
x1=0;
x2=0;
n=2; % max number of iterations
method=2; % method = 1 or 2 or 3
H=jacobian(gradient(func));
conjugate_gradient_method_with_newton_lambda_k(rectified_func,method,n,x1,x2,H)

result;
vpa(result,10) % Keep valid figures
x=result(1,1);
y=result(2,1);
vpa(eval(func),10)