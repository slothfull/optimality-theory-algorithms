syms x
f=@(x)3*x(1).^3 + 2*x(2).^2 - 12*x(1);
% x0=[0,1];
% H=my_hessian(f,x,y,x0);
H1=my_hessian(f,x);