% this newton algorithm is specially for the functions
% which are two-variables functions
function []=newton(func,x1,x2,epsilon,lambda,n)
% func is the objective function
% x is the initial value of the func
% epsilon is the presicion of this algorithm
% lambda is the length of step of this func
global result
result=[inf;inf];
global k
k=0;
% x_k = [x1;x2]
func2(func,x1,x2,epsilon,lambda,n);

function []=func2(func,x1,x2,epsilon,lambda,n)
global result
global k
syms x y
variable_vector=symvar(func);
% using symvar(func) to get the variables x&y in the function 
if norm(subs(gradient(func,[variable_vector(1),variable_vector(2)]),[variable_vector(1),variable_vector(2)],[x1,x2]))<=epsilon || k==n
    result=[x1;x2];
else
    func3(func,x1,x2,epsilon,lambda,n);
end

function []=func3(func,x1,x2,epsilon,lambda,n)
syms x y
variable_vector=symvar(func);
% using symvar(func) to get the variables x&y in the function 
% create the direction of algorithm d_k
gradient_func=subs(gradient(func,[variable_vector(1),variable_vector(2)]),[variable_vector(1),variable_vector(2)],[x1,x2]);
d_k=-subs(inv(jacobian(gradient(func,[variable_vector(1),variable_vector(2)]))),[variable_vector(1),variable_vector(2)],[x1,x2])*gradient_func;
func4(func,x1,x2,epsilon,lambda,d_k,n)

function []=func4(func,x1,x2,epsilon,lambda,d_k,n)
global k
% confirm the new iteration x_k point
% x_k+1 = x_k + d_k
x1=x1+lambda*d_k(1,1);
x2=x2+lambda*d_k(2,1);
k=k+1;
func2(func,x1,x2,epsilon,lambda,n)


