% 'zui su xia jiang fa'
% this algorithm only suitable for quadratic convex problem
function []=steepest_descent_method(func,x1,x2,epsilon,method)
% method parameter can be range 1-3 
% 1='using the best lambda_k'; 2=Armijo-Goldstern rule; 3=Wolfe-Powell rule
% x1&x2 are the initial points of iteration
% epsilon is the prescision of this algorithm
global result
result=[inf;inf];
global n
n=1;
func2(func,x1,x2,epsilon,method);

function []=func2(func,x1,x2,epsilon,method)
global result
variable_vector=symvar(func);
g_k=subs(gradient(func,[variable_vector(1),variable_vector(2)]),[variable_vector(1),variable_vector(2)],[x1,x2]);
if norm(g_k)<=epsilon
    result=[x1;x2];
else
    func3(func,x1,x2,epsilon,method,g_k);
end

function []=func3(func,x1,x2,epsilon,method,g_k)
% define the direction of searches
d_k=-g_k;
func4(func,x1,x2,epsilon,method,g_k,d_k)

function []=func4(func,x1,x2,epsilon,method,g_k,d_k)
variable_vector=symvar(func);
H=jacobian(gradient(func,[variable_vector(1),variable_vector(2)]));
if method==1
    % 'zui you bu chang' use the function to set it
    % the 'formular' works only if the func is 'er ci tu gui hua'
    lambda_k=-((g_k')*d_k)/((d_k')*H*d_k);
elseif method==2
    % call armijo_goldstein function to find the acceptable lambda
elseif method==3
    % call wolfe_powell function to find the acceptable lambda
end
func5(func,x1,x2,epsilon,method,lambda_k,d_k);

function []=func5(func,x1,x2,epsilon,method,lambda_k,d_k)
global n
x1=x1+lambda_k*d_k(1,1);
x2=x2+lambda_k*d_k(2,1);
n=n+1;
func2(func,x1,x2,epsilon,method)

   
