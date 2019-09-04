% 'zui su xia jiang fa'
% this algorithm is 'steepest descent method' 
% with 'one-dimension-search':newton method
function []=steepest_descent_method_with_newton_lambda_k(func,x1,x2,epsilon,method,max)
% method parameter can be range 1-3 
% 1='using the best lambda_k'; 2=Armijo-Goldstern rule; 3=Wolfe-Powell rule
% x1&x2 are the initial points of iteration
% epsilon is the prescision of this algorithm
global result
result=[inf;inf];
global n
n=0;
global k
k=max;
func2(func,x1,x2,epsilon,method);

function []=func2(func,x1,x2,epsilon,method)
global k
global result
global n
variable_vector=symvar(func);
g_k=subs(gradient(func,[variable_vector(1),variable_vector(2)]),[variable_vector(1),variable_vector(2)],[x1,x2]);
if norm(g_k)<=epsilon || n==k
    result=[x1;x2];
else
    func3(func,x1,x2,epsilon,method,g_k);
end

function []=func3(func,x1,x2,epsilon,method,g_k)
% define the direction of searches
d_k=-g_k;
func4(func,x1,x2,epsilon,method,d_k)

function []=func4(func,x1,x2,epsilon,method,d_k)
global newton_result
newton_result=inf;
variable_vector=symvar(func);
% H=jacobian(gradient(func,[variable_vector(1),variable_vector(2)]));
if method==1
    % call one-dimension-search(newton-method) to find the best lambda_k
    % construct the function of lambda
    syms lambda
    x1_for_lambda_k=x1+lambda*d_k(1,1);
    x2_for_lambda_k=x2+lambda*d_k(2,1);
    % variable_vector=symvar(func);
    func_body=subs(func,[variable_vector(1),variable_vector(2)],[x1_for_lambda_k,x2_for_lambda_k]);
    % set the initial number of 'bu chang'=lambda0
    lambda0=1;
    % set the initial number of 'percision'=epsilon0
    epsilon0=10.^(-3);
    ols_newton(func_body,lambda0,epsilon0)
    lambda_k=newton_result;
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


% -------------- newton method ------------------------------------

function []=ols_newton(func,lambda,epsilon)
% the variable of func should be x but not others!!!
% lambda get the initial value of x lambda=x0
global newton_result
newton_result=inf;
newton_func2(func,lambda,epsilon);

function []=newton_func2(func,lambda,epsilon)
global newton_result
variable_vector=symvar(func);
% get the first-order derivative
diff_func_1=diff(func,variable_vector,1);
x=lambda;
% variable_vector=lambda;
phai_diff_func_1=eval(diff_func_1);
% get the second-order derivative
diff_func_2=diff(func,variable_vector,2);
x=lambda;
% variable_vector=lambda;
phai_diff_func_2=eval(diff_func_2);


if abs(phai_diff_func_1)<=epsilon
    newton_result=lambda;
else
    newton_func3(func,lambda,epsilon,phai_diff_func_1,phai_diff_func_2);
end

function []=newton_func3(func,lambda,epsilon,phai_diff_func_1,phai_diff_func_2)
% create the new iteration points
lambda=lambda-(phai_diff_func_1/phai_diff_func_2);
newton_func2(func,lambda,epsilon)
   
