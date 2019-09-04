% this Quasi-Newton Methods algorithm is specially for the functions
% which are two-variables functions
% #####and attention that this function use 'ke jie shou bu chang'####
% 'ke jie shou bu chang' is derived from wolfe_powell method
% 'zui you bu chang' is derived form 'xian cheng de gong shi'
% with the limit :'wu yue shu tu er ci gui hua'
function []=DFP_Quasi_newton(func,x1,x2,H1,epsilon,lambda0,wolfe_powell,n)
% func is the objective function
% x1 & x2 is the initial value of the func
% H_k is the initial matrix of H
% epsilon is the presicion of this algorithm (epsilon>0)
% lambda0 is the length of step of this func
% wolfe_powell is a parameter to control whether to use wolfe_powell method 
% or the best-lambda_k method
% wolfe_powell=1 using wolfe_powell
% wolfe_powell=0 using 'best-lambda_k method'
% n is the number of max iteration times
% global variable H holds the initial matrix for func5() to use
global H
H=H1;
% global variable k defines the number of iteration steps
global k
k=1;
global result
variable_vector=symvar(func);
if norm(subs(gradient(func,[variable_vector(1),variable_vector(2)]),[variable_vector(1),variable_vector(2)],[x1,x2]))<=epsilon
    result=[x1;x2];
else
    g_k_new=subs(gradient(func,[variable_vector(1),variable_vector(2)]),[variable_vector(1),variable_vector(2)],[x1,x2]);
end
func2(func,x1,x2,H1,epsilon,lambda0,g_k_new,wolfe_powell,n)

function []=func2(func,x1_new,x2_new,H_k,epsilon,lambda_k,g_k_new,wolfe_powell,n)
global result
d_k=-H_k*g_k_new;
if d_k==0
    result=[x1_new;x2_new];
else
    func3(func,x1_new,x2_new,H_k,epsilon,lambda_k,g_k_new,d_k,wolfe_powell,n)
end


function []=wolfe_powell_method(func,lambda,sigma1,sigma3,alpha,label)
global a
a=0;
global b
b=inf;
global lambda_result
lambda_result=inf;
wp_func2(func,lambda,sigma1,sigma3,alpha,label);

function []=func3(func,x1_old,x2_old,H_k,epsilon,lambda_k,g_k_old,d_k,wolfe_powell,n)
global lambda_result
lambda_result=inf;
syms x y
variable_vector=symvar(func);
% using symvar(func) to get the variables x&y in the function 
% create the direction of algorithm d_k
gradient_func=subs(gradient(func,[variable_vector(1),variable_vector(2)]),[variable_vector(1),variable_vector(2)],[x1_old,x2_old]);
% determine the lambda_k
if wolfe_powell==0
    Hessian_func=jacobian(gradient(func,[variable_vector(1),variable_vector(2)]));
    lambda_k=-((gradient_func')*d_k)/((d_k')*Hessian_func*d_k);
elseif wolfe_powell==1
    % initial 'bu chang' lambda_k
    initial_lambda=-1.2;
    % acceptable parameters
    sigma1=0.2;
    sigma3=0.8;
    % 'jia bu can shu' alpha
    alpha=1.5;
    % label intend to choose 1 from 3 step-limit conditions
    label=1;
    syms lambda
    % get the (x_k+lambda_k*d_k)
    x1_for_lambda_k=x1_old+lambda*d_k(1,1);
    x2_for_lambda_k=x2_old+lambda*d_k(2,1);
    variable_vector=symvar(func);
    % attention:code below can change the sym function into function handle
    % pay attention to the usage of <<matlabFunction()>>
    func_body=subs(func,[variable_vector(1),variable_vector(2)],[x1_for_lambda_k,x2_for_lambda_k]);
    % func_body=collect(func_body,lambda);
    % phai_func=matlabFunction(func_body);
    
    wolfe_powell_method(func_body,initial_lambda,sigma1,sigma3,alpha,label);
    lambda_k=lambda_result;
end
if lambda_k==inf
    disp('wrong lambda_k derived using wolfe_powell method')
else
    x1_new=x1_old+lambda_k*d_k(1,1);
    x2_new=x2_old+lambda_k*d_k(2,1);
    gradient_func_new=subs(gradient(func,[variable_vector(1),variable_vector(2)]),[variable_vector(1),variable_vector(2)],[x1_new,x2_new]);
    g_k_new=gradient_func_new;
    % gradient_func_old=subs(gradient(func,[variable_vector(1),variable_vector(2)]),[variable_vector(1),variable_vector(2)],[x1_old,x2_old]);
    % g_k_old=gradient_func_old;
end
func4(func,x1_new,x2_new,x1_old,x2_old,H_k,epsilon,lambda_k,g_k_new,g_k_old,wolfe_powell,n);

function []=func4(func,x1_new,x2_new,x1_old,x2_old,H_k,epsilon,lambda_k,g_k_new,g_k_old,wolfe_powell,n)
global result
if norm(g_k_new)<=epsilon
    result=[x1_new;x2_new];
else
    func5(func,x1_new,x2_new,x1_old,x2_old,H_k,epsilon,lambda_k,g_k_new,g_k_old,wolfe_powell,n)
end

function []=func5(func,x1_new,x2_new,x1_old,x2_old,H_k,epsilon,lambda_k,g_k_new,g_k_old,wolfe_powell,n)
global H
global k
if k==n
    H_k=H;
else
    p_k=[x1_new-x1_old;x2_new-x2_old];
    q_k=g_k_new-g_k_old;
    H_k=H_k+(p_k*(p_k'))/((p_k')*q_k)-(H_k*q_k*(q_k')*H_k)/((q_k')*H_k*q_k);
end
k=k+1;
func2(func,x1_new,x2_new,H_k,epsilon,lambda_k,g_k_new,wolfe_powell,n);

%-------------------------------------------------------------------


function []=wp_func2(func,lambdak,sigma1,sigma3,alpha,label)
% syms x
% f=3*x.^4-4*x.^3-12*x.^2;
% y=diff(f,x,1);
% x=0;
% result_diff=eval(y);
variable_vector=symvar(func);
diff_func=diff(func,variable_vector(1),1);
lambda=0;
result_diff=eval(diff_func);
global b
inline_func=inline(func);
if inline_func(lambdak)<=inline_func(0)+sigma1*result_diff*lambdak
    wp_func3(func,lambdak,sigma1,sigma3,alpha,label);
else
    b=lambdak;
    wp_func4(func,lambdak,sigma1,sigma3,alpha,label);
end


function []=wp_func3(func,lambdak,sigma1,sigma3,alpha,label)
% syms x
% f=3*x.^4-4*x.^3-12*x.^2;
% y=diff(f,x,1);
% x=0;
% result_diff=eval(y);
variable_vector=symvar(func);
diff_func=diff(func,variable_vector(1),1);
lambda=0;
result_diff=eval(diff_func);
lambda=lambdak;
result_diff_1=eval(diff_func);
global a
global lambda_result
inline_func=inline(func);
if inline_func(lambdak)>=sigma3*result_diff && label==1
    lambda_result=lambdak;
elseif abs(inline_func(lambdak))<=-sigma3*result_diff && label==2
    lambda_result=lambdak;
elseif 0>=result_diff_1 && result_diff_1>=sigma3*result_diff && label==3
    lambda_result=lambdak;
else
    a=lambdak;
    wp_func4(func,lambdak,sigma1,sigma3,alpha,label);
end
    
function []=wp_func4(func,lambdak,sigma1,sigma3,alpha,label)
global a
global b
if b<inf
    lambda=(a+b)./2;
else
    lambda=alpha*lambdak;
end
wp_func2(func,lambda,sigma1,sigma3,alpha,label)

