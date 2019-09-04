% p280
function []=wolfe_powell(func,lambda,sigma1,sigma3,alpha,label)
global a
a=0;
global b
b=inf;
global result
result=inf;
func2(func,lambda,sigma1,sigma3,alpha,label);


function []=func2(func,lambda,sigma1,sigma3,alpha,label)
% syms x
% f=3*x.^4-4*x.^3-12*x.^2;
% y=diff(f,x,1);
% x=0;
% result_diff=eval(y);
variable_vector=symvar(func);
diff_func=diff(func,variable_vector(1),1);
x=0;
result_diff=eval(diff_func);
global b
inline_func=inline(func);
if inline_func(lambda)<=inline_func(0)+sigma1*result_diff*lambda
    func3(func,lambda,sigma1,sigma3,alpha,label);
else
    b=lambda;
    func4(func,lambda,sigma1,sigma3,alpha,label);
end


function []=func3(func,lambda,sigma1,sigma3,alpha,label)
% syms x
% f=3*x.^4-4*x.^3-12*x.^2;
% y=diff(f,x,1);
% x=0;
% result_diff=eval(y);
variable_vector=symvar(func);
diff_func=diff(func,variable_vector(1),1);
x=0;
result_diff=eval(diff_func);
x=lambda;
result_diff_1=eval(diff_func);
global a
global result
inline_func=inline(func);
if inline_func(lambda)>=sigma3*result_diff && label==1
    result=lambda;
elseif abs(inline_func(lambda))<=-sigma3*result_diff && label==2
    result=lambda;
elseif 0>=result_diff_1 && result_diff_1>=sigma3*result_diff && label==3
    result=lambda;
else
    a=lambda;
    func4(func,lambda,sigma1,sigma3,alpha,label);
end
    
function []=func4(func,lambda,sigma1,sigma3,alpha,label)
global a
global b
if b<inf
    lambda=(a+b)./2;
else
    lambda=alpha*lambda;
end
func2(func,lambda,sigma1,sigma3,alpha,label)