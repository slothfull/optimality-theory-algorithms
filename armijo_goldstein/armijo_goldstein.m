% p280
function []=armijo_goldstein(func,lambda,sigma1,sigma2,alpha)
global a
a=0;
global b
b=inf;
global result
result=inf;
func2(func,lambda,sigma1,sigma2,alpha);


function []=func2(func,lambda,sigma1,sigma2,alpha)
% syms x
% f=3*x.^4-4*x.^3-12*x.^2;
% y=diff(f,x,1);
% x=0;
% result_diff=eval(y);

% the code above can be replaced by the single line below
% take care of the func: polyder & polyval

% result_diff=polyval(polyder([3 -4 -12 0 0]),-1.2);

% another choice of getting the result using one line is:
% take care of the func: subs() & sym
% result=subs(diff(sym('3*x^4-4*x^3-12*x^2')),-1.2)

variable_vector=symvar(func);
diff_func=diff(func,variable_vector(1),1);
x=0;
result_diff=eval(diff_func);

global b
inline_func=inline(func);
if inline_func(lambda)<=inline_func(0)+sigma1*result_diff*lambda
    func3(func,lambda,sigma1,sigma2,alpha);
else
    b=lambda;
    func4(func,lambda,sigma1,sigma2,alpha);
end


function []=func3(func,lambda,sigma1,sigma2,alpha)
variable_vector=symvar(func);
diff_func=diff(func,variable_vector(1),1);
x=0;
result_diff=eval(diff_func);

global a
global result
inline_func=inline(func);
if inline_func(lambda)>=inline_func(0)+sigma2*result_diff*lambda
    result=lambda;
else
    a=lambda;
    func4(func,lambda,sigma1,sigma2,alpha);
end
    
function []=func4(func,lambda,sigma1,sigma2,alpha)
global a
global b
if b<inf
    lambda=(a+b)./2;
else
    lambda=alpha*lambda;
end
func2(func,lambda,sigma1,sigma2,alpha)

