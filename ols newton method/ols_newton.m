function []=ols_newton(func,lambda,epsilon)
%% the variable of func should be x but not others!!!
% lambda get the initial value of x lambda=x0
global result
result=inf;
func2(func,lambda,epsilon);

function []=func2(func,lambda,epsilon)
%%
global result
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
    result=lambda;
else
    func3(func,lambda,epsilon,phai_diff_func_1,phai_diff_func_2);
end

function []=func3(func,lambda,epsilon,phai_diff_func_1,phai_diff_func_2)
%% create the new iteration points
lambda=lambda-(phai_diff_func_1/phai_diff_func_2);
func2(func,lambda,epsilon)