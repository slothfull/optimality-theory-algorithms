%% frank-wolfe method
% this method is to solve "fei xian xing objective func,xian xing yue shu"
% this code is special for f(x1,x2) problem
function []=frank_wolfe(func,x1,x2,epsilon,n,A,b,Aeq,Beq,ib,ub)
%% main function
% func=objective function;[x1;x2]=initial point;epsilon=algorithm accuracy
% n=max iteration steps; A&b=limit funcs: Ax<=b;Aeq&Beq=limit funcs:Aeqx=Beq
% ib is the lower bound of x1x2
% ub is the upper bound of x1x2
% if there no need to transmit certain parameter, get it=inf;
global result
result=[inf;inf];
global k
k=0;
func2(func,x1,x2,epsilon,n,A,b,Aeq,Beq,ib,ub);
function []=func2(func,x1,x2,epsilon,n,A,b,Aeq,Beq,ib,ub)
%% solve a Linear programming problem
variable_vector=symvar(func);
gradient_func=(subs(gradient(func,[variable_vector(1),variable_vector(2)]),[x1,x2]));
% use matlab internal method linprog()
f=[gradient_func(1,1);gradient_func(1,2)];
%[x,fval]=linprog(f,A,b,Aeq,Beq,lb,ub)
new_f=[double(f(1));double(f(2))];
% pay attention to the optimoptions; if not use it call the error:
% Warning: Your current settings will run a different algorithm 
%('dual-simplex') in a future release.
options=optimoptions('linprog','Algorithm','dual-simplex');
[y,fval]=linprog(new_f,A,b,[],[],ib,[],options);
y1=y(1,1);
y2=y(2,1);
% derive the direction of one-dimension search
d_k=[y1-x1;y2-x2];
z_k=(gradient_func')*d_k;
func3(func,x1,x2,d_k,z_k,epsilon,n,A,b,Aeq,Beq,ib,ub);
function []=func3(func,x1,x2,d_k,z_k,epsilon,n,A,b,Aeq,Beq,ib,ub)
%% terminate func
global result
global k
global temporal_result
global optimality_discrimination_criteria
if k==n
    temporal_result=[x1;x2];
end
if abs(z_k)<=epsilon
    optimality_discrimination_criteria=z_k;
    result=[x1;x2];
else
    func4(func,x1,x2,d_k,epsilon,n,A,b,Aeq,Beq,ib,ub)
end
function []=func4(func,x1,x2,d_k,epsilon,n,A,b,Aeq,Beq,ib,ub)
%% call one-dimension-search(newton-method) to find the best lambda_k
% construct the function of lambda
global newton_result
syms lambda
x1_for_lambda_k=x1+lambda*d_k(1,1);
x2_for_lambda_k=x2+lambda*d_k(2,1);
% variable_vector=symvar(func);
variable_vector=symvar(func);
func_body=subs(func,[variable_vector(1),variable_vector(2)],[x1_for_lambda_k,x2_for_lambda_k]);
% set the initial number of 'bu chang'=lambda0
lambda0=1;
% set the initial number of 'percision'=epsilon0
epsilon0=10.^(-3);
ols_newton(func_body,lambda0,epsilon0)
lambda_k=newton_result;
func5(func,x1,x2,d_k,lambda_k,epsilon,n,A,b,Aeq,Beq,ib,ub);
function []=func5(func,x1,x2,d_k,lambda_k,epsilon,n,A,b,Aeq,Beq,ib,ub)
global k
x1=x1+lambda_k*d_k(1,1);
x2=x2+lambda_k*d_k(2,1);
k=k+1;
func2(func,x1,x2,epsilon,n,A,b,Aeq,Beq,ib,ub);
%% -------------- newton method ------------------------------------

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
    if lambda<=1
        newton_result=lambda;
    else
        % limit search step in [0,1]
        newton_result=1;
    end
else
    newton_func3(func,lambda,epsilon,phai_diff_func_1,phai_diff_func_2);
end

function []=newton_func3(func,lambda,epsilon,phai_diff_func_1,phai_diff_func_2)
% create the new iteration points
lambda=lambda-(phai_diff_func_1/phai_diff_func_2);
newton_func2(func,lambda,epsilon)