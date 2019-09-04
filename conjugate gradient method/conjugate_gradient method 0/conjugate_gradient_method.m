% UQP conjugate gradient method
% func takes two variables
% method parameter takes 1 2 3 which stand for Fletcher-Reeves & Dixon &
% Polak-Ribiere-Polyak method
% n is defined as the epochs of iterations
% x1 & x2 is the intial points
% H is derived from func
function []=conjugate_gradient_method(func,method,n,x1,x2,H)
global result
result=inf;
global k
k=1;
x1_new=x1;
% initial x1_old
x1_old=x1;
x2_new=x2;
% initial x2_old
x2_old=x2;
% initial d_k
d_k_old=[inf;inf];
func2(func,method,n,x1_new,x2_new,x1_old,x2_old,H,d_k_old);


function []=func2(func,method,n,x1_new,x2_new,x1_old,x2_old,H,d_k_old)
global result
syms x y
variable_vector=symvar(func);
gradient_func=subs(gradient(func,[variable_vector(1),variable_vector(2)]),[variable_vector(1),variable_vector(2)],[x1_new,x2_new]);
if norm(gradient_func)==0
    result=[x1_new;x2_new];
else
    func3(func,method,n,x1_new,x2_new,x1_old,x2_old,H,d_k_old);
end


function []=func3(func,method,n,x1_new,x2_new,x1_old,x2_old,H,d_k_old)
% d_k_old represents the 'd_k' in the last iteration
% d_k_new represents the 'd_k' in this number'k' iteration
global k
syms x y
variable_vector=symvar(func);
gradient_func=subs(gradient(func,[variable_vector(1),variable_vector(2)]),[variable_vector(1),variable_vector(2)],[x1_new,x2_new]);

% get d_k
if k==1
    d_k_new=-gradient_func;
elseif k>1
    % gey g_k and g_k-1
    gradient_func_old=subs(gradient(func,[variable_vector(1),variable_vector(2)]),[variable_vector(1),variable_vector(2)],[x1_old,x2_old]);
    gradient_func_new=subs(gradient(func,[variable_vector(1),variable_vector(2)]),[variable_vector(1),variable_vector(2)],[x1_new,x2_new]);
    % get beta_k
    if method==1 %FR
        beta_k=norm(gradient_func_new).^2/norm(gradient_func_old).^2;
    elseif method==2 %Dixon
        % beta_k=norm(gradient_func_new).^2/(inv(d_k_old)*gradient_func_old);
        beta_k=-norm(gradient_func_new).^2/((d_k_old')*gradient_func_old);
    elseif method==3 % PRP
        % beta_k=(inv(gradient_func_new)*(gradient_func_new-gradient_func_old))/norm(gradient_func_old).^2;
        beta_k=((gradient_func_new')*(gradient_func_new-gradient_func_old))/norm(gradient_func_old).^2;
    end
    
    d_k_new=-gradient_func+beta_k*d_k_old;
end

func4(func,method,n,x1_new,x2_new,H,d_k_new)


function []=func4(func,method,n,x1_old,x2_old,H,d_k_old)
global k
global result
syms x y
variable_vector=symvar(func);
gradient_func=subs(gradient(func,[variable_vector(1),variable_vector(2)]),[variable_vector(1),variable_vector(2)],[x1_old,x2_old]);
% get lambda_k
lambda_k=-((gradient_func')*d_k_old)/((d_k_old')*H*d_k_old);
x1_new=x1_old+lambda_k*d_k_old(1,1);
x2_new=x2_old+lambda_k*d_k_old(2,1);
if k==n
    result=[x1_new;x2_new];
else
    k=k+1;
    func2(func,method,n,x1_new,x2_new,x1_old,x2_old,H,d_k_old);
end


