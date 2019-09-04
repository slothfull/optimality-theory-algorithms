function H=my_hessian(f,x,y,x0)
% matlab has the built-in function hessian
% while my_hessian is built by @sjtu hanzo
switch  nargin
    % nargin returns the number of input arguments passed in the call to 
    % the currently executing function
    case 1
        error('please input enough variables in f(x)')
    case 2
        H=subhessian(f,x);
    case 3
        H=subhessian(f,x);
        H=subs(H,x,x0);
    otherwise
        error('too many arguments or nothing input')
end

% H=hessian(f,x) retrun the hessian form of the func
function HH=subhessian(f,x)
n=length(x);
% jacobian is the matlab built-in function
J=jacobian(f,[x]);
HH=[];
HH=sym(HH);
for i=1:n
HH(i,:)=jacobian(J(1,i),x);
end
