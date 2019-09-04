function []=gold_ols(func,a,b)
global alpha
alpha=(sqrt(5)-1)/2;
global delta
delta=0.001;
global result
result=-1;

func2(func,a,b);

function []=func2(func,a,b)
global i;
i=1;
global alpha;

lambda=a+(1-alpha)*(b-a);
u=a+alpha*(b-a);

%phai_lambda=feval(func,lambda);
%phai_u=feval(func,u);

func3(func,a,b,lambda,u);

function []=func3(func,a,b,lambda,u)
global delta;
global result;
delta=0.02;
if((b-a)<=delta)
    result=(a+b)/2;
else
    func4(func,a,b,lambda,u);
end

function []=func4(func,a,b,lambda,u)
phai_lambda=feval(func,lambda);
phai_u=feval(func,u);
if(phai_lambda<=phai_u)
    func5(func,a,lambda,u);
else
    func6(func,b,lambda,u);
end

function []=func5(func,a,lambda,u)
new_a=a;
new_b=u;
new_u=lambda;
new_lambda=new_a+(new_b-new_u);
%phai_new_u=feval(func,lambda);
%phai_new_lambda=feval(func,new_lambda);
func7(func,new_a,new_b,new_lambda,new_u);

function []=func6(func,b,lambda,u)
new_a=lambda;
new_b=b;
new_lambda=u;
new_u=new_b-(new_lambda-new_a);
%phai_new_lambda=feval(func,u);
%phai_new_u=feval(func,new_u);
func7(func,new_a,new_b,new_lambda,new_u);

function []=func7(func,new_a,new_b,new_lambda,new_u)
global i;
if(new_u>new_lambda)
    i=i+1;
    a=new_a;
    b=new_b;
    u=new_u;
    lambda=new_lambda;
    func3(func,a,b,lambda,u);
else
    a=new_a;
    b=new_b;
    func2(func,a,b);
end
