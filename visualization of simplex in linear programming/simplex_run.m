% problem_1
% f=[1 -2]';
% A=[-1 -1;1 -1;0 1];
% b=[-2 -1 3]';
% [dcxb,x,fval,exitflag,flag]=simplex(f,A,b);

% problem_2
% for the problem below it can not get the right answer!
f=[-1 3 1]';
A=[3 -1 2;-2 4 0;-4 3 8];
b=[7 12 10]';
[dcxb,x,fval,exitflag,flag]=simplex(f,A,b);
