% this function can get the H matrix from the function input
% hence we can get the H matrix from the defination
% H=jacobian(gradient(func))
func=x.^2+2*y.^2;
syms x y t
H1=jacobian(gradient(func));

% *********code below makes no sense***************
% newfunc=subs(func,x*y,t);
% % collect() Sort the function by the power of x
% rectified_func_x=collect(func,x);
% % coeffs() derive the function coeffients by x
% coefficients_x=coeffs(rectified_func_x,x);
% rectified_func_y=collect(func,y);
% coefficients_y=coeffs(rectified_func_y,y);
% rectified_func_t=collect(newfunc,t);
% coefficients_t=coeffs(rectified_func_t,t);
% 
% % if the func has no 'x' we need code below to debug
% if length(coefficients_x)==1
%     coefficients_x=(0);
% end
% 
% % if the func has no 'y' we need code below to debug
% if length(coefficients_y)==1
%     coefficients_y=(0);
% end
% 
% % if the func has no 't=x*y' we need code below to debug
% if length(coefficients_t)==1
%     coefficients_t=(0);
% end
% % form the H matrix by coeffients derived above
% H=[coefficients_x(end),1/2*coefficients_t(end);
%     1/2*coefficients_t(end),coefficients_y(end)];