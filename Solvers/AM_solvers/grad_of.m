function [ g ] = grad_of( x, y, F, Ft )

% grad_of - computes the gradient of the function .5*||Fx-y||_2^2,
%
% where:
%
%   - F and Ft: Linear operator and its transpose

x=x(:);
N=size(x,1);
n=round(sqrt(N)); m=n;

g=Ft(F(x))-Ft(y);

g=g(:);



end

