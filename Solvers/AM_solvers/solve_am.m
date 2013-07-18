function [sol, of, crit_AM]=solve_am(ymeas, param)
% 
%
% solve_AM - Solves the following minimization problem      
%
%       min_{x,y,z} ||A(x \otimes y \otimes z)||_2^2 s.t. sol \in R_+,
%
% using Alternate Minimization algorithm, as a generalized version of algorithm presented in [1] 
%       
%   INPUT
%       - ymeas: vector of measurements
%       - param: data structure with parameters of the algorithm
%
%   OUTPUT
%       - sol: local solution
%       - of: objective function
%       - crit_AM: stoping criteria used
%
%   param is a Matlab structure containing:
%
%   General parameters:
% 
%   - verbose: 0 no log, 1 print main steps, 2 print all steps.
%
%   - n,m: dimensions of the sought image
%
%   - T: table indicating the probed triplets of frequencies
%
%   - max_iter: max. nb. of iterations (default: 200).
%
%   - rel_obj: minimum relative change of the solution (default:
%   1e-4)
%       The algorithm stops if
%           | ||x(t)||_2 - ||x(t-1)||_2 | / ||x(t)||_2 < rel_obj and
%           | ||y(t)||_2 - ||y(t-1)||_2 | / ||y(t)||_2 < rel_obj and
%           | ||z(t)||_2 - ||z(t-1)||_2 | / ||z(t)||_2 < rel_obj, 
%       where x(t), y(t), z(t) are the partial estimates of the solution at iteration t.
%
%   Parameters for the backtracking procedure in Forward-Backward Algorithm
%
%   - L: initial approximation Lipsitsz constant
%
%   - eta: tuning parameter (to find the optimal step size in the gradient step)
%
%   - max_iter_BT: max. nb. of iterations (backtracking)
%
%
%   [1] J.P. Haldar and D. Hernando, "Rank-Constrained Solutions to Linear Matrix
%   Equations Using PowerFactorization", 2009, IEEE Signal Processing Letters, 16, 584
%



if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'rel_obj'), param.rel_obj = 1e-3; end
if ~isfield(param, 'max_iter'), param.max_iter = 200; end


n=param.n;
m=param.m;

iter=0;prev_x=0;prev_y=0;prev_z=0;

while 1
    
%% define Ayz operators to solve for x
% compute Fourier trnasforms y and z
yhat=1/n*fft2(reshape(param.yinit,n,m));
zhat=1/n*fft2(reshape(param.zinit,n,m));

param.A=@(x)Ayz(x,yhat(:),zhat(:),param.T);
param.At=@(x)Ayz_t(x,yhat(:),zhat(:),param.T);

param.solini=param.xinit(:);


% solve problem for x

[x, crit_ncFBx]=solve_fb_bt(ymeas,param);

% update initial solution
param.xinit=x(:);


%% define Axz operators to solve for y
% compute Fourier trnasforms x and z
xhat=1/n*fft2(reshape(param.xinit,n,m));
zhat=1/n*fft2(reshape(param.zinit,n,m));
param.A=@(x)Axz(x,xhat(:),zhat(:),param.T);
param.At=@(x)Axz_t(x,xhat(:),zhat(:),param.T);


param.solini=param.yinit(:);

% solve problem for y

[y,crit_ncFBy]=solve_fb_bt(ymeas,param);

% update initial solution
param.yinit=y(:);


%% define Axy operators to solve for z
% compute Fourier trnasforms x and y
xhat=1/n*fft2(reshape(param.xinit,n,m));
yhat=1/n*fft2(reshape(param.yinit,n,m));
param.A=@(x)Axy(x,xhat(:),yhat(:),param.T);
param.At=@(x)Axy_t(x,xhat(:),yhat(:),param.T);

param.solini=param.zinit(:);

%solve problem for z
[z, crit_ncFBz]=solve_fb_bt(ymeas,param);

% update initial solution
param.zinit=z(:);



%% Global stopping criterion
    
       
rel_varx = norm(x - prev_x)/norm(prev_x);
rel_vary = norm(y - prev_y)/norm(prev_y);
rel_varz = norm(z - prev_z)/norm(prev_z);

% compute objective
rel_error= norm(A(x,y,z,param.T)-ymeas)/norm(ymeas);
of=.5*norm(A(x,y,z,param.T)-ymeas)^2;

% final solution
sol=1/3*(x+y+z);

if param.verbose >= 1
    fprintf('  ||Ax-b||_2 = %e, rel_var = %e\n', ...
        norm(A(x,y,z,param.T)-ymeas), rel_varx);
end

if (rel_varx < param.rel_obj && rel_vary < param.rel_obj && rel_varz < param.rel_obj || norm(x-prev_x)==0)
    crit_AM = 'REL_NORM';
    break;
elseif (rel_error < 1e-3)
    crit_AM = 'REL_ERROR';
    break;
elseif iter >= param.max_iter
    crit_AM = 'MAX_IT';
    break;
    
elseif strcmp(crit_ncFBx,'MAX_IT') || strcmp(crit_ncFBy,'MAX_IT') ||strcmp(crit_ncFBz,'MAX_IT')
    crit_AM = 'MAX_IT';
end
    
    
    
% Update variables
iter = iter + 1;
prev_x = x;
prev_y = y;
prev_z = z;

end

end

