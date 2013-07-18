function [sol, u] = fast_proj_b2(x, param)
% PROJ_B2 - Projection onto a L2-ball
%
% sol = fast_proj_b2(x, param) solves:
%
%   min_{z} ||x - z||_2^2   s.t.  ||y - A z||_2 < epsilon
%
% param is a Matlab structure containing the following fields:
%
%   - y: measurements (default: 0).
%
%   - A: Forward operator (default: Id).
%
%   - At: Adjoint operator (default: Id).
%
%   - epsilon: Radius of the L2 ball (default = 1e-3).
%
%   - tight: 1 if A is a tight frame or 0 if not (default = 1)
%
%   - nu: bound on the norm of the operator A, i.e.
%       ||A x||^2 <= nu * ||x||^2 (default: 1)
%
%   - tol: tolerance for the projection onto the L2 ball. The algorithms
%   stops if
%       epsilon/(1-tol) <= ||y - A z||_2 <= epsilon/(1+tol)
%   (default: 1e-3)
%
%   - max_iter: max. nb. of iterations (default: 200).
%
%   - verbose: 0 no log, 1 a summary at convergence, 2 print main
%   steps (default: 1)
%
%
%
% References:
% [1] M.J. Fadili and J-L. Starck, "Monotone operator splitting for
% optimization problems in sparse recovery" , IEEE ICIP, Cairo,
% Egypt, 2009.
% [2] Amir Beck and Marc Teboulle, "A Fast Iterative Shrinkage-Thresholding
% Algorithm for Linear Inverse Problems",  SIAM Journal on Imaging Sciences
% 2 (2009), no. 1, 183--202.
%


% Optional input arguments
if ~isfield(param, 'y'), param.y = 0; end
if ~isfield(param, 'A'), param.A = @(x) x; param.At = @(x) x; end
if ~isfield(param, 'epsilon'), param.epsilon = 1e-3; end
if ~isfield(param, 'tight'), param.tight = 1; end
if ~isfield(param, 'tol'), param.tol = 1e-3; end
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'nu'), param.nu = 1; end
if ~isfield(param, 'max_iter'), param.max_iter = 200; end
if ~isfield(param, 'u'), param.u = zeros(size(param.y)); end
if ~isfield(param, 'pos'), param.pos = 0; end
if ~isfield(param, 'real'), param.real = 0; end

% Useful functions for the projection
sc = @(z) z*min(param.epsilon/norm(z(:)), 1); % scaling

% Projection

% TIGHT FRAME CASE
% In NM algorithm this is the used case, the general (non-tight) case is
% included for completeness
if (param.tight && ~(param.pos||param.real))
    
    temp = param.A(x) - param.y;
    sol = x + 1/param.nu * param.At(sc(temp)-temp);
    crit_B2 = 'TOL_EPS'; iter = 0;
    u = 0;
    
else % NON TIGHT FRAME CASE
    
    % Initializations
    sol = x; u = param.u; v = u;
    iter = 1; true = 1; told = 1;
    
    % Tolerance onto the L2 ball
    epsilon_low = param.epsilon_low;
    epsilon_up = param.epsilon_up;
    
    % Check if we are in the L2 ball
    dummy = param.A(sol);
    norm_res = norm(param.y(:)-dummy(:), 2);
    if norm_res <= epsilon_up
        crit_B2 = 'TOL_EPS'; true = 0;
    end
    
    % Projection onto the L2-ball
    % Init
    if param.verbose > 1
        fprintf('  Proj. B2:\n');
    end
    while true
        
        % Residual
        res = param.A(sol) - param.y; norm_res = norm(res(:), 2);
        
        % Scaling for the projection
        res = u*param.nu + res; norm_proj = norm(res(:), 2);
        
        % Log
        if param.verbose>1
            fprintf('   Iter %i, epsilon = %e, ||y - Ax||_2 = %e\n', ...
                iter, param.epsilon, norm_res);
        end
        
        % Stopping criterion
        if (norm_res>=epsilon_low && norm_res<=epsilon_up)
            crit_B2 = 'TOL_EPS'; break;
        elseif iter >= param.max_iter
            crit_B2 = 'MAX_IT'; break;
        end
        
        % Projection onto the L2 ball and FISTA update
        t = (1+sqrt(1+4*told^2))/2;
        ratio = min(1, param.epsilon/norm_proj);
        u = v;
        v = 1/param.nu * (res - res*ratio);
        u = v + (told-1)/t * (v - u);
        
        % Current estimate
        sol = x - param.At(u);
        
        %Projection onto the non-negative orthant (positivity constraint)
        
        if (param.pos)
            sol=real(sol);
            sol(sol<0)=0;
            
        end
        
        %Projection onto the real orthant (reality constraint)
        
        if (param.real)
            sol=real(sol);
            
        end
        
        
        
        % Update number of iteration
        told = t;
        
        % Update number of iterations        
        iter = iter + 1;
        
    end
end

% Log after the projection onto the L2-ball
if param.verbose >= 1
    temp = param.A(sol);
    fprintf(['  Proj. B2: epsilon = %e, ||y-Ax||_2 = %e,', ...
        ' %s, iter = %i\n'], param.epsilon, norm(param.y(:)-temp(:)), ...
        crit_B2, iter);
end


end