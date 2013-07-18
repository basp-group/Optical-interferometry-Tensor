function sol = solve_nm(y, param)
%
% Douglas-Rachford splitting algorithm [1] to the Tensor Recovery problem 
% 
% solve_nm - Solve NM problem. 
%
%
%   min ||CX||_*   s.t.  ||y-AX||_2 < epsilon, CX \succeq 0, X \in R_+
%
% INPUTS
%   - y contains the measurements. 
%
%   - param is a Matlab structure containing the following fields:
%
%   General parameters:
% 
%   - verbose: 0 no log, 1 print main steps, 2 print all steps.
%
%   - N: dimension of the signal
%
%   - max_iter: max. nb. of iterations (default: 200).
%
%   - rel_obj: minimum relative change of the solution (default:
%   1e-4)
%       The algorithm stops if
%           | ||x(t)||_2 - ||x(t-1)||_2 | / ||x(t)||_2 < rel_obj,
%       where x(t) is the estimate of the solution at iteration t.
%
%    
%   Projection onto the L2-ball :
%
%   - A: Forward operator // At: Adjoint operator
%
%   - param.tightB2: 1 if A is a tight frame or 0 if not (default = 1)
% 
%   - nuB2: bound on the norm of the operator A, i.e.
%       ||A x||^2 <= nu * ||x||^2 (default: 1)
%
%   
%   Proximal NM operator:
%
%   - C: Forward operator // Ct: Adjoint operator
%
%   - rel_obj_NM: Used as stopping criterion for the proximal NM
%   operator. Min. relative change of the objective value between two
%   successive estimates.
%
%   - max_iter_NM: Used as stopping criterion for the proximal NM
%   operator. Maximun number of iterations.
% 
%   - param.nu_NM: bound on the norm^2 of the operator C, i.e.
%       ||C x||^2 <= nu * ||x||^2 (default: 1)
% 
%   - param.tight_NM: 1 if Ct is a tight frame or 0 if not (default = 1)          
%
%
% OUTPUT
%   - sol: solution of the problem.
%
%
% The problem is solved thanks to a Douglas-Rachford splitting
% algorithm.
% 
% 
%
% References:
% [1] P. L. Combettes and J-C. Pesquet, "A Douglas-Rachford Splitting 
% Approach to Nonsmooth Convex Variational Signal Recovery", IEEE Journal
% of Selected Topics in Signal Processing, vol. 1, no. 4, pp. 564-574, 2007.



% Optional input arguments
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'rel_obj'), param.rel_obj = 1e-4; end
if ~isfield(param, 'max_iter'), param.max_iter = 200; end

% Input parameters for projection onto the L2 ball
param_B2.A = param.A; param_B2.At = param.At;
param_B2.y = y; param_B2.epsilon = param.epsilon;
param_B2.epsilon_up = param.epsilon_up;
param_B2.epsilon_low = param.epsilon_low;
param_B2.verbose = param.verboseB2;
param_B2.nu = param.nuB2;
param_B2.tol=1e-3; %cal? com?
param_B2.max_iter =0; %cal? com?
param_B2.tight = param.tightB2;
param_B2.pos=0;
param_B2.real=0;

% Input parameters for prox NM 
param_NM.pos = 1;
param_NM.real = 1;
param_NM.C=param.C;
param_NM.Ct=param.Ct;
param_NM.lambda=param.lambda_NM;
param_NM.verbose = param.verbose_NM; 

param_NM.rel_obj = param.rel_obj_NM;
param_NM.N=param.N;
param_NM.nu=param.nu_NM;
param_NM.max_iter = param.max_iter_NM;


% Initialization
if isfield(param,'xinit')
    xhat = param.xinit;
else
    xhat = param.At(y); 
end

iter = 1; prev_sol = 0;

% Main loop
while 1
    
    
    if param.verbose >= 1
        fprintf('Iteration %i:\n', iter);
    end
    
    
    % Projection onto the L2-ball
    [sol, param_B2.u] = fast_proj_b2(xhat, param_B2);
    
    % Global stopping criterion
    rel_var=norm(sol-prev_sol)/norm(prev_sol);
    prev_sol=sol;
   

    if param.verbose >= 1
        fprintf('  rel_var = %e\n', ...
             rel_var);
    end
    if (rel_var < param.rel_obj)
        crit_BPDN = 'REL_NORM';
        break;
    elseif iter >= param.max_iter
        crit_BPDN = 'MAX_IT';
        break;
    end
     
    % Proximal NM function and DR update
    xhat = 2*sol - xhat;
    temp = fast_proj_nm(xhat, param_NM);
    xhat = temp + sol - xhat;
    
    % Update variables
    iter = iter + 1;
    
    
end

% Log
if param.verbose>=1
    
    fprintf('\n Solution found:\n');
    
    % Residual
    dummy = param.A(sol); res = norm(y(:)-dummy(:), 2);
    fprintf(' epsilon = %e, ||y-Ax||_2=%e\n', param.epsilon , res);
    
    % Stopping criterion
    fprintf(' %i iterations\n', iter);
    fprintf(' Stopping criterion: %s \n\n', crit_BPDN);
    
end

end