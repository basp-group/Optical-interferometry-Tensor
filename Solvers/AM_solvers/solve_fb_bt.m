function [sol, crit_ncFB] = solve_fb_bt(y, param)

% solve_fb_bt - solve the following problem
%
% Solves min ||A*sol - y||_2^2 s.t. sol \in R_+,
%
%   based on the FAST ITERATIVE SHRINKAGE-THRESHOLDING ALGORITHM (with
%   a backtracking line search procedure) described in [1]
%
%   INPUT
%       - y: vector of measurements
%       - param: data structure with parameters of the algorithm
%
%   OUTPUT
%       - sol: local solution
%       - crit_ncFB: stoping criteria used
%
%   param is a Matlab structure containing:
%
%   General parameters:
% 
%   - verbose: 0 no log, 1 print main steps, 2 print all steps.
%
%   - n,m: dimensions of the sought image
%
%   - A and At: operator and its adjoint
%
%   - T: table indicating the probed triplets of frequencies
%
%   - max_iter: max. nb. of iterations (default: 200).
%
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
% References
% [1] A. Beck and M. Teboulle, "A Fast Iterative Shrinkage-Thresholding Algorithm for Linear Inverse
% Problems"
% 2009 Society for Industrial and Applied Mathematics Vol. 2, No. 1, pp. 183?202




% Optional input arguments
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'rel_obj'), param.rel_obj = 1e-3; end
if ~isfield(param, 'max_iter'), param.max_iter = 200; end



sol=param.solini;
iter = 0; prev_sol = 0;


% Main loop
while 1
    
    %
    if param.verbose >= 1
        fprintf('Iteration %i:\n', iter);
    end
    
    iter_BT=0;i=0;
    
    %compute gradient objective function
    g=grad_of(sol,y,param.A,param.At);
    
   
    
    % backtracking
    while 1
        
        iter_BT=iter_BT+1;
        i=i+1;
        param.Lhat=(param.eta)^i*param.L;
        param.gamma=1/param.Lhat;
        
        %gradient descent step
        dummy=sol-param.gamma*g;
        
        
        %projection (real)positive quadrant
        dummy=real(dummy);
        dummy(dummy<0)=0;
    
        e_norm=0.5*norm(param.A(dummy)-y).^2;
        of=e_norm;
        
        G=.5*norm(param.A(sol)-y,2)^2 + dot(dummy(:)-sol(:),real(g)) +.5*param.Lhat*norm(dummy-sol,2)^2;
        
        % backtracking criteria
        if of< G
            crit_BT = 'CONV';
            break; 
        end
        if iter_BT>param.max_iter_BT
            crit_BT = 'MAX_ITER';
            break;
        end
        
    end
   
    
    if param.verbose >= 1 && strcmp(crit_BT,'MAX_ITER')
        fprintf(' Stopping criterion BT: %s \n\n', crit_BT);
    end
    
    % Update variables
    sol=dummy;
    param.L=param.Lhat;
        
    rel_var = norm(sol - prev_sol)/norm(prev_sol);
    
    
    if param.verbose >= 1
        fprintf('  0.5*(||Ax-b||_2)^2 = %e, rel_var = %e\n', ...
            e_norm, rel_var);
    end
    
    if (rel_var < param.rel_obj || norm(sol-prev_sol)==0)
        crit_ncFB = 'REL_NORM';
        break;
    elseif iter >= param.max_iter
        crit_ncFB = 'MAX_IT';
        break;
    end
    
    % Update variables
    iter = iter + 1;
    prev_sol = sol;
      
end

% Log
if param.verbose>=1
    
    fprintf('\n Solution found:\n');
    
    % Residual
    dummy = param.A(sol); res = norm(y(:)-dummy(:), 2);
    fprintf('||y-Ax||_2=%e\n', res);
    
    % Stopping criterion
    fprintf(' %i iterations\n', iter);
    fprintf(' Stopping criterion: %s \n\n', crit_ncFB);
    
end

end