function [sol, u] = fast_proj_nm(x, param)
% 
%
% sol = fast_proj_NM(x, param) solves:
%
%   min_{Z} ||CZ||_* + 1/2*||Z-X||^2 s.t. CZ \succeq 0 and Z \in R_+
%
%   
%
% param a Matlab structure containing the following fields:
%  
%
%   - C: Forward operator (default: Id).
%
%   - Ct: Adjoint operator (default: Id).
%
%   - tight: 1 if C is a tight frame or 0 if not (default = 1)
%
%   - nu: bound on the norm of the operator C, i.e.
%       ||C x||^2 <= nu * ||x||^2 (default: 1)
%
%   - lambda: parameter used for the soft-thresholding operator
%
%   - pos: positivity flag 
%
%   - real: reality flag
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

N=param.N;


if ~isfield(param, 'tight'), param.tight = 0; end
if ~isfield(param, 'rel_obj'), param.rel_obj = 1e-3; end
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'nu'), param.nu = 1; end
if ~isfield(param, 'max_iter'), param.max_iter = 200; end
if ~isfield(param, 'u'), param.u = zeros(N*N,1); end
if ~isfield(param, 'pos'), param.pos = 0; end
if ~isfield(param, 'real'), param.real = 0; end



% Useful functions for the projection
sc = @(z) z*min(param.epsilon/norm(z(:)), 1); % scaling

% Compute proximal operator
if (param.tight && ~(param.pos||param.real)) % TIGHT FRAME CASE
    
    temp = param.C(x,N) - param.y;
    sol = x + 1/param.nu * param.Ct(sc(temp)-temp);
    crit_NM = 'TOL_EPS'; iter = 0;
    u = 0;
    
else % NON TIGHT FRAME CASE
    
    % Initializations
    sol = x; u = param.u; v = u;
    iter = 1; true = 1; told = 1;
    of=0; 
   
    while true
        
        % Current estimate
        dummy=x-param.Ct(u);
        
        % projection (real)positive orthant
        sol=real(dummy);
        sol(sol<0)=0;
        
        % Evaluation Objective function
        of_old=of;
        
        [V,Sigma] = svd(reshape(param.C(sol),N,N));
        aux=diag(Sigma); %vector with the singular values
        of=sum(aux(:)) + 0.5*norm(sol(:)-x(:),'fro')^2;
        
        rel_var=abs(of - of_old)/abs(of_old);
        
        % Stopping criterion
        if (rel_var<param.rel_obj)
            crit_NM = 'TOL_OBJ_FUNC'; break;
        elseif iter >= param.max_iter
            crit_NM = 'MAX_IT'; break;
        end
        
        % FISTA update
        t = (1+sqrt(1+4*told^2))/2;
        u = v; 
        z=param.nu*u + param.C(sol);
        
        %prox nuclear norm (+ PSD constraint)
        v=1/param.nu*(z - soft_sgv(z,param.N,param.nu*param.lambda));
        u = v + (told-1)/t * (v - u);
        
        
        % Update number of iteration
        told = t;
        
        % Update number of iterations        
        iter = iter + 1;
        
    end
end

% Log 
if param.verbose >= 1
    
    fprintf(['  Proj. NN: obj func= %e, , crit=%s, iter = %i\n'], of, crit_NM, iter);
end


end