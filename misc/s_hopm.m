
function  u  = s_hopm( T, N, tolu )
%
%   s_hopm computes the best rank-1 tensor approximation in the LS sense using the 
%   Symmetric Higher-Order Power Method (S-HOPM) described in [1]
%
%   INPUTS
%
%   - T: Initial Tensor unfolded (cubic and supersymmetric)
%   - N: dimension of one side of the cubic tensor T
%   - tolu: tolerance on the relative variation of the solution. 
%       The algorithm stops if
%           | ||u(t)||_2 - ||u(t-1)||_2 | / ||u(t)||_2 < tolu,
%       where u(t) is the estimate of the solution at iteration t.
%
%   OUTPUT
%   
%   - u: underpinning vector of the rank-1 tensor approximation of T.
%   
% [1] E. Kofidis and P.A. Regalia, "On the best rank-1 approximation of
% higher-order supersymmetric tensors", 2002, SIAM J. Matrix Anal. Appl,23, 863
%


% compute initial point through higer order svd
[S U sv tol] = hosvd(reshape(T,N,N,N));
u=U{1}(:,1);
u=u./norm(u(:),2);
uold=u;

if nargin < 3
	tolu = 1e-3;
end
iterS_HOPM=1;

%main loop
while(1)
    
    T1=reshape(T,N,N*N);
    u=T1*kron(u,u);
    u=u./norm(u(:),2);
   
    if (norm(uold(:)-u(:))/norm(uold(:)) < tolu)
        break;
    end
    uold=u;
    iterS_HOPM=iterS_HOPM+1;
end

end

