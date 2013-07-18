function [p,r]=soft_sgv(z, N, lambda)
% soft_sgv computes the singular value soft-threasholding and PSD
% projection
%
%   INPUTS
%   
%   - N: dimension of the squared matrix
%   - z: matrix (unfolded)
%   - lambda: soft-thresholding parameter
%
%   OUTPUTS
%   
%   - p: low-rank approximation of the matrix (unfolded)
%   - r: rank of the low-rank approximation matrix
%

%reshape in matrix form
z = reshape(z, [], N);

% symmetrization for rounding errors
z=(z+z.')/2;

[V,Sigma] = eig(z);
aux=diag(Sigma);


% soft-thresholding
aux=max(aux-lambda,0);
% PSD projection
aux(aux<0)=0;
% rank
r=length(find(aux~=0));


%build low-rank approximation
Sigma=diag(aux);
p = V * Sigma * V';

%unfold solution
p = p(:);
end

