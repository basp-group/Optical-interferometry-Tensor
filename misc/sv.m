function w = sv( sol,N )

% sv computes the eigenvalues of matrix C(sol)

sol=reshape(sol,N,N,N);
Z1=C(sol,N);
Z1=reshape(Z1,N,N);
[V D] = eig(Z1);
w=(diag(D));
w(w<0)=0;


end

