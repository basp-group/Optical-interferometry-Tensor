function [z] = Axy_t( Xhat,xhat,yhat,T )
%
%   Axy_t is the adjoint operator of Axy

n=sqrt(size(xhat,1));m=n;

aux=zeros(size(xhat));

for i=1:size(T,1)
    
    
        aux(T(i,3))=aux(T(i,3))+Xhat(i)*(conj(yhat(T(i,2))*xhat(T(i,1))));
        
    
    
end
        
z=n*ifft2(reshape(aux,n,m));       
z=z(:);
end

