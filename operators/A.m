function [ Xhat ] = A( x,y,z,T )
%
% A performs A(x \otimes y \otimes z) : computes triple products of the
% Fourier transforms of x, y and z at frequencies indicated in table T

n=sqrt(size(x(:),1));m=n;

xhat=1/n*fft2(reshape(x,n,m));xhat=xhat(:);
yhat=1/n*fft2(reshape(y,n,m));yhat=yhat(:);
zhat=1/n*fft2(reshape(z,n,m));zhat=zhat(:);

Xhat=zeros(size(T,1),1);


for i=1:size(T,1)
    
    Xhat(i,1)=xhat(T(i,1))*yhat(T(i,2))*zhat(T(i,3));
    
end
        
        

end

