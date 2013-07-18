function [ Xhat ] = Ayz( x,yhat,zhat,T )
%
% Ayz operator A_{yz}(x) \equiv A(x \otimes y \otimes z) fixing variables y
% and z

n=sqrt(size(yhat,1));m=n;
xhat=1/n*fft2(reshape(x,n,m));xhat=xhat(:);

Xhat=zeros(size(T,1),1);


for i=1:size(T,1)
    
    Xhat(i,1)=xhat(T(i,1))*yhat(T(i,2))*zhat(T(i,3));
    %i=i+1;
end
        
       

end

