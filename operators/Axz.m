function [ Xhat ] = Axz( y,xhat,zhat,T )
%
% Axz operator A_{xz}(x) \equiv A(x \otimes y \otimes z) fixing variables x
% and z

n=sqrt(size(xhat,1));m=n;
yhat=1/n*fft2(reshape(y,n,m));yhat=yhat(:);

Xhat=zeros(size(T,1),1);


for i=1:size(T,1)
    
    Xhat(i,1)=xhat(T(i,1))*yhat(T(i,2))*zhat(T(i,3));
    %i=i+1;
end
        
       

end

