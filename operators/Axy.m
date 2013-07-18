function [ Xhat ] = Axy( z,xhat,yhat,T )
%
% Axy operator A_{xy}(x) \equiv A(x \otimes y \otimes z) fixing variables x
% and y
n=sqrt(size(xhat,1));m=n;
zhat=1/n*fft2(reshape(z,n,m));zhat=zhat(:);

Xhat=zeros(size(T,1),1);


for i=1:size(T,1)
    
    Xhat(i,1)=xhat(T(i,1))*yhat(T(i,2))*zhat(T(i,3));
    
end
        
       

end

