function [ X ] = fast_tensor_product( x )
%
%   fast_tensor_product - computes the tensor product of vector x with
%   itself
%
%   OUTPUT: X = x \otimes x \otimes x

N=size(x,1);
x1=x.'; x1=x1(ones(1,N),:,ones(1,N)); %equivalent to x1=repmat(x.',[N,1,N]);

x2=x(:,ones(1,N), ones(1,N)); %equivalent x2=repmat(x,[1,N,N]);

x3=zeros(1,1,N);
x3(1,1,:)=x;

x3=x3(ones(1,N),ones(1,N),:);%equivalent x3=repmat(x3,[N,N,1]);

X=x1.*x2.*x3;

end

