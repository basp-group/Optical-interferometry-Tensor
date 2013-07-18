function [ X ] = build_full_tensor( z,mask )
%
% build_full_tensor builds a symmetric tensor from vector of values z and mask.
% 

X=zeros(size(mask));
X(mask==1)=z;

for i=1:size(mask,1)
    for j=1:i
        for k=1:j
            X(i,k,j)=X(i,j,k);
            X(j,i,k)=X(i,j,k);
            X(j,k,i)=X(i,j,k);
            X(k,i,j)=X(i,j,k);
            X(k,j,i)=X(i,j,k);
        end
    end
end


end

