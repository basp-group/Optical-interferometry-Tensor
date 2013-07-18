function [ m ] = nonredundant_mask( N )
%
% nonredundant_mask creates a cubic mask open only for one point of each
% class of a supersymmetric cubic tensor
% 
% INPUT: N dimension of the tensor
% OUTPUT: m mask
%



m=zeros(N,N,N);


for i=1:N
    for j=1:i
        for k=1:j
                m(i,j,k)=1;
        end
    end
end


m(1,1,1)=1; %mask open for the F0 frequency




end

