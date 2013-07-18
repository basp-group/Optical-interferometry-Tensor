function [ Tr ] = build_redundant_table( T )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
M=size(T,1);

Tr=zeros(6*M,3);
iT=1;

for i=1:M
    
    if(T(i,1)>T(i,2) && T(i,2)>T(i,3))
        Tr(iT,:)=T(i,:);
        iT=iT+1;
        Tr(iT,:)=[T(i,1),T(i,3),T(i,2)];
        
        iT=iT+1;
        Tr(iT,:)=[T(i,2),T(i,1),T(i,3)];
        
        iT=iT+1;
        Tr(iT,:)=[T(i,2),T(i,3),T(i,1)];
        
        iT=iT+1;
        Tr(iT,:)=[T(i,3),T(i,1),T(i,2)];
        
        iT=iT+1;
        Tr(iT,:)=[T(i,3),T(i,2),T(i,1)];
        iT=iT+1;
        
    elseif (T(i,1)>T(i,2) && T(i,2)==T(i,3))
        Tr(iT,:)=T(i,:);
        iT=iT+1;
        Tr(iT,:)=[T(i,2),T(i,1),T(i,3)];
        
        iT=iT+1;
        Tr(iT,:)=[T(i,3),T(i,2),T(i,1)];
        
        iT=iT+1;
        
    elseif (T(i,1)==T(i,2) && T(i,2)>T(i,3))
        Tr(iT,:)=T(i,:);
        iT=iT+1;
        Tr(iT,:)=[T(i,1),T(i,3),T(i,2)];
        
        iT=iT+1;
        Tr(iT,:)=[T(i,3),T(i,1),T(i,2)];
        
        iT=iT+1; 
        
    elseif (T(i,1)==T(i,2) && T(i,2)==T(i,3))
        Tr(iT,:)=T(i,:);
        iT=iT+1;
    end
    
    Tr(all(Tr==0,2),:)=[];
    

end

