function [ yr ] = build_redundant_meas( y, T )
%
% build_redundant_meas - Replicates the measurements (according to table T)


M=size(T,1);

yr=zeros(6*M,1);
iT=1;

for i=1:M
    
    if(T(i,1)>T(i,2) && T(i,2)>T(i,3))
        yr(iT,:)=y(i,:);
        iT=iT+1;
        yr(iT,:)=y(i,:);
        
        iT=iT+1;
        yr(iT,:)=y(i,:);
        
        iT=iT+1;
        yr(iT,:)=y(i,:);
        
        iT=iT+1;
        yr(iT,:)=y(i,:);
        
        iT=iT+1;
        yr(iT,:)=y(i,:);
        iT=iT+1;
        
    elseif (T(i,1)>T(i,2) && T(i,2)==T(i,3))
        yr(iT,:)=y(i,:);
        iT=iT+1;
        yr(iT,:)=y(i,:);
        
        iT=iT+1;
        yr(iT,:)=y(i,:);
        
        iT=iT+1;
        
    elseif (T(i,1)==T(i,2) && T(i,2)>T(i,3))
        yr(iT,:)=y(i,:);
        iT=iT+1;
        yr(iT,:)=y(i,:);
        
        iT=iT+1;
        yr(iT,:)=y(i,:);
        
        iT=iT+1; 
        
    elseif (T(i,1)==T(i,2) && T(i,2)==T(i,3))
        yr(iT,:)=y(i,:);
        iT=iT+1;
    end
    
    yr(all(yr==0,2),:)=[];
    

end

