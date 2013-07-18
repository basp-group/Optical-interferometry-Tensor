function [ T ] = mask_3D_3( n,m,p, sigmau, sigmav, sigmaw )

% mask_3D_3 generates a sampling pattern of the 3D-cube
%
% The sampling pattern is obtained by sampling frequencies along each of the 3 
% tensor dimensions from a bidimensional Gaussian profile in the corresponding Fourier plane. 
% The procedure is repeated until M samples are obtained. 
% Note that the originally continuous frequencies are associated with their nearest discrete neighbour, 
% and if a product is sampled twice the result is discarded.
%
% INPUTS
%
%   N=n*m dimension of the signal
%   p: M/N (ratio measurements/dim signal) required
%   sigmau,sigmav,sigmaw: correspond to the stdev of the Gaussian profile for
%   the 3 dimensions
%
% OUTPUT
%
%   T table M rows (number of sampled triplets), 3 columns (each
%   frequency of the triplet)




nmeas=n*m*p;

[x,y] = meshgrid(linspace(-1, 1, m), linspace(-1, 1, n)); % Cartesian grid
r = sqrt(x.^2+y.^2); r = r/max(r(:)); % Polar grid


meas=0; mask=zeros(n,m);
T=zeros(nmeas,3);

while meas<nmeas

    u=randn(1,2)*sigmau;
    
    % find the nearest neighbor on the discrete grid
    d=2;idx_min=0;
    for idx=1:n
        if abs(u(1,1)-x(1,idx))< d
            idx_min=idx; d=abs(u(1,1)-x(1,idx));
        end
    end
    ux=idx_min;

    d=2;idy_min=0;
    for idy=1:n
        if abs(u(1,2)-y(idy,1))< d
            idy_min=idy; d=abs(u(1,2)-y(idy,1));
        end
    end
    uy=idy_min;
    
    % find frequency index
    mask(ux,uy)=1; mask=ifftshift(mask);
    u=find(mask(:));
    mask=zeros(n,m);
    
    v=randn(1,2)*sigmav;
    
    % find the nearest neighbor on the discrete grid
    d=2;idx_min=0;
    for idx=1:n
        if abs(v(1,1)-x(1,idx))< d
            idx_min=idx; d=abs(v(1,1)-x(1,idx));
        end
    end
    vx=idx_min;

    d=2;idy_min=0;
    for idy=1:n
        if abs(v(1,2)-y(idy,1))< d
            idy_min=idy; d=abs(v(1,2)-y(idy,1));
        end
    end
    vy=idy_min;
    
    % find frequency index
    mask(vx,vy)=1; mask=ifftshift(mask);
    v=find(mask(:));
    mask=zeros(n,m);
    
    w=randn(1,2)*sigmaw;

    % find the nearest neighbor on the discrete grid
    d=2;idx_min=0;
    for idx=1:n
        if abs(w(1,1)-x(1,idx))< d
            idx_min=idx; d=abs(w(1,1)-x(1,idx));
        end
    end
    wx=idx_min;

    d=2;idy_min=0;
    for idy=1:n
        if abs(w(1,2)-y(idy,1))< d
            idy_min=idy; d=abs(w(1,2)-y(idy,1));
        end
    end
    wy=idy_min;
    
    % find frequency index
    mask(wx,wy)=1; mask=ifftshift(mask);
    w=find(mask(:));
    mask=zeros(n,m);
    
    
    V=sort([u,v,w],'descend'); % predefined descending order
    
    if ismember(T,V,'rows')==zeros(size(V,1))
        T(meas+1,:)=V;
        meas=meas+1;
    end
    
end
    
    if isempty(find(all(T==1,2)))
        T(meas+1,:)=[1,1,1]; % open the mask for the F0 frequency
    end



end

