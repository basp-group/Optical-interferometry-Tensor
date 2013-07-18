clear;


addpath(genpath('misc'))
addpath(genpath('operators'));
addpath(genpath('Solvers'));
addpath('Data')

%% Generate image example

% input parameters
n= 16; % n rows, image dimension n x m
m= 16; % m columns
N=n*m;
plots=1; %auxiliary flag for plots
noisy=1; %auxiliary flag indicating presence of noise. 1: noisy case; 0: noisless case
input_snr=30; %input snr



% parameters image
% field type:'image' loads image from file (natural image resized) // 'spikes' creates a synthetic image with k spykes
paramIM.dim=[n m];
paramIM.type='image'; 
paramIM.k=16;
paramIM.F0=1;
paramIM.fileimage='eta-carinae_16x16.mat';

% create an image in x
x=gen_astro_object(paramIM); 
if(plots), figure, imagesc(x),colorbar, axis image, end

        

  
%% Generate mask for the open frequences 

u=1; %undersampling ratio u=M/N;
sigmau=1/4;sigmav=1/4;sigmaw=1/4; %standard deviation of the gaussian profile
T=mask_3D_3(n,m,u,sigmau, sigmav, sigmaw);
Tr=build_redundant_table( T );


%% Add noise only to the non-redundant part of the measurements vector
xhat=1/n*fft2(reshape(x,n,m));
ymeas=Ayz(x(:),xhat(:), xhat(:),T);  %true measurements without noise 
NB=numel(ymeas);


if noisy
    % Add Gaussian i.i.d. noise
    
    eB=sqrt(1/NB*sum(abs(ymeas(2:end)).^2));
    sigma_noise=.5*10^(-input_snr/20)*eB;
    
    nBr=(randn(size(ymeas))); 
    nBi=(randn(size(ymeas)));
    
    nBr=sqrt(NB)*sigma_noise*nBr./norm(nBr(:),2);
    nBi=sqrt(NB)*sigma_noise*nBi./norm(nBi(:),2);
    
    nB=nBr + 1i*nBi;
     
    ymeas=ymeas+nB;
    
    
end

y=build_redundant_meas(ymeas,T);
        

%% Parameters of the algorithm

param.n=n;param.m=m;
param.T=Tr;
param.rel_obj=1e-3;
param.max_iter=1000;
param.verbose=1;
        
        
% compute starting estimate for the Lipschitz constant for the backtracking procedure         
aux=rand(size(x));aux=aux./sum(aux(:));
auxhat=1/n*fft2(reshape(aux,n,m));
A=@(x)Ayz(x,auxhat(:),auxhat(:),T);
At=@(x)Ayz_t(x,auxhat(:),auxhat(:),T);

param.L=norm(grad_of( aux, y, A, At ),2)/norm(aux,2);
param.eta=1.1; 
param.max_iter_BT=100;

max_iter=5; %max number of time we reinitialize the problem with rand point
of_vals=zeros(max_iter,1);
sol_vals=zeros(max_iter,n,m);
prev_sol=zeros(N,1);
iter=0;
crit=0;

%% solve AM problem
t_start=tic;

% loop for different initializations of the problem
for iter=1:max_iter

    %initialize x,y,z 
    param.xinit = rand(size(x(:)));param.xinit=param.xinit(:)./sum(param.xinit(:));
    param.yinit=param.xinit;
    param.zinit=param.xinit;


    [sol, of, crit_AM]=solve_am(y,param);


    of_vals(iter)=of;
    sol_vals(iter,:,:)=reshape(sol,n,m);

end
t_end=toc(t_start);

% choose the solution corresponding to the minimum objective function  
[Max, Max_i]=min(of_vals(:));
sol=sol_vals(Max_i,:,:);sol=sol(:);
% projection (real)positive orthant
sol=real(sol); sol(sol<0)=0;

x_AM=reshape(sol,n,m);

%% Log

maxsnr=calculate_snr(x,x_AM);
% plot the result using the same scale as for the original plot
mx=max(x(:));
if (plots>0), figure, imagesc(x_AM,[0,mx]),colorbar, axis image, end 

fprintf(['\n M/N= %e, SNR0 = %e, ellapsed time= %e\n\n'], u, maxsnr, t_end);



    








