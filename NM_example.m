clear;


addpath(genpath('/tptool')) %change this path to the Tensor Product Toolbox path in your computer
addpath('operators')
addpath('Data')
addpath(genpath('misc'));
addpath(genpath('Solvers'));


%% Generate image example

% input parameters
n= 8; % n rows, image dimension n x m
m= 8; % m columns
N=n*m;
plots=1; %auxiliary flag for plots
noisy=1; %auxiliary flag indicating presence of noise. 1: noisy case; 0: noisless case
input_snr=30; %input snr

% parameters image
% field type:'image' loads image from file (natural image resized) // 'spikes' creates a synthetic image with k spykes
paramIM.dim=[n m];
paramIM.type='spikes'; 
paramIM.k=6;
paramIM.F0=1;
paramIM.fileimage='eta-carinae_16x16.mat';

% create an image in x
x=gen_astro_object(paramIM); 
if(plots), figure, imagesc(x),colorbar, axis image, end
xu=x(:); %unfolded signal
Tx=fast_tensor_product(xu); % outer product of the image vector -> Tx is a (n*m, n*m, n*m) tensor

  
%% Generate mask for the open frequences
u=1;
mask=nonredundant_mask(N);

sigmau=1/4;sigmav=1/4;sigmaw=1/4;
T=mask_3D_3(n,m,u,sigmau, sigmav, sigmaw);
Tr=build_redundant_table( T );

maskB=zeros(N,N,N);
for i=1:size(Tr,1)
  I=Tr(i,:);
  maskB(I(1),I(2),I(3))=1;
end



%% Add noise only to the non-redundant part of the measurements vector

y0=analop_fft_6d(Tx,n,m); %compute 3-dimension 2D-Fourier Transform
yB=y0(maskB&mask); %measured values, only once
NB=numel(yB);



idmask=find(maskB&mask==1);

if noisy
    % Add Gaussian i.i.d. noise
    
    eB=sqrt(1/NB*sum(abs(yB(2:end)).^2));
    param.sigma_noise=.5*10^(-input_snr/20)*eB;
    
    nBr=(randn(size(yB))); 
    nBi=(randn(size(yB)));
    
    nBr=sqrt(NB)*param.sigma_noise*nBr./norm(nBr(:),2);
    nBi=sqrt(NB)*param.sigma_noise*nBi./norm(nBi(:),2);
    
    nB=nBr + 1i*nBi;
    yB=yB+nB;
    
    % noise bound based on a chi-squared model
    param.epsilon=sqrt(6)*sqrt(2*NB + 4*sqrt(NB))*param.sigma_noise ;
    % tolerance on the noise bound
    param.epsilon_up= sqrt(6)*(sqrt(param.sigma_noise^2*(2*NB) + 5*param.sigma_noise^2*sqrt(2*(2*NB))));%expectation + 2.1 std
    param.epsilon_low= sqrt(6)*(sqrt(param.sigma_noise^2*(2*NB)  -5*param.sigma_noise^2*sqrt(2*(2*NB))));
else
    param.epsilon=1e-7; %default 1e-7 or less, depending on the scale of the signal (ideally should be 0)
end

% build measurement vector
maskB=logical(maskB);
Ty=build_full_tensor(yB,maskB&mask); % symmetrization of the tensor (replicate measurements)
y=Ty(maskB);
y=y(:);

  
%% Parameters of the algorithm

% general parameters
param.max_iter = 200;
param.rel_obj=5e-3;
param.verbose = 2;
param.N=N;

% specfic parameters L2-ball projection  
param.A = @(x)analop3dfft(x,n,m, maskB); 
param.At = @(x)synop3dfft(x,n,m, maskB);
param.nuB2=1;
param.tightB2=1;
param.verboseB2=2;
 
% specific parameters NN-prox 
param.C=@(x)C(x,N);
param.Ct=@(x)Ct(x,N);
param.nu_NM=N;
param.tight_NM=1;
param.max_iter_NM=200;
param.rel_obj_NM=5e-3;
param.verbose_NM=2;

%estimate initial point  
Aty=param.At(y);
param.xinit = Aty;
%compute threshold parameter for the NN-prox
aux=sv(Aty,N); %docu!
param.lambda_NM=0.1*max(abs(aux))/param.nu_NM; 



%% Solve NM problem
t_start=tic;

sol=solve_nm(y,param);

% projection (real)positive orthant
sol=real(sol); sol(sol<0)=0;

% rank 1 approximation 
[ x_NM] = s_hopm( sol, N, 1e-3 );
x_NM=x_NM./sum(x_NM(:)); %normalisation
x_NM=reshape(x_NM,n,m);

t_end=toc(t_start);

%% Log

maxsnr=calculate_snr(x,x_NM); 

% plot the result using the same scale as for the original plot
mx=max(x(:));
if (plots), figure, imagesc(x_NM, [0,mx]),colorbar, axis image,  end 

fprintf('\n M/N= %e, SNR0 = %e, Ellapsed time= %e\n\n', u, maxsnr, t_end);





    
