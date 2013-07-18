function x = gen_astro_object( paramIM )
%
% gen_astro_object generates an image with the specifications contained in
% structure paramIM
%
% paramIM is a Matlab structure containing the following fields:
%
%   - dim: dimensions of the image
%   - type: type of image to be generated. It admits the following values:
%       -'image' loads image from file (natural image resized)
%       -'spikes' creates a synthetic image with paramIM.k spykes
%   - k: number of spikes in a synthetic image
%   - filename: file name cointaining the image to be generated
%   - F0: zero frequency (sum of the flux)
%


n=paramIM.dim(1);
m=paramIM.dim(2);


if ~isfield(paramIM, 'type'), paramIM.type='spikes';end
if ~isfield(paramIM, 'k'), paramIM.k = round(0.1*n*m); end


if strcmp(paramIM.type,'spikes')
    %sparse synthetic image %%
    x= zeros(n,m);
    k=paramIM.k/(n*m);
    x=vdsmask_full(n,m,k); x=double(x);
    
    x(x==1)=rand(paramIM.k,1);
    
elseif strcmp(paramIM.type,'image')
    %natural resized image %%
    load(paramIM.fileimage);
    x=eval('x');
    x=real(x);
    x(x<0)=0;
    
    
end

x=paramIM.F0*x/sum(sum(x)); % normalisation (flux =1)


end

