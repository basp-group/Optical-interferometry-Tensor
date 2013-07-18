function mask = genMask(pdf, seed)
% GENMASK - Generate a mask with variable density sampling
% 


if nargin==2
    rand('seed', seed);
end

%
mask = rand(size(pdf))<pdf;

end