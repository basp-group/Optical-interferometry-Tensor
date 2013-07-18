function T = analop_fft_6d(T, n, m)
% 
% analop_fft_6d computes 2D-Fourier Transform along the 3 dimensions of cubic tensor T


N=size(T,1);


T=reshape(T, [n m n m n m]);

for p = 1:(length(size(T)))
    T = fft(T,[],p);
    
end


T=T./(sqrt(n*m))^3;

T=reshape(T, [n*m n*m n*m]);

end

