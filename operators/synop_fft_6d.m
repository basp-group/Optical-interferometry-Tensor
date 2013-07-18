function T = synop_fft_6d(y, n, m)
%
% adjoint operator of ana3dfft_6d
N=n*m;

T=reshape(y,[n m n m n m]); 


for p = 1:(length(size(T)))
    T = ifft(T,[],p);
    
end

T=reshape(T,[n*m n*m n*m]);
T=T.*(sqrt(n*m))^3;





end

