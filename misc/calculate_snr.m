function snr = calculate_snr( x,w )

% calculate_snr computes the SNR between x (original image) and w (reconstruction)  

snr=20*log(norm(x,'fro')/norm(x-w, 'fro'));

end

