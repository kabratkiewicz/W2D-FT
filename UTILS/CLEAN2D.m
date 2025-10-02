function [comp, residuals] = CLEAN2D(signal, NFFT_omega, NFFT_eta, Ncomp)

T = size(signal,1);
R = size(signal,2);
t = 1:T;
t = t - mean(t);
t = t.';
r = 1:R;
r = r - mean(r);
comp = zeros(T, R, Ncomp);

for i = 1:Ncomp

    F = fft2(signal, NFFT_omega,NFFT_eta);
    omega_scale = fftshift((-NFFT_omega/2:NFFT_omega/2-1)./NFFT_omega);
    eta_scale = fftshift((-NFFT_eta/2:NFFT_eta/2-1)./NFFT_eta);
    [~,I] = max(F,[],"all");
    [dim1, dim2] = ind2sub(size(F),I);
    amp = F(dim1, dim2);
    comp(:,:,i) = exp(1j*2*pi.*t*omega_scale(dim1)) .* exp(1j*2*pi.*r*eta_scale(dim2));
    
    Fcomp = fft2(comp(:,:,i), NFFT_omega,NFFT_eta);
    Fcomp = Fcomp./max(Fcomp(:));
    
    FF = F - (amp)*Fcomp;
    
    signal = ifft2(FF,NFFT_omega,NFFT_eta);
    signal = signal(1:T, 1:R);
    
end
residuals = signal;

end
