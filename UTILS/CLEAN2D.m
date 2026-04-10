function [comp, residuals] = CLEAN2D(signal, NFFT_omega, NFFT_eta, Ncomp)

T = size(signal,1);
R = size(signal,2);

t = (1:T); t = t - mean(t); t = t.';
r = (1:R); r = r - mean(r);

comp = zeros(T, R, Ncomp);

omega_scale = fftshift((-NFFT_omega/2:NFFT_omega/2-1)./NFFT_omega);
eta_scale   = fftshift((-NFFT_eta/2:NFFT_eta/2-1)./NFFT_eta);

for i = 1:Ncomp

    F = fft2(signal, NFFT_omega, NFFT_eta);
    [~,I] = max(abs(F),[],"all");  
    [dim1, dim2] = ind2sub(size(F),I);
    atom = exp(1j*2*pi.*t*omega_scale(dim1)) .* exp(1j*2*pi.*r*eta_scale(dim2));
    atom = atom / norm(atom, 'fro');   %
    amp = sum(signal .* conj(atom), 'all');
    comp(:,:,i) = amp * atom;
    signal = signal - comp(:,:,i);

end

residuals = signal;

end
