function [F_wRwT, wRwT_mean] = BivariateMultitaperNoConcentration(signal, M,sigmaT, sigmaR, NFFT_omega, NFFT_eta)

Y = size(signal,1);
X = size(signal,2);

[hT, ~ ,~] = hermf(X, M, sigmaT);
[hR, ~, ~] = hermf(Y, M, sigmaR);

F_wRwT  = zeros(NFFT_eta, NFFT_omega, M);

for i = 1:M
    wRwT = hT(i,:) .* hR(i,:).';
    F_wRwT(:,:,i) = abs(fftshift(fftshift(fft2(signal .* wRwT ,    NFFT_eta, NFFT_omega),1),2)).^2;
end

wRwT_mean =  mean(F_wRwT, 3);

end