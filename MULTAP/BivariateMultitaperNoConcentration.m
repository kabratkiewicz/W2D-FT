% Author: Karol Abratkiewicz
% karol.abratkiewicz@pw.edu.pl  
% Warsaw University of technology
% K. Abratkiewicz, "Multitaper ISAR Noise Suppression," in IEEE 
% Transactions on Geoscience and Remote Sensing, vol. 62, pp. 1-13, 2024, 
% Art no. 5217313, doi: 10.1109/TGRS.2024.3427397.

function [F_wRwT, F_wRwT_mean] = BivariateMultitaperNoConcentration(signal, M, sigmaT, sigmaR, NFFT_omega, NFFT_eta)
% inputs:
% signal - two-dimensional time-domain signal
% M - Hermite function order
% sigmaT - window spread parameter in t
% sigmaR - window spread parameter in r
% NFFT_omega, NFFT_eta - FFT size in the omega and eta frequency domain,
%   respectively
% output:
% F_wRwT - original windowed bivariate spectrum
% F_wRwT_mean - mutitaper bivariate spectrum mean

Y = size(signal,1);
X = size(signal,2);

[hT, ~ ,~] = hermf(X, M, sigmaT);
[hR, ~, ~] = hermf(Y, M, sigmaR);

F_wRwT  = zeros(NFFT_eta, NFFT_omega, M);

for i = 1:M
    wRwT = hT(i,:) .* hR(i,:).';
    F_wRwT(:,:,i) = abs(fftshift(fftshift(fft2(signal .* wRwT ,    NFFT_eta, NFFT_omega),1),2)).^2;
end

F_wRwT_mean =  mean(F_wRwT, 3);

end