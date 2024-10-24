% Author: Karol Abratkiewicz
% karol.abratkiewicz@pw.edu.pl  
% Warsaw University of technology
% K. Abratkiewicz, "Multitaper ISAR Noise Suppression," in IEEE 
% Transactions on Geoscience and Remote Sensing, vol. 62, pp. 1-13, 2024, 
% Art no. 5217313, doi: 10.1109/TGRS.2024.3427397.

function [W2DFT_original, W2DFT_original_mean] = BivariateMultitaperNoConcentration(signal, M, sigmaT, sigmaR, NFFT_omega, NFFT_eta)
% inputs:
% signal - two-dimensional time-domain signal
% M - Hermite function order
% sigmaT - window spread parameter in t
% sigmaR - window spread parameter in r
% NFFT_omega, NFFT_eta - FFT size in the omega and eta frequency domain,
%   respectively
% output:
% F_hRhT - original windowed bivariate spectrum
% F_hRhT_mean - mutitaper bivariate spectrum mean

R = size(signal,2);
T = size(signal,1);

[hT, ~ ,~] = hermf(T, M, sigmaT);
[hR, ~, ~] = hermf(R, M, sigmaR);

W2DFT_hRhT  = zeros(NFFT_omega, NFFT_eta, M);

for i = 1:M
    hRhT = hT(i,:).' .* hR(i,:);
    W2DFT_hRhT(:,:,i) = abs(fftshift(fftshift(fft2(signal .* hRhT ,    NFFT_omega, NFFT_eta),1),2)).^2;
end

W2DFT_original_mean =  mean(W2DFT_hRhT, 3);
W2DFT_original = W2DFT_hRhT;

end