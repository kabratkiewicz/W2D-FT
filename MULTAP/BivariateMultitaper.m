% Author: Karol Abratkiewicz
% karol.abratkiewicz@pw.edu.pl  
% Warsaw University of technology
% K. Abratkiewicz, "Multitaper ISAR Noise Suppression," in IEEE 
% Transactions on Geoscience and Remote Sensing, vol. 62, pp. 1-13, 2024, 
% Art no. 5217313, doi: 10.1109/TGRS.2024.3427397.

function [W2DFT_original, W2DFT_mean, W2DFT_con] = BivariateMultitaper(signal, M, sigmaT, sigmaR, NFFT_omega, NFFT_eta, avgtype)
% inputs:
% signal - two-dimensional time-domain signal
% M - Hermite function order
% sigmaT - window spread parameter in t
% sigmaR - window spread parameter in r
% NFFT_omega, NFFT_eta - FFT size in the omega and eta frequency domain,
%   respectively
% avgtype - averaging type: 'noncoherent' and 'coherent'
% output:
% W2DFT_hRhT - original windowed bivariate spectrum
% W2DFT_mean - mutitaper concentrated bivariate spectrum mean
% W2DFT_con - concentrated bivariate spectra

if ~exist('avgtype',"var")
    avgtype = 'noncoherent';
end

R = size(signal,2);
T = size(signal,1);

[hT, dhT ,~] = hermf(T, M, sigmaT);
[hR, dhR, ~] = hermf(R, M, sigmaR);

W2DFT_hRhT  = zeros(NFFT_omega, NFFT_eta, M);
W2DFT_dhRhT = zeros(NFFT_omega, NFFT_eta, M);
W2DFT_hRdhT = zeros(NFFT_omega, NFFT_eta, M);

omega_shift = zeros(NFFT_omega, NFFT_eta,M);
eta_shift = zeros(NFFT_omega, NFFT_eta,M);

omega_bins = 1:NFFT_omega;
eta_bins = 1:NFFT_eta;

for i = 1:M
    hRhT   = hT(i,:).'  .* hR(i,:);
    W2DFT_hRhT(:,:,i) = fftshift(fftshift(fft2(signal .* hRhT ,   NFFT_omega, NFFT_eta),1),2);

    dhRhT   = hT(i,:).'  .* dhR(i,:);
    W2DFT_dhRhT(:,:,i) = fftshift(fftshift(fft2(signal .* dhRhT ,   NFFT_omega, NFFT_eta),1),2);

    hRdhT   = dhT(i,:).'  .* hR(i,:);
    W2DFT_hRdhT(:,:,i) = fftshift(fftshift(fft2(signal .* hRdhT ,   NFFT_omega, NFFT_eta),1),2);

    omega_shift(:,:,i) = omega_bins.'   - round(imag(NFFT_omega .* W2DFT_hRdhT(:,:,i)./W2DFT_hRhT(:,:,i)./2/pi));
    eta_shift(:,:,i) = eta_bins - round(imag(NFFT_eta .* W2DFT_dhRhT(:,:,i)./W2DFT_hRhT(:,:,i)./2/pi));
end

W2DFT_con = zeros(size(W2DFT_hRhT));
for k = 1 : M
    for i = 1 : NFFT_omega
        for j = 1 : NFFT_eta
            omega_idx = omega_shift(i, j, k);
            eta_idx = eta_shift(i, j, k);
            if eta_idx < 1 || eta_idx > NFFT_eta
                continue;
            end
            if omega_idx < 1 || omega_idx > NFFT_omega 
                continue;
            end
            if strcmp(avgtype,'coherent')
                W2DFT_con(omega_idx, eta_idx, k) = W2DFT_con(omega_idx, eta_idx, k) + W2DFT_hRhT(i, j, k);
            elseif strcmp(avgtype,'noncoherent')
                W2DFT_con(omega_idx, eta_idx, k) = W2DFT_con(omega_idx, eta_idx, k) + abs(W2DFT_hRhT(i, j, k)).^2;
            end

        end
    end
end
if strcmp(avgtype,'coherent')
    W2DFT_mean =  abs(mean(W2DFT_con, 3)).^2;
elseif strcmp(avgtype,'noncoherent')
    W2DFT_mean =  mean(W2DFT_con, 3);
end

W2DFT_original = W2DFT_hRhT;

end