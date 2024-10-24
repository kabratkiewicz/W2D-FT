% Author: Karol Abratkiewicz
% karol.abratkiewicz@pw.edu.pl  
% Warsaw University of Technology
% K. Abratkiewicz, "Multitaper ISAR Noise Suppression," in IEEE 
% Transactions on Geoscience and Remote Sensing, vol. 62, pp. 1-13, 2024, 
% Art no. 5217313, doi: 10.1109/TGRS.2024.3427397.

function [W2DFT_wRwT, W2DFT_mean, W2DFT_con] = BivariateMultitaper(signal, M, sigmaT, sigmaR, NFFT_omega, NFFT_eta, avgtype)
% inputs:
% signal - two-dimensional time-domain signal
% M - Hermite function order
% sigmaT - window spread parameter in t
% sigmaR - window spread parameter in r
% NFFT_omega, NFFT_eta - FFT size in the omega and eta frequency domain,
%   respectively
% avgtype - averaging type: 'noncoherent' and 'coherent'
% output:
% W2DFT_wRwT - original windowed bivariate spectrum
% W2DFT_mean - mutitaper concentrated bivariate spectrum mean
% W2DFT_con - concentrated bivariate spectra

if ~exist('avgtype',"var")
    avgtype = 'noncoherent';
end

R = size(signal,2);
T = size(signal,1);

[hT, dhT ,~] = hermf(T, M, sigmaT);
[hR, dhR, ~] = hermf(R, M, sigmaR);

W2DFT_wRwT  = zeros(NFFT_omega, NFFT_eta, M);
W2DFT_dwRwT = zeros(NFFT_omega, NFFT_eta, M);
W2DFT_wRdwT = zeros(NFFT_omega, NFFT_eta, M);

omega_shift = zeros(NFFT_omega, NFFT_eta,M);
eta_shift = zeros(NFFT_omega, NFFT_eta,M);

omega_bins = 1:NFFT_omega;
eta_bins = 1:NFFT_eta;

for i = 1:M
    wRwT   = hT(i,:).'  .* hR(i,:);
    W2DFT_wRwT(:,:,i) = fftshift(fftshift(fft2(signal .* wRwT ,   NFFT_omega, NFFT_eta),1),2);

    dwRwT   = dhT(i,:).'  .* hR(i,:);
    W2DFT_dwRwT(:,:,i) = fftshift(fftshift(fft2(signal .* dwRwT ,   NFFT_omega, NFFT_eta),1),2);

    wRdwT   = hT(i,:).'  .* dhR(i,:);
    W2DFT_wRdwT(:,:,i) = fftshift(fftshift(fft2(signal .* wRdwT ,   NFFT_omega, NFFT_eta),1),2);

    omega_shift(:,:,i) = omega_bins   - round(imag(NFFT_omega .* W2DFT_wRdwT(:,:,i)./W2DFT_wRwT(:,:,i)./2/pi));
    eta_shift(:,:,i) = eta_bins.' - round(imag(NFFT_eta .* W2DFT_dwRwT(:,:,i)./W2DFT_wRwT(:,:,i)./2/pi));
end

W2DFT_con = zeros(size(W2DFT_wRwT));
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
                W2DFT_con(omega_idx, eta_idx, k) = W2DFT_con(omega_idx, eta_idx, k) + W2DFT_wRwT(i, j, k);
            elseif strcmp(avgtype,'noncoherent')
                W2DFT_con(omega_idx, eta_idx, k) = W2DFT_con(omega_idx, eta_idx, k) + abs(W2DFT_wRwT(i, j, k)).^2;
            end

        end
    end
end
if strcmp(avgtype,'coherent')
    W2DFT_mean =  abs(mean(W2DFT_con, 3)).^2;
elseif strcmp(avgtype,'noncoherent')
    W2DFT_mean =  mean(W2DFT_con, 3);
end

end
