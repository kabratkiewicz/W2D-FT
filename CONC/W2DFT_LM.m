% Author: Karol Abratkiewicz
% karol.abratkiewicz@pw.edu.pl  
% Warsaw University of technology
% K.Abratkiewicz, M. K. Baczyk, P. Samczynski "Software-Oriented Adjustable
% Concentration of ISAR Images -- Towards Super-Resolution", IEEE Sensors
% Journal 2025

function [W2DFT_original, W2DFT_con] = W2DFT_LM(signal, NFFT_omega, NFFT_eta, mu)
% inputs:
% signal - two-dimensional time-domain signal
% NFFT_omega, NFFT_eta - FFT size in the omega and eta frequency domain,
%   respectively
% mu - dump parameter
% output:
% S_con - concentrated bivariate spectrum

R = size(signal,2);
T = size(signal,1);

hT             = blackman_harris_window( T, 0, 0).';
dhT            = blackman_harris_window( T, 1, 0).';
thT            = blackman_harris_window( T, 0, 1).';
dthT           = blackman_harris_window( T, 1, 1).';
hR             = blackman_harris_window( R, 0, 0);
dhR            = blackman_harris_window( R, 1, 0);
rhR            = blackman_harris_window( R, 0, 1);
drhR           = blackman_harris_window( R, 1, 1);

hRhT           = hR   .* hT;
dhRhT          = dhR  .* hT;
hRdhT          = hR   .* dhT;
hRthT          = hR   .* thT;
hRdthT         = hR   .* dthT;
rhRhT          = rhR  .* hT;
drhRhT         = drhR .* hT;
rhRdhT         = rhR  .* dhT;
dhRthT         = dhR  .* thT;

W2DFT_hRhT     = fftshift(fftshift(fft2(signal .* hRhT,   NFFT_omega, NFFT_eta),1),2);
W2DFT_dhRhT    = fftshift(fftshift(fft2(signal .* dhRhT,  NFFT_omega, NFFT_eta),1),2);
W2DFT_hRdhT    = fftshift(fftshift(fft2(signal .* hRdhT,  NFFT_omega, NFFT_eta),1),2);
W2DFT_hRthT    = fftshift(fftshift(fft2(signal .* hRthT,  NFFT_omega, NFFT_eta),1),2);
W2DFT_hRdthT   = fftshift(fftshift(fft2(signal .* hRdthT, NFFT_omega, NFFT_eta),1),2);
W2DFT_rhRhT    = fftshift(fftshift(fft2(signal .* rhRhT,  NFFT_omega, NFFT_eta),1),2);
W2DFT_drhRhT   = fftshift(fftshift(fft2(signal .* drhRhT, NFFT_omega, NFFT_eta),1),2);
W2DFT_rhRdhT   = fftshift(fftshift(fft2(signal .* rhRdhT, NFFT_omega, NFFT_eta),1),2);
W2DFT_dhRthT   = fftshift(fftshift(fft2(signal .* dhRthT, NFFT_omega, NFFT_eta),1),2);

eta_est    = imag(W2DFT_dhRhT ./ W2DFT_hRhT);
omega_est  = imag(W2DFT_hRdhT ./ W2DFT_hRhT);

W2DFT_con = zeros(size(W2DFT_hRhT));

for i = 1 : NFFT_omega
    for j = 1 : NFFT_eta

        deta_deta     = -real((W2DFT_drhRhT(i,j)*W2DFT_hRhT(i,j) - W2DFT_dhRhT(i,j)*W2DFT_rhRhT(i,j))/(W2DFT_hRhT(i,j)^2));
        deta_domega   = -real((W2DFT_dhRthT(i,j)*W2DFT_hRhT(i,j) - W2DFT_dhRhT(i,j)*W2DFT_hRthT(i,j))/(W2DFT_hRhT(i,j)^2));
        domega_deta   = -real((W2DFT_rhRdhT(i,j)*W2DFT_hRhT(i,j) - W2DFT_hRdhT(i,j)*W2DFT_rhRhT(i,j))/(W2DFT_hRhT(i,j)^2));
        domega_domega = -real((W2DFT_hRdthT(i,j)*W2DFT_hRhT(i,j) - W2DFT_hRdhT(i,j)*W2DFT_hRthT(i,j))/(W2DFT_hRhT(i,j)^2));

        R      = [omega_est(i,j); eta_est(i,j)];
        nablaR = [mu+domega_domega, domega_deta; deta_domega  , mu+deta_deta];
        R_hat  = [nablaR(2,2) -nablaR(1,2);-nablaR(2,1) nablaR(1,1)] / (nablaR(1,1)*nablaR(2,2)-nablaR(2,1)*nablaR(1,2)) * R;
        omega_idx = i - (round(NFFT_omega/2/pi*R_hat(1)));
        eta_idx   = j - (round(NFFT_eta/2/pi*R_hat(2)));

        if eta_idx < 1 || eta_idx > NFFT_eta
            continue;
        end
        if omega_idx < 1 || omega_idx > NFFT_omega 
            continue;
        end
        W2DFT_con(omega_idx, eta_idx) = W2DFT_con(omega_idx, eta_idx) + abs(W2DFT_hRhT(i,j))^2;
    end
end
W2DFT_original = W2DFT_hRhT;
end
