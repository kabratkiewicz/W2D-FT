% Author: Karol Abratkiewicz
% karol.abratkiewicz@pw.edu.pl  
% Warsaw University of Technology
% K.Abratkiewicz, M. K. Baczyk, P. Samczynski "Software-Oriented Adjustable
% Concentration of ISAR Images -- Towards Super-Resolution", IEEE Sensors
% Journal 2025

function [S_con] = W2DFT_LM(signal, NFFT_omega, NFFT_eta, mu)
% inputs:
% signal - two-dimensional time-domain signal
% NFFT_omega, NFFT_eta - FFT size in the omega and eta frequency domain,
%   respectively
% mu - dump parameter
% output:
% S_con - concentrated bivariate spectrum

R = size(signal,1);
T = size(signal,2);

wT             = blackman_harris_window( T, 0, 0).';
dwT            = blackman_harris_window( T, 1, 0).';
twT            = blackman_harris_window( T, 0, 1).';
dtwT           = blackman_harris_window( T, 1, 1).';
wR             = blackman_harris_window( R, 0, 0);
dwR            = blackman_harris_window( R, 1, 0);
rwR            = blackman_harris_window( R, 0, 1);
drwR           = blackman_harris_window( R, 1, 1);

wRwT           = wR   .* wT;
dwRwT          = dwR  .* wT;
wRdwT          = wR   .* dwT;
wRtwT          = wR   .* twT;
wRdtwT         = wR   .* dtwT;
rwRwT          = rwR  .* wT;
drwRwT         = drwR .* wT;
rwRdwT         = rwR  .* dwT;
dwRtwT         = dwR  .* twT;

F_wRwT     = fftshift(fftshift(fft2(signal .* wRwT,   NFFT_omega, NFFT_eta),1),2);
F_dwRwT    = fftshift(fftshift(fft2(signal .* dwRwT,  NFFT_omega, NFFT_eta),1),2);
F_wRdwT    = fftshift(fftshift(fft2(signal .* wRdwT,  NFFT_omega, NFFT_eta),1),2);
F_wRtwT    = fftshift(fftshift(fft2(signal .* wRtwT,  NFFT_omega, NFFT_eta),1),2);
F_wRdtwT   = fftshift(fftshift(fft2(signal .* wRdtwT, NFFT_omega, NFFT_eta),1),2);
F_rwRwT    = fftshift(fftshift(fft2(signal .* rwRwT,  NFFT_omega, NFFT_eta),1),2);
F_drwRwT   = fftshift(fftshift(fft2(signal .* drwRwT, NFFT_omega, NFFT_eta),1),2);
F_rwRdwT   = fftshift(fftshift(fft2(signal .* rwRdwT, NFFT_omega, NFFT_eta),1),2);
F_dwRtwT   = fftshift(fftshift(fft2(signal .* dwRtwT, NFFT_omega, NFFT_eta),1),2);

eta_est    = imag(F_dwRwT ./ F_wRwT);
omega_est  = imag(F_wRdwT ./ F_wRwT);

S_con = zeros(size(F_wRwT));

for i = 1 : NFFT_omega
    for j = 1 : NFFT_eta

        deta_deta = -real((F_drwRwT(i,j)*F_wRwT(i,j) - F_dwRwT(i,j)*F_rwRwT(i,j))/(F_wRwT(i,j)^2));
        deta_domega   = -real((F_dwRtwT(i,j)*F_wRwT(i,j) - F_dwRwT(i,j)*F_wRtwT(i,j))/(F_wRwT(i,j)^2));
        domega_deta   = -real((F_rwRdwT(i,j)*F_wRwT(i,j) - F_wRdwT(i,j)*F_rwRwT(i,j))/(F_wRwT(i,j)^2));
        domega_domega     = -real((F_wRdtwT(i,j)*F_wRwT(i,j) - F_wRdwT(i,j)*F_wRtwT(i,j))/(F_wRwT(i,j)^2));

        R      = [omega_est(i,j); eta_est(i,j)];
        nablaR = [mu+domega_domega, +domega_deta; +deta_domega  , mu+deta_deta];
        R_hat  = [nablaR(2,2) -nablaR(1,2);-nablaR(2,1) nablaR(1,1)] / (nablaR(1,1)*nablaR(2,2)-nablaR(2,1)*nablaR(1,2)) * R;
        omega_idx = i - (round(NFFT_omega/2/pi*R_hat(1)));
        eta_idx   = j - (round(NFFT_eta/2/pi*R_hat(2)));

        if eta_idx < 1 || eta_idx > NFFT_eta
            continue;
        end
        if omega_idx < 1 || omega_idx > NFFT_omega 
            continue;
        end
        S_con(omega_idx, eta_idx) = S_con(omega_idx, eta_idx) + abs(F_wRwT(i,j))^2;
    end
end
end
