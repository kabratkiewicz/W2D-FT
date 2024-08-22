function [S_con] = W2DFT_LM(signal, NFFT_omega, NFFT_eta, mu)
N = size(signal,1);
M = size(signal,2);

wR             = blackman_harris_window( M, 0, 0);
dwR            = blackman_harris_window( M, 1, 0);
rwR            = blackman_harris_window( M, 0, 1);
drwR           = blackman_harris_window( M, 1, 1);
wT             = blackman_harris_window( N, 0, 0).';
dwT            = blackman_harris_window( N, 1, 0).';
twT            = blackman_harris_window( N, 0, 1).';
dtwT           = blackman_harris_window( N, 1, 1).';

wRwT           = wR   .* wT;
dwRwT          = dwR  .* wT;
wRdwT          = wR   .* dwT;
wRtwT          = wR   .* twT;
wRdtwT         = wR   .* dtwT;
rwRwT          = rwR  .* wT;
drwRwT         = drwR .* wT;
rwRdwT         = rwR  .* dwT;
dwRtwT         = dwR  .* twT;

F_wRwT     = fftshift(fftshift(fft2(signal .* wRwT,   NFFT_eta, NFFT_omega),1),2);
F_dwRwT    = fftshift(fftshift(fft2(signal .* dwRwT,  NFFT_eta, NFFT_omega),1),2);
F_wRdwT    = fftshift(fftshift(fft2(signal .* wRdwT,  NFFT_eta, NFFT_omega),1),2);
F_wRtwT    = fftshift(fftshift(fft2(signal .* wRtwT,  NFFT_eta, NFFT_omega),1),2);
F_wRdtwT   = fftshift(fftshift(fft2(signal .* wRdtwT, NFFT_eta, NFFT_omega),1),2);
F_rwRwT    = fftshift(fftshift(fft2(signal .* rwRwT,  NFFT_eta, NFFT_omega),1),2);
F_drwRwT   = fftshift(fftshift(fft2(signal .* drwRwT, NFFT_eta, NFFT_omega),1),2);
F_rwRdwT   = fftshift(fftshift(fft2(signal .* rwRdwT, NFFT_eta, NFFT_omega),1),2);
F_dwRtwT   = fftshift(fftshift(fft2(signal .* dwRtwT, NFFT_eta, NFFT_omega),1),2);

omega_est   = imag(F_dwRwT ./ F_wRwT);
eta_est     = imag(F_wRdwT ./ F_wRwT);

S_con = zeros(size(F_wRwT));

for i = 1 : NFFT_omega
    for j = 1 : NFFT_eta

        domega_domega = -real((F_drwRwT(j,i)*F_wRwT(j,i) - F_dwRwT(j,i)*F_rwRwT(j,i))/(F_wRwT(j,i)^2));
        domega_deta   = -real((F_dwRtwT(j,i)*F_wRwT(j,i) - F_dwRwT(j,i)*F_wRtwT(j,i))/(F_wRwT(j,i)^2));
        deta_domega   = -real((F_rwRdwT(j,i)*F_wRwT(j,i) - F_wRdwT(j,i)*F_rwRwT(j,i))/(F_wRwT(j,i)^2));
        deta_deta     = -real((F_wRdtwT(j,i)*F_wRwT(j,i) - F_wRdwT(j,i)*F_wRtwT(j,i))/(F_wRwT(j,i)^2));

        R      = [omega_est(j,i); eta_est(j,i)];
        nablaR = [mu+domega_domega, domega_deta; deta_domega  , mu + deta_deta];
        R_hat  = [nablaR(2,2) -nablaR(1,2);-nablaR(2,1) nablaR(1,1)] / (nablaR(1,1)*nablaR(2,2)-nablaR(2,1)*nablaR(1,2)) * R;
        omega_idx = i - (round(NFFT_omega/2/pi*R_hat(1)));
        eta_idx   = j - (round(NFFT_eta/2/pi*R_hat(2)));

        if eta_idx < 1 || eta_idx > NFFT_eta
            continue;
        end
        if omega_idx < 1 || omega_idx > NFFT_omega 
            continue;
        end
        S_con(eta_idx, omega_idx) = S_con(eta_idx, omega_idx) + abs(F_wRwT(j,i))^2;
    end
end
end