% Author: Karol Abratkiewicz
% karol.abratkiewicz@pw.edu.pl  
% Warsaw University of technology
% K. Abratkiewicz, "Windowed Two-Dimensional Fourier Transform 
% Concentration and Its Application to ISAR Imaging," in IEEE Transactions 
% on Image Processing, vol. 32, pp. 6260-6273, 2023, 
% doi: 10.1109/TIP.2023.3330603. 

function [W2DFT_wRwT, W2DFT_con] = W2DFT(signal, est, NFFT_omega, NFFT_eta)
% inputs:
% signal - two-dimensional time-domain signal
% est - estimator number
% NFFT_omega, NFFT_eta - FFT size in the omega and eta frequency domain,
%   respectively
% output:
% W2DFT_wRwT - original windowed bivariate spectrum
% W2DFT_con - concentrated bivariate spectrum


R = size(signal,1);
T = size(signal,2);
omega_bins = (1 : NFFT_omega);
eta_bins = (1 : NFFT_eta);

switch est
    case 1
        wT             = blackman_harris_window( T, 0, 0);
        dwT            = blackman_harris_window( T, 1, 0);
        wR             = blackman_harris_window( R, 0, 0).';
        dwR            = blackman_harris_window( R, 1, 0).';
        wRwT           = wR  .* wT;
        dwRwT          = dwR .* wT;
        wRdwT          = wR  .* dwT;
        W2DFT_wRwT     = fftshift(fftshift(fft2(signal .* wRwT,   NFFT_eta, NFFT_omega),1),2);
        W2DFT_wRdwT    = fftshift(fftshift(fft2(signal .* wRdwT,  NFFT_eta, NFFT_omega),1),2);
        W2DFT_dwRwT    = fftshift(fftshift(fft2(signal .* dwRwT,  NFFT_eta, NFFT_omega),1),2);

        omega_shift = round(omega_bins -  NFFT_omega*imag(W2DFT_wRdwT ./ W2DFT_wRwT)/2/pi );
        eta_shift   = round(eta_bins.' -  NFFT_eta*imag(W2DFT_dwRwT ./ W2DFT_wRwT)/2/pi );
    case 2
        wT   = blackman_harris_window( T, 0, 0);
        dwT  = blackman_harris_window( T, 1, 0);
        twT  = blackman_harris_window( T, 0, 1);
        dtwT = blackman_harris_window( T, 1, 1);
        d2wT = blackman_harris_window( T, 2, 0);

        wR   = blackman_harris_window( R, 0, 0).';
        dwR  = blackman_harris_window( R, 1, 0).';
        rwR  = blackman_harris_window( R, 0, 1).';
        drwR = blackman_harris_window( R, 1, 1).';
        d2wR = blackman_harris_window( R, 2, 0).';

        wRtwT       = wR   .* twT;
        wRdtwT      = wR   .* dtwT;
        rwRwT       = rwR  .* wT;
        drwRwT      = drwR .* wT;
        d2wRwT      = d2wR .* wT;
        wRd2wT      = wR   .* d2wT;
        wRwT        = wR   .* wT;
        dwRwT       = dwR  .* wT;
        wRdwT       = wR   .* dwT;
        
        W2DFT_wRwT        = fftshift(fftshift(fft2(signal .* wRwT,   NFFT_eta, NFFT_omega),1),2);
        W2DFT_wRdwT       = fftshift(fftshift(fft2(signal .* wRdwT,  NFFT_eta, NFFT_omega),1),2);
        W2DFT_dwRwT       = fftshift(fftshift(fft2(signal .* dwRwT,  NFFT_eta, NFFT_omega),1),2);
        W2DFT_wRtwT       = fftshift(fftshift(fft2(signal .* wRtwT,  NFFT_eta, NFFT_omega),1),2);
        W2DFT_wRdtwT      = fftshift(fftshift(fft2(signal .* wRdtwT, NFFT_eta, NFFT_omega),1),2);
        W2DFT_wRd2wT      = fftshift(fftshift(fft2(signal .* wRd2wT, NFFT_eta, NFFT_omega),1),2);
        W2DFT_rwRwT       = fftshift(fftshift(fft2(signal .* rwRwT,  NFFT_eta, NFFT_omega),1),2);
        W2DFT_drwRwT      = fftshift(fftshift(fft2(signal .* drwRwT, NFFT_eta, NFFT_omega),1),2);
        W2DFT_d2wRwT      = fftshift(fftshift(fft2(signal .* d2wRwT, NFFT_eta, NFFT_omega),1),2);

        omega_shift = round(omega_bins + NFFT_omega * imag((W2DFT_wRdtwT .* W2DFT_wRdwT - W2DFT_wRd2wT .* W2DFT_wRtwT)./(W2DFT_wRdwT .* W2DFT_wRtwT - W2DFT_wRdtwT .* W2DFT_wRwT)./2/pi));
        eta_shift = round(eta_bins.'   + NFFT_eta * imag((W2DFT_drwRwT .* W2DFT_dwRwT - W2DFT_d2wRwT .* W2DFT_rwRwT)./(W2DFT_dwRwT .* W2DFT_rwRwT - W2DFT_drwRwT .* W2DFT_wRwT)./2/pi));
        
    case 3
        wT   = blackman_harris_window( T, 0, 0);
        dwT  = blackman_harris_window( T, 1, 0);
        twT  = blackman_harris_window( T, 0, 1);
        dtwT = blackman_harris_window( T, 1, 1);
        t2wT = blackman_harris_window( T, 0, 2);

        wR   = blackman_harris_window( R, 0, 0).';
        dwR  = blackman_harris_window( R, 1, 0).';
        rwR  = blackman_harris_window( R, 0, 1).';
        drwR = blackman_harris_window( R, 1, 1).';
        r2wR = blackman_harris_window( R, 0, 2).';

        wRtwT             = wR   .* twT;
        wRdtwT            = wR   .* dtwT;
        wRt2wT            = wR   .* t2wT;
        rwRwT             = rwR  .* wT;
        drwRwT            = drwR .* wT;
        r2wRwT            = r2wR .* wT;
        wRwT              = wR   .* wT;
        dwRwT             = dwR  .* wT;
        wRdwT             = wR   .* dwT;

        W2DFT_wRwT        = fftshift(fftshift(fft2(signal .* wRwT,   NFFT_eta, NFFT_omega),1),2);
        W2DFT_wRdwT       = fftshift(fftshift(fft2(signal .* wRdwT,  NFFT_eta, NFFT_omega),1),2);
        W2DFT_dwRwT       = fftshift(fftshift(fft2(signal .* dwRwT,  NFFT_eta, NFFT_omega),1),2);
        W2DFT_wRtwT       = fftshift(fftshift(fft2(signal .* wRtwT,  NFFT_eta, NFFT_omega),1),2);
        W2DFT_wRdtwT      = fftshift(fftshift(fft2(signal .* wRdtwT, NFFT_eta, NFFT_omega),1),2);
        W2DFT_wRt2wT      = fftshift(fftshift(fft2(signal .* wRt2wT, NFFT_eta, NFFT_omega),1),2);
        W2DFT_rwRwT       = fftshift(fftshift(fft2(signal .* rwRwT,  NFFT_eta, NFFT_omega),1),2);
        W2DFT_drwRwT      = fftshift(fftshift(fft2(signal .* drwRwT, NFFT_eta, NFFT_omega),1),2);
        W2DFT_r2wRwT      = fftshift(fftshift(fft2(signal .* r2wRwT, NFFT_eta, NFFT_omega),1),2);

        omega_shift   = round(omega_bins - NFFT_omega * imag((W2DFT_wRdtwT .* W2DFT_wRtwT - W2DFT_wRdwT.*W2DFT_wRt2wT + W2DFT_wRtwT .* W2DFT_wRwT)./(W2DFT_wRtwT.^2 - W2DFT_wRt2wT.*W2DFT_wRwT))./2/pi);
        eta_shift = round(eta_bins.' - NFFT_eta * imag((W2DFT_drwRwT .* W2DFT_rwRwT - W2DFT_dwRwT.*W2DFT_r2wRwT + W2DFT_rwRwT .* W2DFT_wRwT)./(W2DFT_rwRwT.^2 - W2DFT_r2wRwT.*W2DFT_wRwT))./2/pi);
end

W2DFT_con = zeros(size(W2DFT_wRwT));
for i = 1 : NFFT_omega
    for j = 1 : NFFT_eta
        omega_idx = omega_shift(j,i);
        eta_idx = eta_shift(j,i);
        if eta_idx < 1 || eta_idx > NFFT_eta
            continue;
        end
        if omega_idx < 1 || omega_idx > NFFT_omega
            continue;
        end
        W2DFT_con(eta_idx, omega_idx) = W2DFT_con(eta_idx, omega_idx) + abs(W2DFT_wRwT(j,i))^2;
    end
end

end
