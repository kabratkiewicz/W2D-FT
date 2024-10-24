% Author: Karol Abratkiewicz
% karol.abratkiewicz@pw.edu.pl  
% Warsaw University of Technology
% K. Abratkiewicz, "Windowed Two-Dimensional Fourier Transform 
% Concentration and Its Application to ISAR Imaging," in IEEE Transactions 
% on Image Processing, vol. 32, pp. 6260-6273, 2023, 
% doi: 10.1109/TIP.2023.3330603. 

function [W2DFT_original, W2DFT_con] = W2DFT_full2Dcon(signal, est, NFFT_omega, NFFT_eta)
% inputs:
% signal - two-dimensional time-domain signal
% est - estimator number
% NFFT_omega, NFFT_eta - FFT size in the omega and eta frequency domain,
%   respectively
% output:
% W2DFT_con - concentrated 4D distribution

r = 1 : size(signal,2); 
r = r - mean(r)/2; 
Lr = length(r);
t = 1 : size(signal,1); 
t = t - mean(t)/2; 
Lt = length(t);
T = size(signal, 1); 
R = size(signal, 2);

x = zeros(size(signal,1) + length(t),...
    size(signal,2) + length(r));

x(floor(Lt/2) + 1 : floor(Lt/2) + T, floor(Lr/2) + 1 : floor(Lr/2)+R) = signal;

omega_bins = (1 : NFFT_omega);
eta_bins   = (1 : NFFT_eta);

W2DFT_con      = zeros(NFFT_omega, NFFT_eta, T, R);
W2DFT_original = zeros(NFFT_omega, NFFT_eta, T, R);
omega_est      = zeros(NFFT_omega, NFFT_eta, T, R);
eta_est        = zeros(NFFT_omega, NFFT_eta, T, R);

%%%%%%%%%%%%%
switch est
    case 1
        wT             = blackman_harris_window( T, 0, 0).';
        dwT            = blackman_harris_window( T, 1, 0).';
        wR             = blackman_harris_window( R, 0, 0);
        dwR            = blackman_harris_window( R, 1, 0);

        wRwT           = wR  .* wT;
        dwRwT          = dwR .* wT;
        wRdwT          = wR  .* dwT;

        for rshift = 1 : R
            for tshift = 1 : T
                x_tmp = x( tshift:tshift + T - 1, rshift:rshift + R - 1);
                W2DFT_original(:,:,tshift,rshift)  = fftshift(fftshift(fft2(x_tmp .* wRwT,  NFFT_omega,  NFFT_eta),1),2);
                W2DFT_wRdwT                        = fftshift(fftshift(fft2(x_tmp .* wRdwT, NFFT_omega,  NFFT_eta),1),2);
                W2DFT_dwRwT                        = fftshift(fftshift(fft2(x_tmp .* dwRwT, NFFT_omega,  NFFT_eta),1),2);     
        
                omega_est(:,:,tshift,rshift) = round(omega_bins.' - NFFT_omega * imag(W2DFT_wRdwT ./ W2DFT_original(:, :, tshift, rshift))/2/pi );
                eta_est(:,:,tshift,rshift)   = round(eta_bins - NFFT_eta * imag(W2DFT_dwRwT ./ W2DFT_original(:, :, tshift, rshift))/2/pi );
            end
        end

    case 2
        wT   = blackman_harris_window( T, 0, 0).';
        dwT  = blackman_harris_window( T, 1, 0).';
        twT  = blackman_harris_window( T, 0, 1).';
        dtwT = blackman_harris_window( T, 1, 1).';
        d2wT = blackman_harris_window( T, 2, 0).';

        wR   = blackman_harris_window( R, 0, 0);
        dwR  = blackman_harris_window( R, 1, 0);
        rwR  = blackman_harris_window( R, 0, 1);
        drwR = blackman_harris_window( R, 1, 1);
        d2wR = blackman_harris_window( R, 2, 0);

        wRtwT       = wR   .* twT;
        wRdtwT      = wR   .* dtwT;
        rwRwT       = rwR  .* wT;
        drwRwT      = drwR .* wT;
        d2wRwT      = d2wR .* wT;
        wRd2wT      = wR   .* d2wT;
        wRwT        = wR   .* wT;
        dwRwT       = dwR  .* wT;
        wRdwT       = wR   .* dwT;

        for rshift = 1 : R
            for tshift = 1 : T
                x_tmp = x( tshift:tshift + T - 1, rshift:rshift + R - 1);
                W2DFT_original(:,:,tshift,rshift)  = fftshift(fftshift(fft2(x_tmp .* wRwT, NFFT_omega,  NFFT_eta),1),2);
                W2DFT_wRdwT       = fftshift(fftshift(fft2(x_tmp .* wRdwT, NFFT_omega,  NFFT_eta),1),2);
                W2DFT_dwRwT       = fftshift(fftshift(fft2(x_tmp .* dwRwT, NFFT_omega,  NFFT_eta),1),2);     
                W2DFT_wRtwT       = fftshift(fftshift(fft2(x_tmp .* wRtwT, NFFT_omega,  NFFT_eta),1),2);
                W2DFT_wRdtwT      = fftshift(fftshift(fft2(x_tmp .* wRdtwT, NFFT_omega,  NFFT_eta),1),2);
                W2DFT_wRd2wT      = fftshift(fftshift(fft2(x_tmp .* wRd2wT, NFFT_omega,  NFFT_eta),1),2);
                W2DFT_rwRwT       = fftshift(fftshift(fft2(x_tmp .* rwRwT, NFFT_omega,  NFFT_eta),1),2);
                W2DFT_drwRwT      = fftshift(fftshift(fft2(x_tmp .* drwRwT, NFFT_omega,  NFFT_eta),1),2);
                W2DFT_d2wRwT      = fftshift(fftshift(fft2(x_tmp .* d2wRwT, NFFT_omega,  NFFT_eta),1),2);

                omega_est(:,:,tshift,rshift) = round(omega_bins.' + NFFT_omega * imag((W2DFT_wRdtwT .* W2DFT_wRdwT - W2DFT_wRd2wT .* W2DFT_wRtwT)./(W2DFT_wRdwT .* W2DFT_wRtwT - W2DFT_wRdtwT .* W2DFT_original(:,:,tshift,rshift))./2/pi));
                eta_est(:,:,tshift,rshift)   = round(eta_bins + NFFT_eta   * imag((W2DFT_drwRwT .* W2DFT_dwRwT - W2DFT_d2wRwT .* W2DFT_rwRwT)./(W2DFT_dwRwT .* W2DFT_rwRwT - W2DFT_drwRwT .* W2DFT_original(:,:,tshift,rshift))./2/pi));
            end
        end
             
    case 3
        wT   = blackman_harris_window( T, 0, 0).';
        dwT  = blackman_harris_window( T, 1, 0).';
        twT  = blackman_harris_window( T, 0, 1).';
        dtwT = blackman_harris_window( T, 1, 1).';
        t2wT = blackman_harris_window( T, 0, 2).';

        wR   = blackman_harris_window( R, 0, 0);
        dwR  = blackman_harris_window( R, 1, 0);
        rwR  = blackman_harris_window( R, 0, 1);
        drwR = blackman_harris_window( R, 1, 1);
        r2wR = blackman_harris_window( R, 0, 2);

        wRtwT             = wR   .* twT;
        wRdtwT            = wR   .* dtwT;
        wRt2wT            = wR   .* t2wT;
        rwRwT             = rwR  .* wT;
        drwRwT            = drwR .* wT;
        r2wRwT            = r2wR .* wT;
        wRwT              = wR   .* wT;
        dwRwT             = dwR  .* wT;
        wRdwT             = wR   .* dwT;

        for rshift = 1 : R
            for tshift = 1 : T
                x_tmp = x( tshift:tshift + T - 1, rshift:rshift + R - 1);
                W2DFT_original(:,:,tshift,rshift)  = fftshift(fftshift(fft2(x_tmp .* wRwT, NFFT_omega,  NFFT_eta),1),2);
                W2DFT_wRdwT       = fftshift(fftshift(fft2(x_tmp .* wRdwT, NFFT_omega,  NFFT_eta),1),2);
                W2DFT_dwRwT       = fftshift(fftshift(fft2(x_tmp .* dwRwT, NFFT_omega,  NFFT_eta),1),2);
                W2DFT_wRtwT       = fftshift(fftshift(fft2(x_tmp .* wRtwT, NFFT_omega,  NFFT_eta),1),2);
                W2DFT_wRdtwT      = fftshift(fftshift(fft2(x_tmp .* wRdtwT, NFFT_omega,  NFFT_eta),1),2);
                W2DFT_wRt2wT      = fftshift(fftshift(fft2(x_tmp .* wRt2wT, NFFT_omega,  NFFT_eta),1),2);
                W2DFT_rwRwT       = fftshift(fftshift(fft2(x_tmp .* rwRwT, NFFT_omega,  NFFT_eta),1),2);
                W2DFT_drwRwT      = fftshift(fftshift(fft2(x_tmp .* drwRwT, NFFT_omega,  NFFT_eta),1),2);
                W2DFT_r2wRwT      = fftshift(fftshift(fft2(x_tmp .* r2wRwT, NFFT_omega,  NFFT_eta),1),2);

                omega_est(:,:,tshift,rshift) = round(omega_bins.' - NFFT_omega * imag((W2DFT_wRdtwT .* W2DFT_wRtwT - W2DFT_wRdwT.*W2DFT_wRt2wT + W2DFT_wRtwT .* W2DFT_original(:,:,tshift,rshift))./(W2DFT_wRtwT.^2 - W2DFT_wRt2wT.*W2DFT_original(:,:,tshift,rshift)))./2/pi);
                eta_est(:,:,tshift,rshift)   = round(eta_bins - NFFT_eta * imag((W2DFT_drwRwT .* W2DFT_rwRwT - W2DFT_dwRwT.*W2DFT_r2wRwT + W2DFT_rwRwT .* W2DFT_original(:,:,tshift,rshift))./(W2DFT_rwRwT.^2 - W2DFT_r2wRwT.*W2DFT_original(:,:,tshift,rshift)))./2/pi);

            end
        end
end
%%
for tshift = 1 : T
    for rshift = 1 : R
        for i = 1 : NFFT_omega
            for j = 1 : NFFT_eta
                omega_idx = omega_est(i,j, tshift,rshift);
                eta_idx   = eta_est(i,j, tshift,rshift);
                if eta_idx < 1 || eta_idx > NFFT_eta
                    continue;
                end
                if omega_idx < 1 || omega_idx > NFFT_omega
                    continue;
                end
                W2DFT_con(omega_idx, eta_idx,tshift, rshift) = W2DFT_con(omega_idx, eta_idx, tshift, rshift) + abs(W2DFT_original(i,j,tshift,rshift)^2);
            end
        end
    end
end
end
