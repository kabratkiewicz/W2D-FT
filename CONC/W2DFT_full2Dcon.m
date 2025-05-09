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
r = r - mean(r); 
Lr = length(r);
t = 1 : size(signal,1); 
t = t - mean(t); 
Lt = length(t);
T = size(signal, 1); 
R = size(signal, 2);

x = zeros(size(signal,1) + length(t),...
    size(signal,2) + length(r));

x(floor(Lt/2) + 1 : floor(Lt/2) + T, floor(Lr/2) + 1 : floor(Lr/2)+R) = signal;

omega_bins = (1 : NFFT_omega);
eta_bins   = (1 : NFFT_eta);

W2DFT_con      = zeros(NFFT_omega, NFFT_eta, T, R);
W2DFT_hRhT = zeros(NFFT_omega, NFFT_eta, T, R);
omega_est      = zeros(NFFT_omega, NFFT_eta, T, R);
eta_est        = zeros(NFFT_omega, NFFT_eta, T, R);

%%%%%%%%%%%%%
switch est
    case 1
        hT             = blackman_harris_window( T, 0, 0).';
        dhT            = blackman_harris_window( T, 1, 0).';
        hR             = blackman_harris_window( R, 0, 0);
        dhR            = blackman_harris_window( R, 1, 0);

        hRhT           = hR  .* hT;
        dhRhT          = dhR .* hT;
        hRdhT          = hR  .* dhT;

        for rshift = 1 : R
            for tshift = 1 : T
                x_tmp = x( tshift:tshift + T - 1, rshift:rshift + R - 1);
                W2DFT_hRhT(:,:,tshift,rshift)  = fftshift(fftshift(fft2(x_tmp .* hRhT,  NFFT_omega,  NFFT_eta),1),2);
                W2DFT_hRdhT                        = fftshift(fftshift(fft2(x_tmp .* hRdhT, NFFT_omega,  NFFT_eta),1),2);
                W2DFT_dhRhT                        = fftshift(fftshift(fft2(x_tmp .* dhRhT, NFFT_omega,  NFFT_eta),1),2);     
        
                omega_est(:,:,tshift,rshift) = round(omega_bins.' - NFFT_omega * imag(W2DFT_hRdhT ./ W2DFT_hRhT(:, :, tshift, rshift))/2/pi );
                eta_est(:,:,tshift,rshift)   = round(eta_bins - NFFT_eta * imag(W2DFT_dhRhT ./ W2DFT_hRhT(:, :, tshift, rshift))/2/pi );
            end
        end

    case 2
        hT   = blackman_harris_window( T, 0, 0).';
        dhT  = blackman_harris_window( T, 1, 0).';
        thT  = blackman_harris_window( T, 0, 1).';
        dthT = blackman_harris_window( T, 1, 1).';
        d2hT = blackman_harris_window( T, 2, 0).';

        hR   = blackman_harris_window( R, 0, 0);
        dhR  = blackman_harris_window( R, 1, 0);
        rhR  = blackman_harris_window( R, 0, 1);
        drhR = blackman_harris_window( R, 1, 1);
        d2hR = blackman_harris_window( R, 2, 0);

        hRthT       = hR   .* thT;
        hRdthT      = hR   .* dthT;
        rhRhT       = rhR  .* hT;
        drhRhT      = drhR .* hT;
        d2hRhT      = d2hR .* hT;
        hRd2hT      = hR   .* d2hT;
        hRhT        = hR   .* hT;
        dhRhT       = dhR  .* hT;
        hRdhT       = hR   .* dhT;

        for rshift = 1 : R
            for tshift = 1 : T
                x_tmp = x( tshift:tshift + T - 1, rshift:rshift + R - 1);
                W2DFT_hRhT(:,:,tshift,rshift)  = fftshift(fftshift(fft2(x_tmp .* hRhT, NFFT_omega,  NFFT_eta),1),2);
                W2DFT_hRdhT       = fftshift(fftshift(fft2(x_tmp .* hRdhT, NFFT_omega,  NFFT_eta),1),2);
                W2DFT_dhRhT       = fftshift(fftshift(fft2(x_tmp .* dhRhT, NFFT_omega,  NFFT_eta),1),2);     
                W2DFT_hRthT       = fftshift(fftshift(fft2(x_tmp .* hRthT, NFFT_omega,  NFFT_eta),1),2);
                W2DFT_hRdthT      = fftshift(fftshift(fft2(x_tmp .* hRdthT, NFFT_omega,  NFFT_eta),1),2);
                W2DFT_hRd2hT      = fftshift(fftshift(fft2(x_tmp .* hRd2hT, NFFT_omega,  NFFT_eta),1),2);
                W2DFT_rhRhT       = fftshift(fftshift(fft2(x_tmp .* rhRhT, NFFT_omega,  NFFT_eta),1),2);
                W2DFT_drhRhT      = fftshift(fftshift(fft2(x_tmp .* drhRhT, NFFT_omega,  NFFT_eta),1),2);
                W2DFT_d2hRhT      = fftshift(fftshift(fft2(x_tmp .* d2hRhT, NFFT_omega,  NFFT_eta),1),2);

                omega_est(:,:,tshift,rshift) = round(omega_bins.' + NFFT_omega * imag((W2DFT_hRdthT .* W2DFT_hRdhT - W2DFT_hRd2hT .* W2DFT_hRthT)./(W2DFT_hRdhT .* W2DFT_hRthT - W2DFT_hRdthT .* W2DFT_hRhT(:,:,tshift,rshift))./2/pi));
                eta_est(:,:,tshift,rshift)   = round(eta_bins + NFFT_eta   * imag((W2DFT_drhRhT .* W2DFT_dhRhT - W2DFT_d2hRhT .* W2DFT_rhRhT)./(W2DFT_dhRhT .* W2DFT_rhRhT - W2DFT_drhRhT .* W2DFT_hRhT(:,:,tshift,rshift))./2/pi));

            end
        end
        
        
    case 3
        hT   = blackman_harris_window( T, 0, 0).';
        dhT  = blackman_harris_window( T, 1, 0).';
        thT  = blackman_harris_window( T, 0, 1).';
        dthT = blackman_harris_window( T, 1, 1).';
        t2hT = blackman_harris_window( T, 0, 2).';

        hR   = blackman_harris_window( R, 0, 0);
        dhR  = blackman_harris_window( R, 1, 0);
        rhR  = blackman_harris_window( R, 0, 1);
        drhR = blackman_harris_window( R, 1, 1);
        r2hR = blackman_harris_window( R, 0, 2);

        hRthT             = hR   .* thT;
        hRdthT            = hR   .* dthT;
        hRt2hT            = hR   .* t2hT;
        rhRhT             = rhR  .* hT;
        drhRhT            = drhR .* hT;
        r2hRhT            = r2hR .* hT;
        hRhT              = hR   .* hT;
        dhRhT             = dhR  .* hT;
        hRdhT             = hR   .* dhT;

        for rshift = 1 : R
            for tshift = 1 : T
                x_tmp = x( tshift:tshift + T - 1, rshift:rshift + R - 1);
                W2DFT_hRhT(:,:,tshift,rshift)  = fftshift(fftshift(fft2(x_tmp .* hRhT, NFFT_omega,  NFFT_eta),1),2);
                W2DFT_hRdhT       = fftshift(fftshift(fft2(x_tmp .* hRdhT, NFFT_omega,  NFFT_eta),1),2);
                W2DFT_dhRhT       = fftshift(fftshift(fft2(x_tmp .* dhRhT, NFFT_omega,  NFFT_eta),1),2);
                W2DFT_hRthT       = fftshift(fftshift(fft2(x_tmp .* hRthT, NFFT_omega,  NFFT_eta),1),2);
                W2DFT_hRdthT      = fftshift(fftshift(fft2(x_tmp .* hRdthT, NFFT_omega,  NFFT_eta),1),2);
                W2DFT_hRt2hT      = fftshift(fftshift(fft2(x_tmp .* hRt2hT, NFFT_omega,  NFFT_eta),1),2);
                W2DFT_rhRhT       = fftshift(fftshift(fft2(x_tmp .* rhRhT, NFFT_omega,  NFFT_eta),1),2);
                W2DFT_drhRhT      = fftshift(fftshift(fft2(x_tmp .* drhRhT, NFFT_omega,  NFFT_eta),1),2);
                W2DFT_r2hRhT      = fftshift(fftshift(fft2(x_tmp .* r2hRhT, NFFT_omega,  NFFT_eta),1),2);

                omega_est(:,:,tshift,rshift) = round(omega_bins.' - NFFT_omega * imag((W2DFT_hRdthT .* W2DFT_hRthT - W2DFT_hRdhT.*W2DFT_hRt2hT + W2DFT_hRthT .* W2DFT_hRhT(:,:,tshift,rshift))./(W2DFT_hRthT.^2 - W2DFT_hRt2hT.*W2DFT_hRhT(:,:,tshift,rshift)))./2/pi);
                eta_est(:,:,tshift,rshift)   = round(eta_bins - NFFT_eta * imag((W2DFT_drhRhT .* W2DFT_rhRhT - W2DFT_dhRhT.*W2DFT_r2hRhT + W2DFT_rhRhT .* W2DFT_hRhT(:,:,tshift,rshift))./(W2DFT_rhRhT.^2 - W2DFT_r2hRhT.*W2DFT_hRhT(:,:,tshift,rshift)))./2/pi);

            end
        end
end



%%%%%%%%%%%%%%


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
                W2DFT_con(omega_idx, eta_idx,tshift, rshift) = W2DFT_con(omega_idx, eta_idx, tshift, rshift) + abs(W2DFT_hRhT(i,j,tshift,rshift))^2;
            end
        end
    end
end

W2DFT_original = W2DFT_hRhT;
end
