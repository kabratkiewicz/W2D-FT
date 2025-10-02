% Author: Karol Abratkiewicz
% karol.abratkiewicz@pw.edu.pl
% Warsaw University of Technology
% K. Abratkiewicz, "Windowed Two-Dimensional Fourier Transform
% Concentration and Its Application to ISAR Imaging," in IEEE Transactions
% on Image Processing, vol. 32, pp. 6260-6273, 2023,
% doi: 10.1109/TIP.2023.3330603.

function [W2DFT_original, W2DFT_con] = W2DFT_full2Dcon_Inv(signal, est, NFFT_omega, NFFT_eta, sigma_t, sigma_r)
% inputs:
% signal - two-dimensional time-domain signal
% est - estimator number
% NFFT_omega, NFFT_eta - FFT size in the omega and eta frequency domain,
%   respectively
% sigma_t, sigma_r - window time and space spread of the Gaussian window
%   (normalized)
% output:
% W2DFT_original - 4D distribution of the W2DFT
% W2DFT_con - concentrated W2DFT distribution

T = size(signal, 1);
R = size(signal, 2);

x = signal;

omega_bins = (1 : NFFT_omega);
eta_bins   = (1 : NFFT_eta);

W2DFT_con      = zeros(NFFT_omega, NFFT_eta, T, R);
W2DFT_hRhT     = zeros(NFFT_omega, NFFT_eta, T, R);
omega_est      = zeros(NFFT_omega, NFFT_eta, T, R);
eta_est        = zeros(NFFT_omega, NFFT_eta, T, R);

nn = -NFFT_omega/2:NFFT_omega/2-1;
mm = -NFFT_eta/2:NFFT_eta/2-1;

%%%%%%%%%%%%%
switch est
    case 1

        for rshift = 1 : R
            for tshift = 1 : T

                t_min = min(tshift-1, T);
                t_max = min(T-tshift, T);
                t = (-t_min):t_max;

                r_min = min(rshift-1, R);
                r_max = min(R-rshift, R);
                r = (-r_min):r_max;

                hT             = Gab_Gaussian_Window(t, sigma_t, 0, 0, 0);
                dhT            = Gab_Gaussian_Window(t, sigma_t, 1, 0, 0);
                hR             = Gab_Gaussian_Window(r, sigma_r, 0, 0, 0).';
                dhR            = Gab_Gaussian_Window(r, sigma_r, 1, 0, 0).';

                hRhT           = hR  .* hT;
                dhRhT          = dhR .* hT;
                hRdhT          = hR  .* dhT;

                W2DFT_hRhT(:,:,tshift,rshift)  = exp(1j*2*pi/NFFT_omega .* (tshift-1) .* nn') .* exp(1j*2*pi/NFFT_eta .* (rshift-1) .* mm) .* fftshift(fftshift(fft2(x .* hRhT,  NFFT_omega,  NFFT_eta),1),2);
                W2DFT_hRdhT                    = exp(1j*2*pi/NFFT_omega .* (tshift-1) .* nn') .* exp(1j*2*pi/NFFT_eta .* (rshift-1) .* mm) .* fftshift(fftshift(fft2(x .* hRdhT, NFFT_omega,  NFFT_eta),1),2);
                W2DFT_dhRhT                    = exp(1j*2*pi/NFFT_omega .* (tshift-1) .* nn') .* exp(1j*2*pi/NFFT_eta .* (rshift-1) .* mm) .* fftshift(fftshift(fft2(x .* dhRhT, NFFT_omega,  NFFT_eta),1),2);

                omega_est(:,:,tshift,rshift) = round(omega_bins.' - NFFT_omega * imag(W2DFT_hRdhT ./ W2DFT_hRhT(:, :, tshift, rshift))/2/pi );
                eta_est(:,:,tshift,rshift)   = round(eta_bins - NFFT_eta * imag(W2DFT_dhRhT ./ W2DFT_hRhT(:, :, tshift, rshift))/2/pi );
            end
        end

    case 2
        for rshift = 1 : R
            for tshift = 1 : T

                t_min = min(tshift-1, T);
                t_max = min(T-tshift, T);
                t = (-t_min):t_max;

                r_min = min(rshift-1, R);
                r_max = min(R-rshift, R);
                r = (-r_min):r_max;

                hT             = Gab_Gaussian_Window(t, sigma_t, 0, 0, 0);
                dhT            = Gab_Gaussian_Window(t, sigma_t, 1, 0, 0);
                thT            = Gab_Gaussian_Window(t, sigma_t, 0, 1, 0);
                dthT           = Gab_Gaussian_Window(t, sigma_t, 1, 1, 0);
                d2hT           = Gab_Gaussian_Window(t, sigma_t, 2, 0, 0);

                hR             = Gab_Gaussian_Window(r, sigma_r, 0, 0, 0).';
                dhR            = Gab_Gaussian_Window(r, sigma_r, 1, 0, 0).';
                rhR            = Gab_Gaussian_Window(r, sigma_r, 0, 1, 0).';
                drhR           = Gab_Gaussian_Window(r, sigma_r, 1, 1, 0).';
                d2hR           = Gab_Gaussian_Window(r, sigma_r, 2, 0, 0).';

                hRthT       = hR   .* thT;
                hRdthT      = hR   .* dthT;
                rhRhT       = rhR  .* hT;
                drhRhT      = drhR .* hT;
                d2hRhT      = d2hR .* hT;
                hRd2hT      = hR   .* d2hT;
                hRhT        = hR   .* hT;
                dhRhT       = dhR  .* hT;
                hRdhT       = hR   .* dhT;

                W2DFT_hRhT(:,:,tshift,rshift)  = exp(1j*2*pi/NFFT_omega .* (tshift-1) .* nn') .* exp(1j*2*pi/NFFT_eta .* (rshift-1) .* mm) .* fftshift(fftshift(fft2(x .* hRhT, NFFT_omega,  NFFT_eta),1),2);
                W2DFT_hRdhT       = exp(1j*2*pi/NFFT_omega .* (tshift-1) .* nn') .* exp(1j*2*pi/NFFT_eta .* (rshift-1) .* mm) .* fftshift(fftshift(fft2(x .* hRdhT,  NFFT_omega,  NFFT_eta),1),2);
                W2DFT_dhRhT       = exp(1j*2*pi/NFFT_omega .* (tshift-1) .* nn') .* exp(1j*2*pi/NFFT_eta .* (rshift-1) .* mm) .* fftshift(fftshift(fft2(x .* dhRhT,  NFFT_omega,  NFFT_eta),1),2);
                W2DFT_hRthT       = exp(1j*2*pi/NFFT_omega .* (tshift-1) .* nn') .* exp(1j*2*pi/NFFT_eta .* (rshift-1) .* mm) .* fftshift(fftshift(fft2(x .* hRthT,  NFFT_omega,  NFFT_eta),1),2);
                W2DFT_hRdthT      = exp(1j*2*pi/NFFT_omega .* (tshift-1) .* nn') .* exp(1j*2*pi/NFFT_eta .* (rshift-1) .* mm) .* fftshift(fftshift(fft2(x .* hRdthT, NFFT_omega,  NFFT_eta),1),2);
                W2DFT_hRd2hT      = exp(1j*2*pi/NFFT_omega .* (tshift-1) .* nn') .* exp(1j*2*pi/NFFT_eta .* (rshift-1) .* mm) .* fftshift(fftshift(fft2(x .* hRd2hT, NFFT_omega,  NFFT_eta),1),2);
                W2DFT_rhRhT       = exp(1j*2*pi/NFFT_omega .* (tshift-1) .* nn') .* exp(1j*2*pi/NFFT_eta .* (rshift-1) .* mm) .* fftshift(fftshift(fft2(x .* rhRhT,  NFFT_omega,  NFFT_eta),1),2);
                W2DFT_drhRhT      = exp(1j*2*pi/NFFT_omega .* (tshift-1) .* nn') .* exp(1j*2*pi/NFFT_eta .* (rshift-1) .* mm) .* fftshift(fftshift(fft2(x .* drhRhT, NFFT_omega,  NFFT_eta),1),2);
                W2DFT_d2hRhT      = exp(1j*2*pi/NFFT_omega .* (tshift-1) .* nn') .* exp(1j*2*pi/NFFT_eta .* (rshift-1) .* mm) .* fftshift(fftshift(fft2(x .* d2hRhT, NFFT_omega,  NFFT_eta),1),2);

                omega_est(:,:,tshift,rshift) = round(omega_bins.' + NFFT_omega * imag((W2DFT_hRdthT .* W2DFT_hRdhT - W2DFT_hRd2hT .* W2DFT_hRthT)./(W2DFT_hRdhT .* W2DFT_hRthT - W2DFT_hRdthT .* W2DFT_hRhT(:,:,tshift,rshift))./2/pi));
                eta_est(:,:,tshift,rshift)   = round(eta_bins     + NFFT_eta   * imag((W2DFT_drhRhT .* W2DFT_dhRhT - W2DFT_d2hRhT .* W2DFT_rhRhT)./(W2DFT_dhRhT .* W2DFT_rhRhT - W2DFT_drhRhT .* W2DFT_hRhT(:,:,tshift,rshift))./2/pi));
            end
        end

    case 3
        for rshift = 1 : R
            for tshift = 1 : T

                t_min = min(tshift-1, T);
                t_max = min(T-tshift, T);
                t = (-t_min):t_max;

                r_min = min(rshift-1, R);
                r_max = min(R-rshift, R);
                r = (-r_min):r_max;

                hT             = Gab_Gaussian_Window(t, sigma_t, 0, 0, 0);
                dhT            = Gab_Gaussian_Window(t, sigma_t, 1, 0, 0);
                thT            = Gab_Gaussian_Window(t, sigma_t, 0, 1, 0);
                dthT           = Gab_Gaussian_Window(t, sigma_t, 1, 1, 0);
                t2hT           = Gab_Gaussian_Window(t, sigma_t, 0, 2, 0);

                hR             = Gab_Gaussian_Window(r, sigma_r, 0, 0, 0).';
                dhR            = Gab_Gaussian_Window(r, sigma_r, 1, 0, 0).';
                rhR            = Gab_Gaussian_Window(r, sigma_r, 0, 1, 0).';
                drhR           = Gab_Gaussian_Window(r, sigma_r, 1, 1, 0).';
                r2hR           = Gab_Gaussian_Window(r, sigma_r, 0, 2, 0).';

                hRthT             = hR   .* thT;
                hRdthT            = hR   .* dthT;
                hRt2hT            = hR   .* t2hT;
                rhRhT             = rhR  .* hT;
                drhRhT            = drhR .* hT;
                r2hRhT            = r2hR .* hT;
                hRhT              = hR   .* hT;
                dhRhT             = dhR  .* hT;
                hRdhT             = hR   .* dhT;

                W2DFT_hRhT(:,:,tshift,rshift)  = exp(1j*2*pi/NFFT_omega .* (tshift-1) .* nn') .* exp(1j*2*pi/NFFT_eta .* (rshift-1) .* mm) .* fftshift(fftshift(fft2(x .* hRhT, NFFT_omega,  NFFT_eta),1),2);
                W2DFT_hRdhT       = exp(1j*2*pi/NFFT_omega .* (tshift-1) .* nn') .* exp(1j*2*pi/NFFT_eta .* (rshift-1) .* mm) .* fftshift(fftshift(fft2(x .* hRdhT,  NFFT_omega,  NFFT_eta),1),2);
                W2DFT_dhRhT       = exp(1j*2*pi/NFFT_omega .* (tshift-1) .* nn') .* exp(1j*2*pi/NFFT_eta .* (rshift-1) .* mm) .* fftshift(fftshift(fft2(x .* dhRhT,  NFFT_omega,  NFFT_eta),1),2);
                W2DFT_hRthT       = exp(1j*2*pi/NFFT_omega .* (tshift-1) .* nn') .* exp(1j*2*pi/NFFT_eta .* (rshift-1) .* mm) .* fftshift(fftshift(fft2(x .* hRthT,  NFFT_omega,  NFFT_eta),1),2);
                W2DFT_hRdthT      = exp(1j*2*pi/NFFT_omega .* (tshift-1) .* nn') .* exp(1j*2*pi/NFFT_eta .* (rshift-1) .* mm) .* fftshift(fftshift(fft2(x .* hRdthT, NFFT_omega,  NFFT_eta),1),2);
                W2DFT_hRt2hT      = exp(1j*2*pi/NFFT_omega .* (tshift-1) .* nn') .* exp(1j*2*pi/NFFT_eta .* (rshift-1) .* mm) .* fftshift(fftshift(fft2(x .* hRt2hT, NFFT_omega,  NFFT_eta),1),2);
                W2DFT_rhRhT       = exp(1j*2*pi/NFFT_omega .* (tshift-1) .* nn') .* exp(1j*2*pi/NFFT_eta .* (rshift-1) .* mm) .* fftshift(fftshift(fft2(x .* rhRhT,  NFFT_omega,  NFFT_eta),1),2);
                W2DFT_drhRhT      = exp(1j*2*pi/NFFT_omega .* (tshift-1) .* nn') .* exp(1j*2*pi/NFFT_eta .* (rshift-1) .* mm) .* fftshift(fftshift(fft2(x .* drhRhT, NFFT_omega,  NFFT_eta),1),2);
                W2DFT_r2hRhT      = exp(1j*2*pi/NFFT_omega .* (tshift-1) .* nn') .* exp(1j*2*pi/NFFT_eta .* (rshift-1) .* mm) .* fftshift(fftshift(fft2(x .* r2hRhT, NFFT_omega,  NFFT_eta),1),2);

                omega_est(:,:,tshift,rshift) = round(omega_bins.' - NFFT_omega * imag((W2DFT_hRdthT .* W2DFT_hRthT - W2DFT_hRdhT.*W2DFT_hRt2hT + W2DFT_hRthT .* W2DFT_hRhT(:,:,tshift,rshift))./(W2DFT_hRthT.^2 - W2DFT_hRt2hT.*W2DFT_hRhT(:,:,tshift,rshift)))./2/pi);
                eta_est(:,:,tshift,rshift)   = round(eta_bins - NFFT_eta * imag((W2DFT_drhRhT .* W2DFT_rhRhT - W2DFT_dhRhT.*W2DFT_r2hRhT + W2DFT_rhRhT .* W2DFT_hRhT(:,:,tshift,rshift))./(W2DFT_rhRhT.^2 - W2DFT_r2hRhT.*W2DFT_hRhT(:,:,tshift,rshift)))./2/pi);

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
                W2DFT_con(omega_idx, eta_idx,tshift, rshift) = W2DFT_con(omega_idx, eta_idx, tshift, rshift) + W2DFT_hRhT(i,j,tshift,rshift);
            end
        end
    end
end

W2DFT_original = W2DFT_hRhT;
end