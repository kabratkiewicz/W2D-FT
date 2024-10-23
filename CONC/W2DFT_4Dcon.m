function [W2DFT_original, W2DFT_con] = W2DFT_4Dcon(signal, NFFT_omega, NFFT_eta)
% inputs:
% signal - two-dimensional time-domain signal
% est - estimator number
% NFFT_omega, NFFT_eta - FFT size in the omega and eta frequency domain,
%   respectively
% output:
% W2DFT_con - concentrated 4D distribution

r = 1 : size(signal,1); r = r - mean(r)/2; Lr = length(r);
t = 1 : size(signal,2); t = t - mean(t)/2; Lt = length(t);
T = size(signal,2); R = size(signal,1);

hT   = blackman_harris_window( T, 0, 0);
dhT  = blackman_harris_window( T, 1, 0);
thT  = blackman_harris_window( T, 0, 1);
hR   = blackman_harris_window( R, 0, 0).';
dhR  = blackman_harris_window( R, 1, 0).';
thR  = blackman_harris_window( R, 0, 1).';
%%
hRhT   = hR  .* hT;
dhRhT  = dhR .* hT;
thRhT  = thR .* hT;
hRdhT  = hR  .* dhT;
hRthT  = hR  .* thT;

x = zeros(size(signal,1) + length(r),...
    size(signal,2) + length(t));

x(floor(Lr/2) + 1 : floor(Lr/2) + R, floor(Lt/2) + 1 : floor(Lt/2)+T) = signal;
t = 1 : size(x,2);  
r = 1 : size(x,1); 

omega_bins = (1 : NFFT_omega);
eta_bins   = (1 : NFFT_eta);

W2DFT_con      = zeros(NFFT_omega,NFFT_eta, T, R);
W2DFT_original = zeros(NFFT_omega,NFFT_eta,T,R);
omega_est      = zeros(NFFT_omega,NFFT_eta,T,R);
eta_est        = zeros(NFFT_omega,NFFT_eta,T,R);
t_est          = zeros(NFFT_omega,NFFT_eta,T,R);
r_est          = zeros(NFFT_omega,NFFT_eta,T,R);

for rshift = 1 : R
    for tshift = 1 : T
        x_tmp = x( rshift:rshift + R - 1, tshift:tshift + T - 1);
        W2DFT_original(:,:,tshift,rshift)  = fftshift(fftshift(fft2(x_tmp .* hRhT,  NFFT_omega, NFFT_eta),1),2);
        W2DFT_wRdwT = fftshift(fftshift(fft2(x_tmp .* hRdhT, NFFT_omega, NFFT_eta),1),2);
        W2DFT_dwRwT = fftshift(fftshift(fft2(x_tmp .* dhRhT, NFFT_omega, NFFT_eta),1),2);
        W2DFT_wRtwT = fftshift(fftshift(fft2(x_tmp .* hRthT, NFFT_omega, NFFT_eta),1),2);
        W2DFT_twRwT = fftshift(fftshift(fft2(x_tmp .* thRhT, NFFT_omega, NFFT_eta),1),2);

        omega_est(:,:,tshift,rshift) = round(omega_bins - NFFT_omega*imag(W2DFT_wRdwT ./ W2DFT_original(:,:,tshift,rshift))/2/pi );
        eta_est(:,:,tshift,rshift)   = round(eta_bins.' - NFFT_eta*imag(W2DFT_dwRwT ./ W2DFT_original(:,:,tshift,rshift))/2/pi );
        t_est(:,:,tshift,rshift)     = round(tshift   +  real(W2DFT_wRtwT ./ W2DFT_original(:,:,tshift,rshift)) );
        r_est(:,:,tshift,rshift)     = round(rshift   +  real(W2DFT_twRwT ./ W2DFT_original(:,:,tshift,rshift)) );
    end
end
%%
for tshift = 1 : T
    for rshift = 1 : R
        for i = 1 : NFFT_omega
            for j = 1 : NFFT_eta
                t_idx     = t_est(i,j, tshift,rshift);
                r_idx     = r_est(i,j, tshift,rshift);
                omega_idx = omega_est(i,j, tshift,rshift);
                eta_idx   = eta_est(i,j, tshift,rshift);
                if eta_idx < 1 || eta_idx > NFFT_eta
                    continue;
                end
                if omega_idx < 1 || omega_idx > NFFT_omega
                    continue;
                end
                if t_idx < 1 || t_idx > T
                    continue;
                end
                if r_idx < 1 || r_idx > R
                    continue;
                end

                W2DFT_con(omega_idx, eta_idx,t_idx, r_idx) = W2DFT_con(omega_idx, eta_idx,t_idx, r_idx) + abs(W2DFT_original(i,j,tshift,rshift)^2);
            end
        end
    end
end
end
