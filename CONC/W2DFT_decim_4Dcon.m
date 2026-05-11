% Author: Karol Abratkiewicz
% karol.abratkiewicz@pw.edu.pl  
% Warsaw University of Technology
% K. Abratkiewicz, "Downsampled Windowed 2D 
% Fourier Transform in Signal 
% Decomposition and Concentration," 
% in IEEE Transactions on Signal Processing


function [W2DFT_original, W2DFT_con] = W2DFT_decim_4Dcon(signal, NFFT_omega, NFFT_eta, wT, wR, sigma_t, sigma_r, hop_t, hop_r, shift_t, shift_r)
% Inputs:
% signal       - Input 2D signal or image matrix.
% NFFT_omega   - Number of frequency bins in the temporal frequency direction.
% NFFT_eta     - Number of frequency bins in the spatial frequency direction.
% wT           - Length of the analysis window in the temporal dimension.
% wR           - Length of the analysis window in the spatial dimension.
% sigma_t      - Standard deviation of the Gaussian window in the temporal dimension.
% sigma_r      - Standard deviation of the Gaussian window in the spatial dimension.
% hop_t        - Hop size (downsampling factor) in the temporal dimension.
% hop_r        - Hop size (downsampling factor) in the spatial dimension.
% shift_t      - Initial temporal shift of the analysis window.
% shift_r      - Initial spatial shift of the analysis window.
%
% Outputs:
% W2DFT_original - Original downsampled W2DFT representation.
% W2DFT_con      - Four-dimensionally concentrated W2DFT obtained
%                  using reassignment in the time-space-frequency-
%                  spatial frequency domain.

[W2DFT_h, ~]     = W2DFT_decim2D_Inv(signal, NFFT_omega, NFFT_eta, wT, wR, sigma_t, sigma_r, hop_t, hop_r, shift_t, shift_r, 0, 0, 0, 0);
[W2DFT_dth, ~] = W2DFT_decim2D_Inv(signal, NFFT_omega, NFFT_eta, wT, wR, sigma_t, sigma_r, hop_t, hop_r, shift_t, shift_r, 1, 0, 0, 0);
[W2DFT_drh, ~] = W2DFT_decim2D_Inv(signal, NFFT_omega, NFFT_eta, wT, wR, sigma_t, sigma_r, hop_t, hop_r, shift_t, shift_r, 0, 1, 0, 0);
[W2DFT_th, ~]   = W2DFT_decim2D_Inv(signal, NFFT_omega, NFFT_eta, wT, wR, sigma_t, sigma_r, hop_t, hop_r, shift_t, shift_r, 0, 0, 1, 0);
[W2DFT_rh, ~]   = W2DFT_decim2D_Inv(signal, NFFT_omega, NFFT_eta, wT, wR, sigma_t, sigma_r, hop_t, hop_r, shift_t, shift_r, 0, 0, 0, 1);

omega_bins = (1 : NFFT_omega);
eta_bins   = (1 : NFFT_eta);

T_tmp = 1 + shift_t : hop_t : size(signal, 1);
R_tmp = 1 + shift_r : hop_r : size(signal, 2);


omega_est      = zeros(size(W2DFT_h));
eta_est        = zeros(size(W2DFT_h));
t_est          = zeros(size(W2DFT_h));
r_est          = zeros(size(W2DFT_h));

for rshift = 1 : length(R_tmp)
    for tshift = 1 : length(T_tmp)
        omega_est(:,:,tshift,rshift) = round(omega_bins.' - NFFT_omega*imag(W2DFT_dth(:,:,tshift,rshift) ./ W2DFT_h(:,:,tshift,rshift))/2/pi );
        eta_est(:,:,tshift,rshift)   = round(eta_bins - NFFT_eta*imag(W2DFT_drh(:,:,tshift,rshift) ./ W2DFT_h(:,:,tshift,rshift))/2/pi );
        t_est(:,:,tshift,rshift)     = round((T_tmp(tshift)   +  real(W2DFT_th(:,:,tshift,rshift) ./ W2DFT_h(:,:,tshift,rshift)) )./hop_t);
        r_est(:,:,tshift,rshift)     = round((R_tmp(rshift)   +  real(W2DFT_rh(:,:,tshift,rshift) ./ W2DFT_h(:,:,tshift,rshift)) )./hop_r );
    end
end

W2DFT_con      = zeros(size(W2DFT_h));

for tshift = 1 : length(T_tmp)
    for rshift = 1 : length(R_tmp)
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
                if t_idx < 1 || t_idx > numel(T_tmp)
                    continue;
                end
                if r_idx < 1 || r_idx > numel(R_tmp)
                    continue;
                end

                W2DFT_con(omega_idx, eta_idx,t_idx, r_idx) = W2DFT_con(omega_idx, eta_idx,t_idx, r_idx) + abs(W2DFT_h(i,j,tshift,rshift))^2;
            end
        end
    end
end
W2DFT_original = W2DFT_h;

end