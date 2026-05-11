% Author: Karol Abratkiewicz
% karol.abratkiewicz@pw.edu.pl  
% Warsaw University of Technology
% K. Abratkiewicz, "Downsampled Windowed 2D 
% Fourier Transform in Signal 
% Decomposition and Concentration," 
% in IEEE Transactions on Signal Processing

function x = IW2DFT_decim_unit(W2DFT_original, LT, LR, w, hop_t, hop_r, shift_t, shift_r)
% Inputs:
% W2DFT_original - Input downsampled W2DFT representation in the
%                  time-space-frequency-spatial frequency domain.
% LT             - Length of the reconstructed signal in the temporal dimension.
% LR             - Length of the reconstructed signal in the spatial dimension.
% w              - Two-dimensional analysis window used in the forward transform.
% hop_t          - Hop size (downsampling factor) in the temporal dimension.
% hop_r          - Hop size (downsampling factor) in the spatial dimension.
% shift_t        - Initial temporal shift of the analysis window.
% shift_r        - Initial spatial shift of the analysis window.
%
% Outputs:
% x              - Reconstructed 2D signal obtained using the inverse
%                  downsampled W2DFT.

[NFFT_omega, NFFT_eta, T, R] = size(W2DFT_original);
x = zeros(LT, LR);

[wT_len, wR_len] = size(w);
wT = floor(wT_len/2);
wR = floor(wR_len/2);

nn = (-NFFT_omega/2:NFFT_omega/2-1).';
mm = (-NFFT_eta/2:NFFT_eta/2-1);

for i = 1:LT
    
    idxT = 1+(floor(((i-1)-wT-shift_t)/hop_t):ceil(((i-1)+wT+shift_t)/hop_t));
    
    for j = 1:LR
        
        idxR = 1+(floor(((j-1)-wR-shift_r)/hop_r):ceil(((j-1)+wR+shift_r)/hop_r));
        
        acc = 0;
        cnt = 0;
        
        for mt = 1:length(idxT)
            for mr = 1:length(idxR)
                
                t0 = idxT(mt);
                r0 = idxR(mr);
                
                tau_t = (i-1) - hop_t*(t0-1) - shift_t;
                tau_r = (j-1) - hop_r*(r0-1) - shift_r;
                
                if (tau_t < -wT) || (tau_t > wT) || ...
                   (tau_r < -wR) || (tau_r > wR)
                    continue
                end
                
                if (i >= wT+1) && (i <= LT-wT) && ...
                    (j >= wR+1) && (j <= LR-wR)
                    
                    WT = W2DFT_original(:,:,t0,r0);
                
                elseif ((i <= wT) || (i > LT-wT)) && ...
                       ((j >= wR+1) && (j <= LR-wR))
                    
                    t_wrap = 1 + rem((t0-1)+T, T);
                    WT = W2DFT_original(:,:,t_wrap,r0);
                
                elseif ((i >= wT+1) && (i <= LT-wT)) && ...
                       ((j <= wR) || (j > LR-wR))
                    
                    r_wrap = 1 + rem((r0-1)+R, R);
                    WT = W2DFT_original(:,:,t0,r_wrap);
                
                else
                    
                    t_wrap = 1 + rem((t0-1)+T, T);
                    r_wrap = 1 + rem((r0-1)+R, R);
                    WT = W2DFT_original(:,:,t_wrap,r_wrap);
                end
                
                phase = exp(1j*2*pi*nn*tau_t/NFFT_omega) .* ...
                        exp(1j*2*pi*mm*tau_r/NFFT_eta);
                
                w_val = w(wT+1+tau_t, wR+1+tau_r);
                
                Y = mean(mean(WT .* phase* w_val))./ w_val.^2 ;
                
                acc = acc + Y;
                cnt = cnt + 1;
            end
        end
        
        if cnt > 0
            x(i,j) = acc / cnt;
        end
    end
end

end
