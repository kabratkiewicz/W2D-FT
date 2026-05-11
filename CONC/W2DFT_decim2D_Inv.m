% Author: Karol Abratkiewicz
% Warsaw University of Technology
% e-mail: karol.abratkiewicz@pw.edu.pl
%
% Reference:
% K. Abratkiewicz,
% "Downsampled Windowed 2D Fourier Transform in Signal Decomposition
% and Concentration,"
% IEEE Transactions on Signal Processing.

function [W2DFT_original, hRhT] = W2DFT_decim2D_Inv(signal, NFFT_omega, NFFT_eta, wT, wR, sigma_t, sigma_r, hop_t, hop_r, shift_t, shift_r, dt, dr, tt, rr)
% Computes the downsampled Windowed 2D Fourier Transform (W2DFT) of a
% two-dimensional signal using separable Gaussian analysis windows and
% periodic boundary extension.
%
% INPUTS:
% signal        - input 2D signal in the time-space domain
%
% NFFT_omega    - FFT size in the temporal frequency domain
% NFFT_eta      - FFT size in the spatial frequency domain
%
% wT            - half-length of the temporal window
% wR            - half-length of the spatial window
%
% sigma_t       - temporal spread parameter of the Gaussian window
% sigma_r       - spatial spread parameter of the Gaussian window
%
% hop_t         - temporal hop size (time-domain decimation factor)
% hop_r         - spatial hop size (space-domain decimation factor)
%
% shift_t       - temporal starting offset
% shift_r       - spatial starting offset
%
% dt            - temporal sampling step
% dr            - spatial sampling step
%
% tt            - temporal support vector for Gaussian window generation
% rr            - spatial support vector for Gaussian window generation
%
% OUTPUTS:
% W2DFT_original - complex-valued 4D W2DFT distribution
%                  of size:
%                  [NFFT_omega x NFFT_eta x T x R]
%
% hRhT           - separable 2D Gaussian analysis window
%
% NOTES:
% - The transform is evaluated only at decimated time-space positions
%   defined by hop_t and hop_r.
%
% - Circular indexing is applied near signal boundaries to ensure periodic
%   extension in both dimensions.
%
% - The implementation returns the complex W2DFT coefficients before any
%   concentration or reassignment procedure.
% =========================================================================

x = signal;

nn = -NFFT_omega/2:NFFT_omega/2-1;
mm = -NFFT_eta/2:NFFT_eta/2-1;
LT = size(signal, 1);
LR = size(signal, 2);

tau = -wT:wT;
rho = -wR:wR;
hT = Gab_Unit_Gaussian_Window(tau, sigma_t, dt, tt, 0);
hR = Gab_Unit_Gaussian_Window(rho, sigma_r, dr, rr, 0).';
hRhT           = hR  .* hT;

t = 1 + shift_t : hop_t : size(signal, 1);
r = 1 + shift_r : hop_r : size(signal, 2);
T = numel(t);
R = numel(r);
W2DFT_hRhT     = zeros(NFFT_omega, NFFT_eta, T, R);

for rshift = 1:R
    for tshift = 1:T
        
        t0 = t(tshift);
        r0 = r(rshift);
        
        if (t0 > wT) && (t0 <= LT - wT) && ...
           (r0 > wR) && (r0 <= LR - wR)
            
            W2DFT_hRhT(:,:,tshift,rshift) = ...
                exp(1j*2*pi/NFFT_omega*(wT).*nn') .* ...
                exp(1j*2*pi/NFFT_eta*(wR).*mm) .* ...
                fftshift(fftshift( ...
                    fft2(x(t0+tau, r0+rho).*hRhT, NFFT_omega, NFFT_eta),1),2);
        
        elseif ((t0 <= wT) || (t0 > LT - wT)) && ...
               ((r0 > wR) && (r0 <= LR - wR))
            
            W2DFT_hRhT(:,:,tshift,rshift) = ...
                exp(1j*2*pi/NFFT_omega*(wT).*nn') .* ...
                exp(1j*2*pi/NFFT_eta*(wR).*mm) .* ...
                fftshift(fftshift( ...
                    fft2(x(1+mod((t0-1)+tau, LT), r0+rho).*hRhT, ...
                    NFFT_omega, NFFT_eta),1),2);
        
        elseif ((t0 > wT) && (t0 <= LT - wT)) && ...
               ((r0 <= wR) || (r0 > LR - wR))
            
            W2DFT_hRhT(:,:,tshift,rshift) = ...
                exp(1j*2*pi/NFFT_omega*(wT).*nn') .* ...
                exp(1j*2*pi/NFFT_eta*(wR).*mm) .* ...
                fftshift(fftshift( ...
                    fft2(x(t0+tau, 1+mod((r0-1)+rho, LR)).*hRhT, ...
                    NFFT_omega, NFFT_eta),1),2);
        
        else
            
            W2DFT_hRhT(:,:,tshift,rshift) = ...
                exp(1j*2*pi/NFFT_omega*(wT).*nn') .* ...
                exp(1j*2*pi/NFFT_eta*(wR).*mm) .* ...
                fftshift(fftshift( ...
                    fft2(x( ...
                        1+mod((t0-1)+tau, LT), ...
                        1+mod((r0-1)+rho, LR) ...
                    ).*hRhT, ...
                    NFFT_omega, NFFT_eta),1),2);
        end
    end
end

W2DFT_original = W2DFT_hRhT;
end