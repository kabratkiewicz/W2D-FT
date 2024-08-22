function [W2DFT_wRwT, W2DFT_mean, W2DFT_con] = BivariateMultitaper(signal, M, sigmaT, sigmaR, NFFT_omega, NFFT_eta, avgtype)

if ~exist('avgtype',"var")
    avgtype = 'noncoherent';
end

Y = size(signal,1);
X = size(signal,2);

[hT, dhT ,~] = hermf(X, M, sigmaT);
[hR, dhR, ~] = hermf(Y, M, sigmaR);

W2DFT_wRwT  = zeros(NFFT_eta, NFFT_omega, M);
W2DFT_dwRwT = zeros(NFFT_eta, NFFT_omega, M);
W2DFT_wRdwT = zeros(NFFT_eta, NFFT_omega, M);

omega_shift = zeros(NFFT_eta, NFFT_omega,M);
eta_shift = zeros(NFFT_eta, NFFT_omega,M);

omega_bins = 1:NFFT_omega;
eta_bins = 1:NFFT_eta;

for i = 1:M
    wRwT   = hT(i,:)  .* hR(i,:).';
    W2DFT_wRwT(:,:,i) = fftshift(fftshift(fft2(signal .* wRwT ,   NFFT_eta, NFFT_omega),1),2);

    dwRwT   = dhT(i,:)  .* hR(i,:).';
    W2DFT_dwRwT(:,:,i) = fftshift(fftshift(fft2(signal .* dwRwT ,   NFFT_eta, NFFT_omega),1),2);

    wRdwT   = hT(i,:)  .* dhR(i,:).';
    W2DFT_wRdwT(:,:,i) = fftshift(fftshift(fft2(signal .* wRdwT ,   NFFT_eta, NFFT_omega),1),2);

    omega_shift(:,:,i) = omega_bins   - round(imag(NFFT_omega .* W2DFT_dwRwT(:,:,i)./W2DFT_wRwT(:,:,i)./2/pi));
    eta_shift(:,:,i) = eta_bins.' - round(imag(NFFT_eta .* W2DFT_wRdwT(:,:,i)./W2DFT_wRwT(:,:,i)./2/pi));
end

W2DFT_con = zeros(size(W2DFT_wRwT));
for k = 1 : M
    for i = 1 : NFFT_omega
        for j = 1 : NFFT_eta
            t_idx = omega_shift(j, i, k);
            r_idx = eta_shift(j, i, k);
            if r_idx < 1 || r_idx > NFFT_eta
                continue;
            end
            if t_idx < 1 || t_idx > NFFT_omega 
                continue;
            end
            if strcmp(avgtype,'coherent')
                W2DFT_con(r_idx, t_idx, k) = W2DFT_con(r_idx, t_idx, k) + W2DFT_wRwT(j, i, k);
            elseif strcmp(avgtype,'noncoherent')
                W2DFT_con(r_idx, t_idx, k) = W2DFT_con(r_idx, t_idx, k) + abs(W2DFT_wRwT(j, i, k)).^2;
            end

        end
    end
end
if strcmp(avgtype,'coherent')
    W2DFT_mean =  abs(mean(W2DFT_con, 3)).^2;
elseif strcmp(avgtype,'noncoherent')
    W2DFT_mean =  mean(W2DFT_con, 3);
end

end