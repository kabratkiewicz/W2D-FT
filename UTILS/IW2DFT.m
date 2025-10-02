function rec = IW2DFT(W2DFT_distribution, sigma_t, sigma_r)

NFFT_omega = size(W2DFT_distribution, 1);
NFFT_eta   = size(W2DFT_distribution, 2);
T          = size(W2DFT_distribution, 3);
R          = size(W2DFT_distribution, 4);

rec = zeros(T, R);
for i = 1:T
    for j = 1:R
        rec(i,j) = 2 .* pi .* sum(sum(W2DFT_distribution(:,:,i,j))) ./ NFFT_eta ./ NFFT_omega .* sigma_t .* sigma_r;
    end
end

end