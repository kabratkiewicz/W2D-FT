function [Y] = W2DFT_comp_ext4D(C, mask, sigma_t, sigma_r)
% W2DFT_comp_ext4D : reconstructs components from the concentrated 4D W2D-FT representation
%
% INPUTS:
%   C      : 4D SST, size [Nf, Nk, Nt, Nr]
%   mask   : 5D logical mask, size [Nf, Nk, Nt, Nr, nb_comp]
%   sigma_t, sigma_r      : window time/range spread
%
% OUTPUT:
%   Y      : reconstructed components, size [nb_comp, Nt, Nr]



[Nf, Nk, Nt, Nr] = size(C);
n_comp = size(mask, 5);

Y = zeros(n_comp, Nt, Nr);

for i = 1:n_comp
    % apply mask
    tmp_stfr = zeros(Nf, Nk, Nt, Nr);
    tmp_stfr(mask(:,:,:,:,i)) = C(mask(:,:,:,:,i));

    Y(i,:,:) = IW2DFT(tmp_stfr, sigma_t, sigma_r);
end

end