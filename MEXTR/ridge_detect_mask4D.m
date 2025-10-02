function [mask, Cs] = ridge_detect_mask4D(C, fs, ks, nr, lambda, clwin_f, clwin_k, Kf, Kk)
% ridge_detect_mask4D : builds masks around detected ridges in 4D
%
% INPUTS:
%   C         : concentrated W2DFT transform [Nf, Nk, Nt, Nr]
%   fs        : frequency axis
%   ks        : spatial frequency axis
%   nr        : number of ridges to detect
%   lambda    : regularization parameter
%   clwin_f   : clearing window (frequency)
%   clwin_k   : clearing window (spatial-frequency)
%   Kf        : mask half-width in frequency
%   Kk        : mask half-width in spatial frequency
%
% OUTPUTS:
%   mask      : 5D logical mask [Nf, Nk, Nt, Nr, nr]
%   Cs        : cell array of ridge trajectories

[Nf, Nk, Nt, Nr] = size(C);

[Cs, ~] = ridge4D_mult(C, fs, ks, nr, lambda, clwin_f, clwin_k);

mask = false(Nf, Nk, Nt, Nr, nr);

for j = 1:nr
    Cj = Cs{j};
    for t = 1:Nt
        for r = 1:Nr
            f_idx = Cj(t,r,1);
            k_idx = Cj(t,r,2);

            if f_idx > 1 && k_idx > 1 && abs(C(f_idx,k_idx,t,r)) > lambda
                
                f_range = max(1, f_idx - Kf) : min(Nf, f_idx + Kf);
                k_range = max(1, k_idx - Kk) : min(Nk, k_idx + Kk);

                mask(f_range, k_range, t, r, j) = true;
            end
        end
    end
end

end
