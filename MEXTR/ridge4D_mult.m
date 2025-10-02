function [Cs, Es] = ridge4D_mult(C, fs, ks, nr, lambda, clwin_f, clwin_k)
% ridge4D_mult : extracts the ridges of a multicomponent 4D signal
%
% INPUTS:
%   C         : Concentrated W2DFT, size [Nf, Nk, Nt, Nr]
%   fs        : frequency axis (length Nf)
%   ks        : spatial frequency axis (length Nk)
%   nr        : number of ridges to extract
%   lambda    : regularization parameter
%   clwin_f   : clearing window in frequency direction (indices)
%   clwin_k   : clearing window in spatial-frequency direction (indices)
%
% OUTPUTS:
%   Cs        : cell array, Cs{j} = ridge trajectory of component j
%               size [Nt, Nr, 2], (:,:,1) = freq indices, (:,:,2) = spatial-freq indices
%   Es        : energies of the ridges

if nargin < 6
    clwin_f = 1;
end
if nargin < 7
    clwin_k = 1;
end

Txs = C;
[Nf, Nk, Nt, Nr] = size(C);

Cs = cell(nr,1);
Es = zeros(nr,1);

for j = 1:nr
    [Cj, Ej] = ridge4D(Txs, fs, ks, lambda);
    Cs{j} = Cj;
    Es(j) = Ej;

    for t = 1:Nt
        for r = 1:Nr
            f_idx = Cj(t,r,1);
            k_idx = Cj(t,r,2);

            fmin = max(1, f_idx - clwin_f);
            fmax = min(Nf, f_idx + clwin_f);
            kmin = max(1, k_idx - clwin_k);
            kmax = min(Nk, k_idx + clwin_k);

            Txs(fmin:fmax, kmin:kmax, t, r) = 0;
        end
    end
end

end
