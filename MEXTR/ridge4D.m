function [c, e] = ridge4D(C, fs, ks, lambda)
% ridge4D : greedy extraction of a 4D ridge (f, k, t, r)
% INPUTS:
%   C      - 4D transform [Nf x Nk x Nt x Nr]
%   fs     - frequency vector (length Nf)
%   ks     - spatial frequency vector (length Nk)
%   lambda - regularization weight
%
% OUTPUTS:
%   c - [Nt x Nr x 2], indices of the ridge:
%         c(:,:,1) = index in frequency dimension (f)
%         c(:,:,2) = index in spatial frequency dimension (k)
%   e - scalar energy of the ridge

Et = log(abs(C) + eps^0.25);

domega = fs(2) - fs(1);
dk     = ks(2) - ks(1);

[Nf, Nk, Nt, Nr] = size(C);

% --- parameters ---
da_f = 10;     % maximal jump allowed in f
da_k = 10;     % maximal jump allowed in k
ng   = 100;    % number of random initializations
aux_f = 0.5 * lambda * (domega^2);
aux_k = 0.5 * lambda * (dk^2);

c = zeros(Nt, Nr, 2);
e = -Inf;

% seed times for random initialization
t_seeds = unique(round(linspace( max(1,ceil(Nt/(ng+1))), min(Nt, Nt-ceil(Nt/(ng+1))), ng )));

for ts = t_seeds
    [~, idx] = max(reshape(Et(:,:,ts,:), [], 1));
    [f0, k0, r0] = ind2sub([Nf, Nk, Nr], idx);
    
    ccur = zeros(Nt, Nr, 2);   % (f_idx, k_idx)
    ecur = 0;

    ccur(ts, r0, 1) = f0;
    ccur(ts, r0, 2) = k0;
    ecur = ecur + Et(f0, k0, ts, r0);
    
    % -------------------------
    for r = (r0+1):Nr
        prev = squeeze(ccur(ts, r-1, :))';
        if all(prev==0), prev = [f0, k0]; end
        bestv = -Inf; bestfk = prev;
        frange = max(1,prev(1)-da_f):min(Nf, prev(1)+da_f);
        krange = max(1,prev(2)-da_k):min(Nk, prev(2)+da_k);
        for fi = frange
            for ki = krange
                val = Et(fi, ki, ts, r);
                val = val - (aux_f*(fi-prev(1))^2 + aux_k*(ki-prev(2))^2);
                if val > bestv
                    bestv = val;
                    bestfk = [fi, ki];
                end
            end
        end
        ccur(ts, r, :) = bestfk;
        ecur = ecur + bestv;
    end
    
    for r = (r0-1):-1:1
        prev = squeeze(ccur(ts, r+1, :))';
        if all(prev==0), prev = [f0, k0]; end
        bestv = -Inf; bestfk = prev;
        frange = max(1,prev(1)-da_f):min(Nf, prev(1)+da_f);
        krange = max(1,prev(2)-da_k):min(Nk, prev(2)+da_k);
        for fi = frange
            for ki = krange
                val = Et(fi, ki, ts, r);
                val = val - (aux_f*(fi-prev(1))^2 + aux_k*(ki-prev(2))^2);
                if val > bestv
                    bestv = val;
                    bestfk = [fi, ki];
                end
            end
        end
        ccur(ts, r, :) = bestfk;
        ecur = ecur + bestv;
    end
    
    % -------------------------
    % forward in time (t > ts)
    for t = ts+1:Nt
        for r = 1:Nr
            prev = squeeze(ccur(t-1, r, :))';
            if all(prev==0), prev = [f0, k0]; end
            bestv = -Inf; bestfk = prev;
            frange = max(1,prev(1)-da_f):min(Nf, prev(1)+da_f);
            krange = max(1,prev(2)-da_k):min(Nk, prev(2)+da_k);
            for fi = frange
                for ki = krange
                    val = Et(fi, ki, t, r);
                    val = val - (aux_f*(fi-prev(1))^2 + aux_k*(ki-prev(2))^2);
                    if val > bestv
                        bestv = val;
                        bestfk = [fi, ki];
                    end
                end
            end
            ccur(t, r, :) = bestfk;
            ecur = ecur + bestv;
        end
    end
    
    % backward in time (t < ts)
    for t = ts-1:-1:1
        for r = 1:Nr
            prev = squeeze(ccur(t+1, r, :))';
            if all(prev==0), prev = [f0, k0]; end
            bestv = -Inf; bestfk = prev;
            frange = max(1,prev(1)-da_f):min(Nf, prev(1)+da_f);
            krange = max(1,prev(2)-da_k):min(Nk, prev(2)+da_k);
            for fi = frange
                for ki = krange
                    val = Et(fi, ki, t, r);
                    val = val - (aux_f*(fi-prev(1))^2 + aux_k*(ki-prev(2))^2);
                    if val > bestv
                        bestv = val;
                        bestfk = [fi, ki];
                    end
                end
            end
            ccur(t, r, :) = bestfk;
            ecur = ecur + bestv;
        end
    end
    
    % update best
    if ecur > e
        e = ecur;
        c = ccur;
    end
end

end
