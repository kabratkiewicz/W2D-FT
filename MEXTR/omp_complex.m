function x = omp_complex(D, y, K)
residual = y;
idx = [];
x = zeros(size(D,2),1);

for k = 1:K
    proj = abs(D' * residual);   
    [~, i] = max(proj);
    idx = unique([idx i]);

    Ds = D(:,idx);
    xs = pinv(Ds) * y;
    residual = y - Ds*xs;

    if norm(residual) < 1e-6
        break
    end
end

x(idx) = xs;
end
