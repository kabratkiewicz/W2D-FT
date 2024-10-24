function [ W ] = blackman_harris_window( len, D, T) 
% len - window length, 
% D - derivative order
% T - time ramp order

    a_0 = 0.3232153788877343;
    a_1 = -0.4714921439576260;
    a_2 = 0.175534129901972;
    a_3 = -0.02849699010614994;
    a_4 = 0.001261357088292677;
    
    n = 0:len-1;
    W = (1.0 * 2 * pi)^D .* a_1 .* cos(1.0 * 2 * pi * n / len + 0.5 * pi * D) + ...
        (2.0 * 2 * pi)^D .* a_2 .* cos(2.0 * 2 * pi * n / len + 0.5 * pi * D) + ...
        (3.0 * 2 * pi)^D .* a_3 .* cos(3.0 * 2 * pi * n / len + 0.5 * pi * D) + ...
        (4.0 * 2 * pi)^D .* a_4 .* cos(4.0 * 2 * pi * n / len + 0.5 * pi * D);
    if D == 0
        W = W + a_0;
    end

    if D > 0
        W = W ./ (len^D);
    end

    if T > 0
        t = -len/2:len/2-1;
        W = W.* (t.^T);
    end
    
end
