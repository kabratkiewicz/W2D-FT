function w = Gab_Gaussian_Window(k, sigma, d_order, t_order, CR)
%% Input args:
% k - time samples of the window
% sigma - time spread of the Gaussian window
% d_order - window derivative order
% t_order - time ramp order
% CR - chirp rate of the window

if ~exist('sigma', 'var')
 sigma = 10;
end
if ~exist('k', 'var')
 error('k range not defined');
end
if ~exist('d_order', 'var')
 d_order = 0;
end
if ~exist('t_order', 'var')
 t_order = 0;
end
if ~exist('CR', 'var')
 CR = 0;
end
%%%%%%%%%%%
A     = 1/(sqrt(2*pi)*sigma);
B     = -1 / (2*sigma^2);
C     = -1j * CR * pi;
k2    = k.^2;
w     = A * exp(B * k2);
if CR ~= 0
    w = w .* exp(C * k2);
end

if d_order == 1
    w = -k./sigma^2 .* w;
    if CR ~= 0
        w = A * exp(B * k2) .* exp(C * k2) .* (2*pi.*CR*1j*k+k./sigma^2);
    end
elseif d_order == 2
    w = -w./sigma^2 + k2 .* w./(sigma^4); 
    if CR ~=0 
        w = A * exp(B * k2) .* exp(C * k2) .* (-2*pi.*CR*1j*k-k./sigma^2).^2 + A * exp(B * k2) .* exp(C * k2) .* (-2*pi.*CR*1j-1/sigma^2);
    end
elseif d_order == 3
    w = 3.*k./sigma^4 .* w - k.^3./sigma^6 .*w;
    if CR~=0
        init_w =  A * exp(B * k2);    
        win_CR = exp(-1j.*pi*CR.*k2);
        w = init_w .* win_CR.*(2*pi*CR*1j*k + k./sigma^2).^3+3.*(-2*pi*CR*1j-1/sigma^2) .* init_w.*win_CR .* (2*pi*CR*1j.*k+k/sigma^2);
    end
elseif d_order == 4
    w = w .* (3*sigma^4-6.*sigma^2.*k2 + k2.^2)./sigma^8;
    if CR ~= 0
        % NOT TESTED !!!!!!!!!
        init_w =  A * exp(B * k2);    
        win_CR = exp(-1j.*pi*CR.*k2);
        w = init_w .* win_CR.*( (-2*pi*CR.*k*1j - k./sigma^2).^4 + 6.*( -2*pi*CR*1j - 1/sigma^2 ) .* (-2*pi*CR.*k*1j - k./sigma^2).^2 + 3.*(-2*pi*CR*1j - 1/sigma^2).^2 );
    end
elseif d_order > 4
    error('Window derivative cannot be defined for the order > 4')
end

if t_order > 0
     w = w .* (k).^t_order;
end

w = w(:);


end


