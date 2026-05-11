% Author: Karol Abratkiewicz
% karol.abratkiewicz@pw.edu.pl  
% Warsaw University of Technology
% K. Abratkiewicz, "Downsampled Windowed 2D 
% Fourier Transform in Signal 
% Decomposition and Concentration," 
% in IEEE Transactions on Signal Processing

close all
clear
clc

fontsize = 50;
img_max_size = 1;
threshold = 30;

addpath("UTILS/")
addpath("CONC/")
addpath("MEXTR/")
Init_Env(fontsize,img_max_size);

fs = 256;
t = linspace(-0.5, 0.5, fs);

x = exp(1j*2*pi*(100.*t.^2));
y = exp(1j*2*pi*(100.*t.^2));
signal = x.' .* y;

NFFT_omega = 128;  % FFT size in omega
NFFT_eta   = 128;  % FFT size in eta

wT = 25;
wR = 25;
sigma_t = 10;
sigma_r = 10;
hop_t = 4;
hop_r = 4;
shift_t = 0;
shift_r = 0;
%% processing
[W2DFT_original, w] = W2DFT_decim2D_Inv(signal, NFFT_omega, NFFT_eta, wT, wR, sigma_t, sigma_r, hop_t, hop_r, shift_t, shift_r, 0, 0, 0, 0);
%% plotting the results
T    = size(signal,1);
R    = size(signal,2);
cut_r = ceil(R/2/hop_r);
cut_t = ceil(T/2/hop_t);
omega_bins = linspace(-0.5, 0.5, NFFT_omega);
eta_bins = linspace(-0.5, 0.5, NFFT_eta);

%%
figure;
imagesc(omega_bins, eta_bins, db(abs(W2DFT_original(:,:,cut_t, cut_r))),'Interpolation','bilinear');
set(gca,'ydir','normal');
xlabel('Normalized freq. $\eta$')
ylabel('Normalized freq. $\omega$')
colormap(turbo)
c = colorbar;
c.Label.String = 'Magnitude [dB]';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
clim([max(max(db(abs(W2DFT_original(:,:,cut_t,cut_r))))) - threshold, max(max(db(abs(W2DFT_original(:,:,cut_t,cut_r)))))])

clwin_f = 10;  
clwin_k = 10;
lambda = 0;
n_comp = 1; % number of components!
Kf = 10; 
Kk = 10;
[mask, Cs] = ridge_detect_mask4D(W2DFT_original, omega_bins, eta_bins, n_comp, lambda, clwin_f, clwin_k, Kf, Kk);
t = 0:T-1;
r = 0:R-1;
x_rec = IW2DFT_decim_unit(W2DFT_original.*mask, T, R, w, hop_t, hop_r, shift_t, shift_r);
%%

figure;
imagesc(t,r,real(signal),'Interpolation','bilinear');
set(gca,'YDir','normal')
xlabel('$r$ samples')
ylabel('$t$ samples')
colormap('turbo')
c = colorbar;
c.Label.String = 'Amplitude';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';

figure;
imagesc(t,r,imag(signal),'Interpolation','bilinear');
set(gca,'YDir','normal')
xlabel('$r$ samples')
ylabel('$t$ samples')
colormap('turbo')
c = colorbar;
c.Label.String = 'Amplitude';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';

figure;
imagesc(t,r,real(x_rec),'Interpolation','bilinear');
set(gca,'YDir','normal')
xlabel('$r$ samples')
ylabel('$t$ samples')
colormap('turbo')
c = colorbar;
c.Label.String = 'Amplitude';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';

figure;
imagesc(t,r,imag(x_rec),'Interpolation','bilinear');
set(gca,'YDir','normal')
xlabel('$r$ samples')
ylabel('$t$ samples')
colormap('turbo')
c = colorbar;
c.Label.String = 'Amplitude';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';

