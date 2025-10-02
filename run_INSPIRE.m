% Author: Karol Abratkiewicz
% karol.abratkiewicz@pw.edu.pl  
% Warsaw University of Technology
% K. Abratkiewicz, "A Robust CLEAN-like Approach to 
% ISAR Image Decomposition with Non-stationary Signal 
% Components," in IEEE Transactions 
% on Geoscience and Remote Sensing

close all
clear
clc

fontsize = 20;
img_max_size = 1;
threshold = 50;

addpath("UTILS/")
addpath("CONC/")
addpath("MEXTR/")
Init_Env(fontsize,img_max_size);

load("two_chirps.mat");
signal = x;

est = 3; % estimator

NFFT_omega = 128;  % FFT size in omega
NFFT_eta   = 128;  % FFT size in eta

signal = awgn(signal,20,"measured");

sigma_t = 8;
sigma_r = 8;
%% processing
[W2DFT_original, W2DFT_con] = W2DFT_full2Dcon_Inv(signal, est, NFFT_omega, NFFT_eta, sigma_t, sigma_r); 
%% plotting the results
T    = size(signal,1);
R    = size(signal,2);
cut_r = ceil(R/2);
cut_t = ceil(T/2);
omega_bins = linspace(-0.5, 0.5, NFFT_omega);
eta_bins = linspace(-0.5, 0.5, NFFT_eta);

%%
figure;
imagesc(omega_bins, eta_bins, db(abs(W2DFT_original(:,:,cut_t,cut_r))),'Interpolation','bilinear');
set(gca,'ydir','normal');
xlabel('Normalized freq. $\eta$')
ylabel('Normalized freq. $\omega$')
colormap("turbo")
c = colorbar;
c.Label.String = 'Magnitude [dB]';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
clim([max(max(db(abs(W2DFT_original(:,:,cut_t,cut_r))))) - threshold, max(max(db(abs(W2DFT_original(:,:,cut_t,cut_r)))))])
%%
figure;
imagesc(omega_bins, eta_bins, db(abs(W2DFT_con(:,:,cut_t,cut_r))),'Interpolation','bilinear');
set(gca,'ydir','normal');
xlabel('Normalized freq. $\eta$')
ylabel('Normalized freq. $\omega$')
colormap("turbo")
c = colorbar;
c.Label.String = 'Magnitude [dB]';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
clim([max(max(db(abs(W2DFT_con(:,:,cut_t,cut_r))))) - threshold, max(max(db(abs(W2DFT_con(:,:,cut_t,cut_r)))))])
% clim([-50 0])

%%
clwin_f = 4;  
clwin_k = 4;
lambda = 0;
n_comp = 2; % number of components!
Kf = 4; 
Kk = 4;
[mask, Cs] = ridge_detect_mask4D(W2DFT_con, omega_bins, eta_bins, n_comp, lambda, clwin_f, clwin_k, Kf, Kk);
[Y] = W2DFT_comp_ext4D(W2DFT_con, mask, sigma_t, sigma_r);
%%
residuals = signal;
rec_signal = zeros(size(signal));

for i = 1:n_comp
    rec_comp = squeeze(Y(i,:,:));
    residuals = residuals - squeeze(Y(i,:,:));
    rec_signal = rec_signal + squeeze(Y(i,:,:));
    [W2DFT_original, ~] = W2DFT_full2Dcon_Inv(rec_comp, 1, NFFT_omega, NFFT_eta, sigma_t, sigma_r); 
    figure;
    imagesc(omega_bins, eta_bins, db(abs(W2DFT_original(:,:,cut_t,cut_r))),'Interpolation','bilinear');
    set(gca,'ydir','normal');
    xlabel('Normalized freq. $\eta$')
    ylabel('Normalized freq. $\omega$')
    colormap("turbo")
    c = colorbar;
    c.Label.String = 'Magnitude [dB]';
    c.Label.Interpreter = 'latex';
    c.TickLabelInterpreter = 'latex';
    clim([max(max(db(abs(W2DFT_original(:,:,cut_t,cut_r))))) - threshold, max(max(db(abs(W2DFT_original(:,:,cut_t,cut_r)))))])

    [W2DFT_original, ~] = W2DFT_full2Dcon_Inv(residuals, 1, NFFT_omega, NFFT_eta, sigma_t, sigma_r); 
    figure;
    imagesc(omega_bins, eta_bins, db(abs(W2DFT_original(:,:,cut_t,cut_r))),'Interpolation','bilinear');
    set(gca,'ydir','normal');
    xlabel('Normalized freq. $\eta$')
    ylabel('Normalized freq. $\omega$')
    colormap("turbo")
    c = colorbar;
    c.Label.String = 'Magnitude [dB]';
    c.Label.Interpreter = 'latex';
    c.TickLabelInterpreter = 'latex';
    clim([max(max(db(abs(W2DFT_original(:,:,cut_t,cut_r))))) - threshold, max(max(db(abs(W2DFT_original(:,:,cut_t,cut_r)))))])

    [W2DFT_original, ~] = W2DFT_full2Dcon_Inv(rec_signal, 1, NFFT_omega, NFFT_eta, sigma_t, sigma_r); 
    figure;
    imagesc(omega_bins, eta_bins, db(abs(W2DFT_original(:,:,cut_t,cut_r))),'Interpolation','bilinear');
    set(gca,'ydir','normal');
    xlabel('Normalized freq. $\eta$')
    ylabel('Normalized freq. $\omega$')
    colormap("turbo")
    c = colorbar;
    c.Label.String = 'Magnitude [dB]';
    c.Label.Interpreter = 'latex';
    c.TickLabelInterpreter = 'latex';
    clim([max(max(db(abs(W2DFT_original(:,:,cut_t,cut_r))))) - threshold, max(max(db(abs(W2DFT_original(:,:,cut_t,cut_r)))))])
    drawnow
end

%%
t = 0:T-1;
r = 0:R-1;
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
imagesc(t,r,real(squeeze(Y(1,:,:))),'Interpolation','bilinear');
set(gca,'YDir','normal')
xlabel('$r$ samples')
ylabel('$t$ samples')
colormap('turbo')
c = colorbar;
c.Label.String = 'Amplitude';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';

figure;
imagesc(t,r,imag(squeeze(Y(1,:,:))),'Interpolation','bilinear');
set(gca,'YDir','normal')
xlabel('$r$ samples')
ylabel('$t$ samples')
colormap('turbo')
c = colorbar;
c.Label.String = 'Amplitude';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';

figure;
imagesc(t,r,real(squeeze(Y(2,:,:))),'Interpolation','bilinear');
set(gca,'YDir','normal')
xlabel('$r$ samples')
ylabel('$t$ samples')
colormap('turbo')
c = colorbar;
c.Label.String = 'Amplitude';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';

figure;
imagesc(t,r,imag(squeeze(Y(2,:,:))),'Interpolation','bilinear');
set(gca,'YDir','normal')
xlabel('$r$ samples')
ylabel('$t$ samples')
colormap('turbo')
c = colorbar;
c.Label.String = 'Amplitude';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';

