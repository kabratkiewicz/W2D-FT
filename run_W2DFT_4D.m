% Author: Karol Abratkiewicz
% karol.abratkiewicz@pw.edu.pl  
% Warsaw University of Technology
% K. Abratkiewicz, "Four-Dimensional Reassignment," in IEEE Signal Processing Letters, 2024 
% 
% This script transforms the bivariate time-domain signal to the four-dimensional
% spatial time-space-frequency-wavenumber distribution

close all
clear
clc

fontsize = 20;
img_max_size = 1;
threshold = 50;

addpath("UTILS/")
addpath("CONC/")
Init_Env(fontsize,img_max_size);

% load("data.mat");  % load example signal - 2D time domain signal
% signal = awgn(signal, 50,"measured"); % add some noise

NFFT_omega = 256;  % FFT size in omega
NFFT_eta   = 256;  % FFT size in eta

%% processing
[W2DFT_original, W2DFT_con] = W2DFT_4Dcon(signal, NFFT_omega, NFFT_eta); 
%% plotting the results
T    = size(signal,2);
R    = size(signal,1);
t    = 1:T;
r    = 1:R;
cut_t = ceil(T/2);
cut_r = ceil(R/2);
cut_omega = NFFT_omega/2;
cut_eta = NFFT_eta/2;
omega_bins = linspace(-0.5, 0.5, NFFT_omega);
eta_bins = linspace(-0.5, 0.5, NFFT_eta);

figure;
imagesc(omega_bins, eta_bins, db(abs(W2DFT_original(:,:,cut_t,cut_r))));
set(gca,'ydir','normal');
ylabel('Normalized frequency $\eta$')
xlabel('Normalized frequency $\omega$')
colormap("turbo")
c = colorbar;
c.Label.String = 'Magnitude [dB]';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
clim([max(max(db(abs(W2DFT_original(:,:,cut_t,cut_r))))) - threshold, max(max(db(abs(W2DFT_original(:,:,cut_t,cut_r)))))])

figure;
imagesc(t, r, db(abs(squeeze(W2DFT_original(cut_eta,cut_omega,:,:)))));
set(gca,'ydir','normal');
ylabel('Space $r$ samples')
xlabel('Time $t$ samples')
colormap("turbo")
c = colorbar;
c.Label.String = 'Magnitude [dB]';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
clim([max(max(db(abs(W2DFT_original(cut_eta,cut_omega,:,:))))) - threshold, max(max(db(abs(W2DFT_original(cut_eta,cut_omega,:,:)))))])

%% plotting the results
figure;
imagesc(omega_bins, eta_bins, db(abs(W2DFT_con(:,:,cut_t,cut_r))));
set(gca,'ydir','normal');
ylabel('Normalized frequency $\eta$')
xlabel('Normalized frequency $\omega$')
colormap("turbo")
c = colorbar;
c.Label.String = 'Magnitude [dB]';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
clim([max(max(db(abs(W2DFT_con(:,:,cut_t,cut_r))))) - threshold, max(max(db(abs(W2DFT_con(:,:,cut_t,cut_r)))))])

figure;
imagesc(t, r, db(abs(squeeze(W2DFT_con(cut_eta,cut_omega,:,:)))));
set(gca,'ydir','normal');
ylabel('Space $r$ samples')
xlabel('Time $t$ samples')
colormap("turbo")
c = colorbar;
c.Label.String = 'Magnitude [dB]';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
clim([max(max(db(abs(W2DFT_con(cut_eta,cut_omega,:,:))))) - threshold, max(max(db(abs(W2DFT_con(cut_eta,cut_omega,:,:)))))])

%% concenration only in the frequency domains
[~, W2DFT_concentrated_distribution] = W2DFT_full2Dcon(signal, 1, NFFT_omega, NFFT_eta); 
figure;
imagesc(omega_bins,eta_bins,db(abs(W2DFT_concentrated_distribution(:,:,cut_t,cut_r))))
set(gca,'ydir','normal');
ylabel('Normalized frequency $\eta$')
xlabel('Normalized frequency $\omega$')
colormap("turbo")
c = colorbar;
c.Label.String = 'Magnitude [dB]';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
clim([max(max(db(abs(W2DFT_concentrated_distribution(:,:,cut_t,cut_r))))) - threshold, max(max(db(abs(W2DFT_concentrated_distribution(:,:,cut_t,cut_r)))))])

figure;
imagesc(t, r, db(abs(squeeze(W2DFT_concentrated_distribution(cut_eta,cut_omega,:,:)))));
set(gca,'ydir','normal');
ylabel('Space $r$ samples')
xlabel('Time $t$ samples')
colormap("turbo")
c = colorbar;
c.Label.String = 'Magnitude [dB]';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
clim([max(max(db(abs(W2DFT_concentrated_distribution(cut_eta,cut_omega,:,:))))) - threshold, max(max(db(abs(W2DFT_concentrated_distribution(cut_eta,cut_omega,:,:)))))])




% [W2DFT_distribution, W2DFT_concentrated_distribution] = W2DFT(signal, 2, NFFT_omega, NFFT_eta); 
% figure;
% imagesc(omega_bins,eta_bins,db(abs(W2DFT_concentrated_distribution)))
% set(gca,'ydir','normal');
% xlabel('Normalized freq. $\eta$')
% ylabel('Normalized freq. $\omega$')
% colormap("turbo")
% c = colorbar;
% c.Label.String = 'Magnitude [dB]';
% c.Label.Interpreter = 'latex';
% c.TickLabelInterpreter = 'latex';
% clim([max(max(db(abs(W2DFT_concentrated_distribution)))) - threshold, max(max(db(abs(W2DFT_concentrated_distribution))))])
% 
