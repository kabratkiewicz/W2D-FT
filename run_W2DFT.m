close all
clear
clc

fontsize = 20;
img_max_size = 1;
threshold = 40;

addpath("UTILS/")
addpath("CONC/")
Init_Env(fontsize,img_max_size);

load("data.mat");  % load example signal - 2D time domain signal

signal = awgn(signal, 10); % add some noise

estimator  = 3;     % estimator number from 1 to 3
NFFT_omega = 1024;  % FFT size in omega
NFFT_eta   = 1024;  % FFT size in eta

%% processing
[W2DFT_distribution, W2DFT_concentrated_distribution] = W2DFT(signal, estimator, NFFT_omega, NFFT_eta); 

%% plotting the results
omega_bins = linspace(-0.5, 0.5, NFFT_omega);
eta_bins = linspace(-0.5, 0.5, NFFT_eta);

figure;
imagesc(omega_bins,eta_bins,db(abs(W2DFT_distribution)))
xlabel('Normalized freq. $\omega$')
ylabel('Normalized freq. $\eta$')
colormap("turbo")
c = colorbar;
c.Label.String = 'Magnitude [dB]';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
clim([max(max(db(abs(W2DFT_distribution)))) - threshold, max(max(db(abs(W2DFT_distribution))))])

figure;
imagesc(omega_bins,eta_bins,db(abs(W2DFT_concentrated_distribution)))
xlabel('Normalized freq. $\omega$')
ylabel('Normalized freq. $\eta$')
colormap("turbo")
c = colorbar;
c.Label.String = 'Magnitude [dB]';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
clim([max(max(db(abs(W2DFT_concentrated_distribution)))) - threshold, max(max(db(abs(W2DFT_concentrated_distribution))))])