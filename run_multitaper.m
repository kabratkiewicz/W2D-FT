close all
clear
clc

fontsize = 20;
img_max_size = 1;
threshold = 40;

addpath("UTILS/")
addpath("MULTAP/")
addpath("HERM/")
Init_Env(fontsize,img_max_size);

load("data.mat");  % load example signal - 2D time domain signal

signal = awgn(signal, 5); % add some noise

NFFT_omega = 1024;  % FFT size in omega
NFFT_eta   = 1024;  % FFT size in eta

M = 4 ; % Hermite function order 
sigmaT = 5; % window time spread parameter for t 
sigmaR = 5; % window time spread parameter for r

%% processing - mutitaper with concentration
[W2DFT_distributions, W2DFT_hermite_mean, W2DFT_concentrated_distributions] = ...
    BivariateMultitaper(signal, M, sigmaT, sigmaR, NFFT_omega, NFFT_eta, 'noncoherent');


%% plotting the results
omega_bins = linspace(-0.5, 0.5, NFFT_omega);
eta_bins = linspace(-0.5, 0.5, NFFT_eta);

% inital image (without concentration)
figure;
imagesc(omega_bins,eta_bins,db(abs(W2DFT_distributions(:,:,1))))
xlabel('Normalized freq. $\omega$')
ylabel('Normalized freq. $\eta$')
colormap("turbo")
c = colorbar;
c.Label.String = 'Magnitude [dB]';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
clim([max(max(db(abs(W2DFT_distributions(:,:,1))))) - threshold, max(max(db(abs(W2DFT_distributions(:,:,1)))))])

% multitaper concetrated image
figure;
imagesc(omega_bins,eta_bins,db(abs(W2DFT_hermite_mean)))
xlabel('Normalized freq. $\omega$')
ylabel('Normalized freq. $\eta$')
colormap("turbo")
c = colorbar;
c.Label.String = 'Magnitude [dB]';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
clim([max(max(db(abs(W2DFT_hermite_mean)))) - threshold, max(max(db(abs(W2DFT_hermite_mean))))])


%% processing - multitaper without concentration
[~, wRwT_mean] = BivariateMultitaperNoConcentration(signal, M,sigmaT, sigmaR, NFFT_omega, NFFT_eta);

figure;
imagesc(omega_bins,eta_bins,db(abs(wRwT_mean)))
xlabel('Normalized freq. $\omega$')
ylabel('Normalized freq. $\eta$')
colormap("turbo")
c = colorbar;
c.Label.String = 'Magnitude [dB]';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
clim([max(max(db(abs(wRwT_mean)))) - threshold, max(max(db(abs(wRwT_mean))))])
%%
