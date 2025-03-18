% Author: Karol Abratkiewicz
% karol.abratkiewicz@pw.edu.pl  
% Warsaw University of Technology
% K. Abratkiewicz, "Windowed Two-Dimensional Fourier Transform 
% Concentration and Its Application to ISAR Imaging," in IEEE Transactions 
% on Image Processing, vol. 32, pp. 6260-6273, 2023, 
% doi: 10.1109/TIP.2023.3330603. 

close all
clear
clc

fontsize = 20;
img_max_size = 1;
threshold = 100;

addpath("UTILS/")
addpath("CONC/")
Init_Env(fontsize,img_max_size);

load("data.mat");  % load example signal - 2D time domain signal
signal = awgn(signal, 50,"measured"); % add some noise

est = 1; %estimator

NFFT_omega = 256;  % FFT size in omega
NFFT_eta   = 256;  % FFT size in eta

%% processing
[W2DFT_original, W2DFT_con] = W2DFT_full2Dcon(signal, est, NFFT_omega, NFFT_eta); 
%% plotting the results
T    = size(signal,2);
R    = size(signal,1);
cut_t = ceil(T/2);
cut_r = ceil(R/2);
omega_bins = linspace(-0.5, 0.5, NFFT_omega);
eta_bins = linspace(-0.5, 0.5, NFFT_eta);

figure;
imagesc(omega_bins, eta_bins, db(abs(W2DFT_original(:,:,cut_t,cut_r))));
set(gca,'ydir','normal');
xlabel('Normalized freq. $\omega$')
ylabel('Normalized freq. $\eta$')
colormap("turbo")
c = colorbar;
c.Label.String = 'Magnitude [dB]';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
clim([max(max(db(abs(W2DFT_original(:,:,cut_t,cut_r))))) - threshold, max(max(db(abs(W2DFT_original(:,:,cut_t,cut_r)))))])

figure;
imagesc(omega_bins, eta_bins, db(abs(W2DFT_con(:,:,cut_t,cut_r))));
set(gca,'ydir','normal');
xlabel('Normalized freq. $\omega$')
ylabel('Normalized freq. $\eta$')
colormap("turbo")
c = colorbar;
c.Label.String = 'Magnitude [dB]';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
clim([max(max(db(abs(W2DFT_con(:,:,cut_t,cut_r))))) - threshold, max(max(db(abs(W2DFT_con(:,:,cut_t,cut_r)))))])
%%

[~, W2DFT_concentrated_distribution] = W2DFT(signal, est, NFFT_omega, NFFT_eta); 
figure;
imagesc(omega_bins,eta_bins,db(abs(W2DFT_concentrated_distribution)))
set(gca,'ydir','normal');
xlabel('Normalized freq. $\omega$')
ylabel('Normalized freq. $\eta$')
colormap("turbo")
c = colorbar;
c.Label.String = 'Magnitude [dB]';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
clim([max(max(db(abs(W2DFT_concentrated_distribution)))) - threshold, max(max(db(abs(W2DFT_concentrated_distribution))))])
