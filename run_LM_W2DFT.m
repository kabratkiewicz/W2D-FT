% Author: Karol Abratkiewicz
% karol.abratkiewicz@pw.edu.pl  
% Warsaw University of technology
% K.Abratkiewicz, M. K. Baczyk, P. Samczynski "Software-Oriented Adjustable
% Concentration of ISAR Images -- Towards Super-Resolution", RadarConf 2025
% 
% This script runs the Lavenberg-Marquad 2D spectrum concentration. The
% input signal is transformed into the bivariate spectrum. Next, the
% distribution is concentrated with a dump parameter mu.

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

mu = [0.001, 0.2 0.5, 1, 1.5 2 5 10];
omega_bins = linspace(-0.5, 0.5, NFFT_omega);
eta_bins = linspace(-0.5, 0.5, NFFT_eta);

for m = 1:numel(mu)
    [S_con] = W2DFT_LM(signal, NFFT_omega, NFFT_eta, mu(m));
    figure;
    imagesc(omega_bins, eta_bins, db(abs(S_con)))
    set(gca, 'YDir','normal')
    xlabel('Normalized freq. $\omega$')
    ylabel('Normalized freq. $\eta$')
    colormap("turbo")
    c = colorbar;
    c.Label.String = 'Magnitude [dB]';
    c.Label.Interpreter = 'latex';
    c.TickLabelInterpreter = 'latex';
    clim([max(max(db(abs(S_con)))) - threshold, max(max(db(abs(S_con))))])
end