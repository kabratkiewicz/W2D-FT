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

addpath('UTILS/')
addpath("CONC/")
Init_Env(30,1)

%% parameters
fs = 256;

NFFT_omega = 128;
NFFT_eta   = 128;

threshold = 40;

wT = 25;
wR = 25;

sigma_t = 10;
sigma_r = 10;

shift_t = 0;
shift_r = 0;

snr = 0;
hop_t = 4;
hop_r = 4;


%% generate clean signal
t = linspace(-0.5, 0.5, fs);

x = exp(1j*2*pi*(-10.*t + 100.*t.^2));
y = exp(1j*2*pi*(-10.*t + 100.*t.^2));
signal1 = x.' .* y;

x = exp(1j*2*pi*(10.*t - 100.*t.^2));
y = exp(1j*2*pi*(10.*t - 100.*t.^2));
signal2 = x.' .* y;

signal_reference = signal1 + signal2;

signal = awgn(signal_reference, snr, 'measured');

[W2DFT_original, W2DFT_con] = ...
    W2DFT_decim_4Dcon(signal, NFFT_omega,  NFFT_eta, wT, wR,  sigma_t, sigma_r, hop_t, hop_r, shift_t, shift_r);


tr_cut = [round(size(W2DFT_original,3)/4), round(size(W2DFT_original,3)/2), round(size(W2DFT_original,3)*0.75) ];
wn_cut = [20 64 108];

T    = size(signal,1);
R    = size(signal,2);
t    = 1:hop_t:T;
r    = 1:hop_r:R;
cut_r = ceil(size(W2DFT_original,3)/2);
cut_t = ceil(size(W2DFT_original,4)/2);
cut_omega = NFFT_omega/2;
cut_eta = NFFT_eta/2;
omega_bins = linspace(-0.5, 0.5, NFFT_omega);
eta_bins = linspace(-0.5, 0.5, NFFT_eta);

for i = 1:3
    figure;
    subplot(2,2,1)
    imagesc(omega_bins, eta_bins, db(abs(squeeze(W2DFT_original(:,:,tr_cut(i),tr_cut(i)))).^2,"power"),'interpolation','bilinear');
    set(gca,'ydir','normal');
    ylabel('Normalized $\omega$')
    xlabel('Normalized $\eta$')
    colormap("turbo")
    c = colorbar;
    c.Label.String = 'Magnitude [dB]';
    c.Label.Interpreter = 'latex';
    c.TickLabelInterpreter = 'latex';
    clim([max(max(db(abs(W2DFT_original(:,:,tr_cut(i),tr_cut(i))).^2,"power"))) - threshold, max(max(db(abs(W2DFT_original(:,:,tr_cut(i),tr_cut(i))).^2,"power")))])

    subplot(2,2,2)
    imagesc(t, eta_bins, db(abs(squeeze(W2DFT_original(wn_cut(i), :, :, tr_cut(i)))).^2,"power"),'interpolation','bilinear');
    set(gca,'ydir','normal');
    ylabel('Normalized $\omega$')
    xlabel('Space $r$ samples')
    colormap("turbo")
    c = colorbar;
    c.Label.String = 'Magnitude [dB]';
    c.Label.Interpreter = 'latex';
    c.TickLabelInterpreter = 'latex';
    clim([max(max(db(abs(W2DFT_original(wn_cut(i), :, :, tr_cut(i))).^2,"power"))) - threshold, max(max(db(abs(W2DFT_original(wn_cut(i), :, :, tr_cut(i))).^2,"power")))])

    % figure;
    subplot(2,2,3)
    imagesc(t, r, db(abs(squeeze(W2DFT_original(wn_cut(i), wn_cut(i), :, :))).^2,"power"),'interpolation','bilinear');
    set(gca,'ydir','normal');
    xlabel('Space $r$ samples')
    ylabel('Time $t$ samples')
    colormap("turbo")
    c = colorbar;
    c.Label.String = 'Magnitude [dB]';
    c.Label.Interpreter = 'latex';
    c.TickLabelInterpreter = 'latex';
    clim([max(max(db(abs(W2DFT_original(wn_cut(i), wn_cut(i), :, :)).^2,"power"))) - threshold, max(max(db(abs(W2DFT_original(wn_cut(i), wn_cut(i), :, :)).^2,"power")))])

    % figure;
    subplot(2,2,4)
    imagesc(r, omega_bins, db(abs(squeeze(W2DFT_original(:, wn_cut(i), tr_cut(i), :))).^2,"power"),'interpolation','bilinear');
    set(gca,'ydir','normal');
    ylabel('Normalized $\omega$')
    xlabel('Space $r$ samples')
    colormap("turbo")
    c = colorbar;
    c.Label.String = 'Magnitude [dB]';
    c.Label.Interpreter = 'latex';
    c.TickLabelInterpreter = 'latex';
    clim([max(max(db(abs(W2DFT_original(:, wn_cut(i), tr_cut(i), :)).^2,"power"))) - threshold, max(max(db(abs(W2DFT_original(:, wn_cut(i), tr_cut(i), :)).^2,"power")))])

    %%%%%%%%%%%%%%%%%%%

    figure;
    subplot(2,2,1)
    imagesc(omega_bins, eta_bins, db(abs(squeeze(W2DFT_con(:,:,tr_cut(i),tr_cut(i))))),'interpolation','bilinear');
    set(gca,'ydir','normal');
    xlabel('Normalized $\eta$')
    ylabel('Normalized $\omega$')
    colormap("turbo")
    c = colorbar;
    c.Label.String = 'Magnitude [dB]';
    c.Label.Interpreter = 'latex';
    c.TickLabelInterpreter = 'latex';
    clim([max(max(db(abs(W2DFT_con(:,:,tr_cut(i),tr_cut(i)))))) - threshold, max(max(db(abs(W2DFT_con(:,:,tr_cut(i),tr_cut(i))))))])



    % figure;
    subplot(2,2,2)
    imagesc(t, eta_bins, db(abs(squeeze(W2DFT_con(wn_cut(i), :, :, tr_cut(i))))),'interpolation','bilinear');
    set(gca,'ydir','normal');
    ylabel('Normalized $\omega$')
    xlabel('Space $r$ samples')
    colormap("turbo")
    c = colorbar;
    c.Label.String = 'Magnitude [dB]';
    c.Label.Interpreter = 'latex';
    c.TickLabelInterpreter = 'latex';
    clim([max(max(db(abs(W2DFT_con(wn_cut(i), :, :, tr_cut(i)))))) - threshold, max(max(db(abs(W2DFT_con(wn_cut(i), :, :, tr_cut(i))))))])


    % figure;
    subplot(2,2,3)
    imagesc(t, r, db(abs(squeeze(W2DFT_con(wn_cut(i), wn_cut(i), :, :)))),'interpolation','bilinear');
    set(gca,'ydir','normal');
    ylabel('Time $t$ samples')
    xlabel('Space $r$ samples')
    colormap("turbo")
    c = colorbar;
    c.Label.String = 'Magnitude [dB]';
    c.Label.Interpreter = 'latex';
    c.TickLabelInterpreter = 'latex';
    clim([max(max(db(abs(W2DFT_con(wn_cut(i), wn_cut(i), :, :))))) - threshold, max(max(db(abs(W2DFT_con(wn_cut(i), wn_cut(i), :, :)))))])



    % figure;
    subplot(2,2,4)
    imagesc(r, omega_bins, db(abs(squeeze(W2DFT_con(:, wn_cut(i), tr_cut(i), :)))),'interpolation','bilinear');
    set(gca,'ydir','normal');
    ylabel('Normalized $\omega$')
    xlabel('Space $r$ samples')
    colormap("turbo")
    c = colorbar;
    c.Label.String = 'Magnitude [dB]';
    c.Label.Interpreter = 'latex';
    c.TickLabelInterpreter = 'latex';
    clim([max(max(db(abs(W2DFT_con(:, wn_cut(i), tr_cut(i), :))))) - threshold, max(max(db(abs(W2DFT_con(:, wn_cut(i), tr_cut(i), :)))))])

end