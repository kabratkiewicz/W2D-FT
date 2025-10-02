close all
clear
clc

addpath("UTILS\")
NFFT_omega = 1024;
NFFT_eta = 1024;
fontsize = 20;
Init_Env(fontsize,1);
Ncomp = [1,2,3,4,5];
threshold = 50;

load("data.mat");

T = size(signal,1);
R = size(signal,2);
signal = signal./max(abs(signal(:)));

t = 0:T-1;
t = t-mean(t);
r = 0:R-1;
r = r-mean(r);

sigma_t = 8;
sigma_r = 8;
hT             = Gab_Gaussian_Window(t, sigma_t, 0, 0, 0);
hR             = Gab_Gaussian_Window(r, sigma_r, 0, 0, 0).';
hRhT           = hR  .* hT;

omega_bins = linspace(-0.5, 0.5, NFFT_omega);
eta_bins = linspace(-0.5, 0.5, NFFT_eta);

signal_original = signal;
F_original = fftshift(fftshift(fft2(signal_original.*hRhT , NFFT_omega,NFFT_eta),1),2);
F_original_max = max(abs(F_original(:)));

figure;
imagesc(omega_bins, eta_bins, db(abs(F_original)),'Interpolation','bilinear');
set(gca,'ydir','normal');
xlabel('Normalized freq. $\eta$')
ylabel('Normalized freq. $\omega$')
colormap("turbo")
c = colorbar;
c.Label.String = 'Magnitude [dB]';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
clim([max(max(db(abs(F_original)))) - threshold, max(max(db(abs(F_original))))])

for i = 1:length(Ncomp)
    [comp, residuals] = CLEAN2D(signal_original, NFFT_omega, NFFT_eta, Ncomp(i));

    aaa = sum(comp,3);
    F_comp = fftshift(fftshift(fft2(aaa.*hRhT , NFFT_omega,NFFT_eta),1),2);
    max_aaa = max(abs(F_comp(:)));
    F_comp = F_comp./max(abs(F_comp(:)));
    F_comp = F_comp.*F_original_max;

    figure;
    imagesc(omega_bins, eta_bins, db(abs(F_comp)),'Interpolation','bilinear');
    set(gca,'ydir','normal');
    xlabel('Normalized freq. $\eta$')
    ylabel('Normalized freq. $\omega$')
    colormap("turbo")
    c = colorbar;
    c.Label.String = 'Magnitude [dB]';
    c.Label.Interpreter = 'latex';
    c.TickLabelInterpreter = 'latex';
    clim([max(max(db(abs(F_original)))) - threshold, max(max(db(abs(F_original))))])

    F_res = fftshift(fftshift(fft2(residuals.*hRhT, NFFT_omega,NFFT_eta),1),2);

    figure;
    imagesc(omega_bins, eta_bins, db(abs(F_res)) ,'Interpolation','bilinear');
    set(gca,'ydir','normal');
    xlabel('Normalized freq. $\eta$')
    ylabel('Normalized freq. $\omega$')
    colormap("turbo")
    c = colorbar;
    c.Label.String = 'Magnitude [dB]';
    c.Label.Interpreter = 'latex';
    c.TickLabelInterpreter = 'latex';
    clim([max(max(db(abs(F_original)))) - threshold, max(max(db(abs(F_original))))])
end