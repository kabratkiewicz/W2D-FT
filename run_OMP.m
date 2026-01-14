close all
clear
clc

addpath("MEXTR\")
addpath("UTILS\")

load("two_chirps.mat");
signal = x;

T = size(signal,1);
R = size(signal,2);
signal = signal./max(abs(signal(:)));
signal = awgn(signal, 20,"measured"); % SNR = 20 dB
y2D = signal;
Nx = size(signal,2);
Ny = size(signal,1);
NFFT_omega = 128;  % FFT size in omega
NFFT_eta   = 128;  % FFT size in eta
omega_bins = linspace(-0.5, 0.5, NFFT_omega);
eta_bins = linspace(-0.5, 0.5, NFFT_eta);

y = y2D(:);

n1 = (0:Nx-1).'; n1 = n1 - mean(n1);
k1 = 0:Nx-1; k1 = k1 - mean(k1);
D1 = exp(1j*2*pi*n1*k1/Nx) / sqrt(Nx);

n2 = (0:Ny-1).'; n2 = n2 - mean(n2);
k2 = 0:Ny-1; k2 = k2 - mean(k2);
D2 = exp(1j*2*pi*n2*k2/Ny) / sqrt(Ny);

D = kron(D2, D1);

K = 500;   % number of atoms
x_hat = omp_complex(D, y, K);

y_hat = D * x_hat;
y2D_hat = reshape(y_hat, Ny, Nx);

F_original = db(fftshift(fftshift(fft2(signal , NFFT_omega,NFFT_eta),1),2));
F_original = F_original - max(F_original(:));

F_rec = db(fftshift(fftshift(fft2(y2D_hat , NFFT_omega,NFFT_eta),1),2));
F_rec = F_rec - max(F_rec(:));

figure;
imagesc(omega_bins, eta_bins, F_original,'Interpolation','bilinear');
set(gca,'ydir','normal');
xlabel('Normalized freq. $\eta$')
ylabel('Normalized freq. $\omega$')
colormap("turbo")
c = colorbar;
c.Label.String = 'Magnitude [dB]';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
clim([-50 0])
drawnow

figure;
imagesc(omega_bins, eta_bins, F_rec,'Interpolation','bilinear');
set(gca,'ydir','normal');
xlabel('Normalized freq. $\eta$')
ylabel('Normalized freq. $\omega$')
colormap("turbo")
c = colorbar;
c.Label.String = 'Magnitude [dB]';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
clim([-50 0])
drawnow

