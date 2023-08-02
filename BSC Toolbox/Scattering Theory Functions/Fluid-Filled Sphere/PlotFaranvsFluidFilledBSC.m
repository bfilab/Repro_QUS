%% Plot Faran vs Fluid Filled Sphere BSC for Phantoms
% 09/18/2020 (THL): Created

% Add BSC toolbox
addpath(genpath('C:\Users\tlye\OneDrive - Riverside Research\Documents\MATLAB\BSC Toolbox'));

%% Inputs

base_dir = ['C:\Users\tlye\OneDrive - Riverside Research\Documents',...
    '\MATLAB\BSC Toolbox\Scattering Theory Functions'];
load('ph_60');
load('ph_15');
load('ph_18');

%% Compute BSC's

f = linspace(0,60*1e6,5000);

fs_60_bsc = computeBSC_mono(f,ph_60,ph_60.Nume);
ffs_60_bsc = computeBSC_mono_FFS(f,ph_60);

fs_15_bsc = computeBSC_mono(f,ph_15,ph_15.Nume);
ffs_15_bsc = computeBSC_mono_FFS(f,ph_15);

fs_18_bsc = computeBSC_mono(f,ph_18,ph_18.Nume);
ffs_18_bsc = computeBSC_mono_FFS(f,ph_18);

%% Plot
subplot(1,3,1);
plot(f*1e-6,10*log10(fs_60_bsc),'LineWidth',2); hold on;
plot(f*1e-6,10*log10(ffs_60_bsc),'LineWidth',2); hold off;
xlabel('MHz'); ylabel('dB/str/m');
title('60 \mum Phantom');

subplot(1,3,2);
plot(f*1e-6,10*log10(fs_18_bsc),'LineWidth',2); hold on;
plot(f*1e-6,10*log10(ffs_18_bsc),'LineWidth',2); hold off;
xlabel('MHz'); ylabel('dB/str/m');
title('18 \mum Phantom');

subplot(1,3,3);
plot(f*1e-6,10*log10(fs_15_bsc),'LineWidth',2); hold on;
plot(f*1e-6,10*log10(ffs_15_bsc),'LineWidth',2); hold off;
xlabel('MHz'); ylabel('dB/str/m');
title('15 \mum Phantom');