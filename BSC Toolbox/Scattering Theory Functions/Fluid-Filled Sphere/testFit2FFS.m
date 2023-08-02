%% Test Fitting Function for Fluid Filled Sphere
% 09/16/2020 (THL):

clear; close all; clc;

%% Inputs

% Directory with phantom parameters
ph_dir = ['C:\Users\tlye\OneDrive - Riverside Research\',...
    'Documents\MATLAB\BSC Toolbox\Scattering Theory Functions'];
% Choose phantom to use
ph_fname = 'ph_18';

% Range of radii to test
arange = [5:30]*1e-6;

%% Setup phantom parameters

% Phantom parameters
load(fullfile(ph_dir,ph_fname));
if strcmp(ph_fname,'ph_60')
    ph = ph_60;
elseif strcmp(ph_fname,'ph_18')
    ph = ph_18;
elseif strcmp(ph_fname,'ph_15')
    ph = ph_15;
end

% Frequency
f = linspace(25*1e6,60*1e6,1000);

%% Compute Faran Sphere BSC
bsc = computeBSC_mono(f,ph,ph.Nume);

%% Fit to Fluid Filled Sphere BSC

[bsc_ffs,params_ffs] = fit2FFS_simple(bsc,f,arange);

%% Plot

figure(1);

subplot(2,1,1)
plot(f*1e-6,10*log10(bsc),'LineWidth',2);
hold on;
plot(f*1e-6,10*log10(bsc_ffs),'LineWidth',2);
hold off;
axis tight;
xlabel('MHz'); ylabel('dB');

subplot(2,1,2)
plot(f*1e-6,10*log10(bsc)-10*log10(bsc_ffs),'g','LineWidth',2);
axis tight;
xlabel('MHz'); ylabel('dB');

sgtitle(['a = ', num2str(params_ffs.a*1e6), ' um']);

%{
a_range = (1:100)*1e-6;
error = nan(size(a_range));

for i = 1:length(a_range)
    
    % Compute Fluid Filled Sphere BSC
    ph_new = ph;
    ph_new.a = a_range(i);
    bsc_new = computeBSC_mono_FFS(f,ph_new);
    
    % Error
    error(i) = mean((bsc_new - bsc).^2);

    % Plot
    figure(1)
    plot(f*1e-6,10*log10(bsc),'LineWidth',2); hold on;
    plot(f*1e-6,10*log10(bsc_new),'LineWidth',2); hold off;
    pause(0.1)
    
end
%}

