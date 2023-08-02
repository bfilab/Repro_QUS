%% Test Fitting of Faran Sphere BSC
% 09/23/2020 (THL): Created

clear; close all; clc;
addpath(genpath(['C:\Users\tlye\OneDrive - Riverside Research',...
    '\Documents\MATLAB\BSC Toolbox']));


%% Inputs

ph_fname = 'ph_18';
arange = [5:30]*1e-6;
f = linspace(25*1e6,60*1e6,2000);

%% Create Faran Sphere BSC
ph_dir = ['C:\Users\tlye\OneDrive - Riverside Research\Documents\',...
    'MATLAB\BSC Toolbox\Scattering Theory Functions'];
load(ph_fname);

if strcmp(ph_fname,'ph_18')
    ph = ph_18;
elseif strcmp(ph_fname,'ph_60')
    ph = ph_60;
elseif strcmp(ph_fname,'ph_15')
    ph = ph_15;
end

bsc = computeBSC_mono(f,ph,ph.Nume);

%% Fit
[bsc_faran,a_star] = fitFaranSphere_simple(bsc,f,arange,ph);