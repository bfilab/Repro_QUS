function [total_bsc] = computeBSC_poly(f,ph,dist_params,Nume)
%COMPUTEBSC Compute theoretical BSC (polydisperse).

% Reference: Franceschini, E., et al. (2016). "Quantitative
% Characterization of Tissue Microstructure in Concentrated Cell Pellet
% Biophantoms Based on the Structure Factor Model." IEEE Trans Ultrason
% Ferroelectr Freq Control 63(9): 1321-1334.

% INPUTS:
% f = frequency vector
% ph.
%    a = scatterer radius [m]
%    a_std = scatterer radius standard deviation [m] 
%    c1 = speed of sound in scatterer [m/s]
%    c0 = speed of sound (background) [m/s]
%    n = average particle number density
%    gamma = relative acoustic impedance between scatterer and surrounding
%            medium
%    poisson = poisson for scatterer
%    pwire = scatterer density [g/cm^2]
%    ph2o = medium density [g/cm^2]
% dist_params.
%    n_a_dist = number of radius values to generate for distribution
%    dist_params.n_bins = number of bins for histogram
% Nume = number iterations

%% Get distribution of scatterer radius
rnd_a_values        = normrnd(ph.a, ph.a_std, [dist_params.n_a_dist,1]);
[hist_N,hist_edges] = histcounts(rnd_a_values, dist_params.n_bins);
a_weights           = hist_N/length(rnd_a_values);                          % weights for each radius value
a_values            = hist_edges(1:(end-1)) + mean(diff(hist_edges))/2;     % range of radius values

%% Compute BSC for each radius

for i = 1:length(a_values)
   
    ph.a = a_values(i);
    bsc = computeBSC_mono(f,ph,Nume);
    weighted_bsc = a_weights(i)*bsc;
    if i == 1
        total_bsc = weighted_bsc;
    else
        total_bsc = total_bsc + weighted_bsc;
    end
    
end

end

