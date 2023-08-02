function [F] = calcFluidFilledSphereFF(f,a,c)
%CALCFLUIDFILLEDSPHEREFF Calculate fluid-filled sphere form factor.
% INPUTS:
%   f = frequency
%   a = scatterer radius
%   c = speed of sound in medium
% OUTPUTS:
%   F = form factor function

% 09/16/2020 (THL): Created

% Reference: Franceschini, E., et al. (2016). "Quantitative
% Characterization of Tissue Microstructure in Concentrated Cell Pellet
% Biophantoms Based on the Structure Factor Model." IEEE Trans Ultrason
% Ferroelectr Freq Control 63(9): 1321-1334.

% Convert to wavenumber
k = freq2wavenumber(f,c);

% Compute form factor
F = ( (3./(2*k*a)) .* sph_bessel_1(2*k*a,1) ).^2;

end

