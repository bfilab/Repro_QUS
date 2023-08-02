function [bsc] = computeBSC_mono(f,ph,Nume)
%COMPUTEBSC Compute theoretical BSC (monodisperse)

% Reference: Franceschini, E., et al. (2016). "Quantitative
% Characterization of Tissue Microstructure in Concentrated Cell Pellet
% Biophantoms Based on the Structure Factor Model." IEEE Trans Ultrason
% Ferroelectr Freq Control 63(9): 1321-1334.

% INPUTS:
% f = frequency vector
% ph.
%    a = scatterer radius [m]
%    c1 = speed of sound in scatterer [m/s]
%    c0 = speed of sound (background) [m/s]
%    n = average particle number density [m/s]
%    gamma = relative acoustic impedance between scatterer and surrounding
%            medium
%    poisson = poisson for scatterer
%    pwire = scatterer density [g/cm^2]
%    ph2o = medium density [g/cm^2]
% Nume = number iterations

% OUTPUTS:
% bsc = backscatter coefficient in linear scale

%% Compute BSC term

%Compute Volume of Sphere
a = ph.a;
Vs = (4/3)*pi*a^3;

% Compute k
c0 = ph.c0;
k = (2*pi*f)./c0;

% Compute final term
n = ph.n;
gamma = ph.gamma;
BSC_term = n*(k.^4*Vs^2*gamma^2)/(4*pi^2);

%% Compute Faran Sphere Form Factor
poisson = ph.poisson;
c1 = ph.c1;
pwire = ph.pwire;
ph2o = ph.ph2o;
F = calcFF(f',poisson,a,c0,c1,pwire,ph2o,Nume);
F = F';

%% Compute BSC
bsc = BSC_term.*F;

end

