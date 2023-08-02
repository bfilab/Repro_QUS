function [outputArg1,outputArg2] = fit2FFS_fminsearch(inputArg1,inputArg2)
%FIT2FFS_FMINSEARCH Summary of this function goes here
%   Detailed explanation goes here

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

%% Compute Fluid-Filled Sphere Form Factor
F = calcFluidFilledSphereFF(f,ph.a,ph.c0);
F(isnan(F)) = 0;

%% Compute BSC
bsc = BSC_term.*F;

end

