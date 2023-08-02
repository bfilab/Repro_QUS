function [err] = fluidspherefit_objFcn(x,bsc,f)
%FARANVSFLUID_OBJFCN Error between given BSC and fluid-filled sphere BSC.
% INPUTS:
%   x(1) = scatterer radius [m]
%    (2) = scatterer number density [scatterers/m^3]
%    (3) = relative acoustic impedance between scatterer and medium
%   bsc = backscatter coefficient in 1/str/m (NOT dB).
%   f = frequency vector [Hz]
% OUTPUTS:
%   err = error between fluid sphere BSC and input BSC

% 09/18/2020 (THL): Created

%% Setup
ph.a = x(1);
ph.c0 = 1500;
ph.n = x(2);
ph.gamma = x(3);

%% Compute fluid sphere BSC
bsc_ffs = computeBSC_mono_FFS(f,ph);

%% Compute error
bsc_ffs = reshape(bsc_ffs,size(bsc));
err = mean((bsc - bsc_ffs).^2);

end

