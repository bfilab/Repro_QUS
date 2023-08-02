function [k] = freq2wavenumber(f,c)
%FREQ2WAVENUMBER Convert frequency to wavenumber.
% INPUTS:
%   f = frequency [Hz]
%   c = speed of sound [m/s]
% OUTPUTS:
%   k = wavenumber [rad/m]

% 09/16/2020 (THL): Created

k = (2*pi*f)/c;    

end

