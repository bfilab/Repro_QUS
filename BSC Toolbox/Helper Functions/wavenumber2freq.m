function [f] = wavenumber2freq(k,c)
%WAVENUMBER2FREQ Convert wavenumber to frequency.
% INPUTS:
%   k = wavenumber [rad/m]
%   c = speed of sound [m/s]
% OUTPUTS:
%   f = frequency [Hz]

% 09/22/2020 (THL): Created

f = (k*c)/(2*pi);  

end

