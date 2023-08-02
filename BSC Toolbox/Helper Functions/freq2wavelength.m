function [l] = freq2wavelength(f,c)
%FREQ2WAVELENGTH Convert frequency to wavelength.
% 7/1/2020 (THL): Created

% INPUTS:
%   f = frequency [Hz]
%   c = speed [m/s]
% OUTPUTS:
%   l = wavelength [m] 

l = (1/f)*c;

end

