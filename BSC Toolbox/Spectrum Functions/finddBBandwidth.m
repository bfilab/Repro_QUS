function [dB_i_1,dB_i_2] = finddBBandwidth(psd,n)
%FINDDBBANDWIDTH Find the -n dB bandwidth. Assumes the psd has a generally
%concave shape.
% INPUTS:
%   psd: power spectral density, in dB
%   n: -n dB bandwidth
% OUTPUTS:
%   dB_i_1: index into frequency vector of psd corresponding to lower limit
%   of -n dB bandwidth.
%   dB_i_2: index into frequency vector of psd corresponding to upper limit
%   of -n dB bandwidth.
% MODIFICATION HISTORY:
%   2/7/2020 (THL): Created
%   3/5/2020 (THL): Moved conversion to dB to the top of the code to make
%                   it easy to remove or add back in.

% Convert to dB
% psd = 10*log10(psd);

% Find peak value of the psd
[peak_val,peak_i] = max(psd);

% Get -n dB value
n = abs(n);
dB_val = peak_val - n;

% Find indices corresponding to -n dB value
[~,dB_i_1] = min(abs(dB_val-psd(1:peak_i)));
[~,dB_i_2] = min(abs(dB_val-psd(peak_i:end)));
dB_i_2 = dB_i_2 + peak_i - 1;

end

