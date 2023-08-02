function [bw_idx,bw_f] = finddBBandwidthv2(psd,f,minDB,frange)
%FINDDBBANDWIDTHV2 Find the -[minDB] dB bandwidth. Based on Daniel's code.
% INPUTS:
%   psd: power spectral density (dB)
%   f: frequency (Hz)
%   minDB: bandwidth value (e.g. minDB=6 for -6 dB bandwidth)
%   frange: [f1, f2] frequency range in which to search for peak (Hz)
% OUTPUTS:
%   bw_idx: [idx1, idx2] indices corresponding to bandwidth
%   bw_f: [f1, f2] frequency (Hz) corresponding to bandwidth

% 07/09/2020 (THL): Created

%% Get minDB value
minDB = -abs( minDB );

%% Normalize and Crop Spectrum
psd = psd - max(psd);
[~,n1] = min( abs(f-frange(1)) );
[~,n2] = min( abs(f-frange(2)) );
psd_crop = psd(n1:n2);

%% Get bandwidth

% Get indices corresponding to bandwidth borders
psd_binary( psd_crop >= minDB) = 1;
psd_binary( psd_crop <= minDB) = 0;
[~, indF1] = max( psd_binary );
[~, indF2] = max( flip( psd_binary ) ); 
indF2 = length(psd_binary) - indF2 + 1;
indF1 = indF1 + n1 - 1;
indF2 = indF2 + n1 - 1;
bw_idx = [indF1,indF2];

% Get bandwidth frequency values
fmin = f( indF1 );
fmax = f( indF2 );
bw_f = [fmin,fmax];

end

