function [bw_f_adapt] = findAdaptBandwidth(psd,f,opt)
%FINDADAPTBANDWIDTH Find the adaptive bandwidth.
% INPUTS:
%   psd = power spectral density, in dB
%   f = frequency, in Hz
%   opt
%      .numMx = number of points to take as psd max value
%      .noise_bw = bandwidth of noise (Hz) [freq1, freq2]
%      .minDB_default = default dB bandwidth value (dB)
%      .minDB_SNR = adaptive dB bandwidth value added to SNR (dB), should
%                   be negative
%      .fRange = range in which to detect peak (Hz)
% OUTPUTS:
%   bw_f_adapt = frequencies corresponding to adaptive bandwidth
%                [f1,f2] (Hz)

% 07/09/2020 (THL): Created

%% Unpack inputs
numMx = opt.numMx;
noise_bw = opt.noise_bw;
minDB_default = opt.minDB_default;
minDB_SNR = opt.minDB_SNR;
fRange = opt.fRange;
psd = mean(psd(:,:),2);

%% Get SNR

% Power of signal
[sorted_psd, sorted_psd_idx] = sort(psd,1,'descend');
max_psd_vals = sorted_psd(1:numMx);
max_psd_idx = sorted_psd_idx(1:numMx);
pSig = nanmean(max_psd_vals); 
max_psd_idx = round(nanmean(max_psd_idx, 1));

% Power of noise
[~,noiseIdx] = min(abs(noise_bw-f'),[],1);
pNoise = nanmean( psd(noiseIdx) );

% Get SNR
snr = pSig - pNoise;

%% Get default dB bandwidth
[bw_idx_default,bw_f_default] = finddBBandwidthv2(psd,f,minDB_default,fRange);

%% Get adaptive dB bandwidth
minDB_adapt = snr + minDB_SNR;
[bw_idx_adapt,bw_f_adapt] = finddBBandwidthv2(psd,f,minDB_adapt,fRange);

%% Return appropriate bandwidth


end

