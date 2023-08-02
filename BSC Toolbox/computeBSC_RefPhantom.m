function [samp,ref] = computeBSC_RefPhantom(samp,ref,roi)
%COMPUTEBSC_REFPHANTOM Compute BSC using the reference phantom approach.

% INPUTS:
%   samp/ref
%          .fs = sampling frequency [Hz]
%          .N_fft = zero-padding for FFT
%          .sub_vol.data = 3D RF data volume within roi
%          .surf_pos = surface position [meters]
%          .attn_params = attenuation parameters
%          .transm_params = transmission membrane parameters
%          .transm_db = transmission membrane curve [dB]
%          .theo_bsc_db (ref only) = theoretical BSC
%   roi
%      .len_x/y/z = length of ROI [meters]
%      .pos_x/y/z = position of ROI [meters]

% OUTPUTS:
%   samp/ref
%          .f = frequency [Hz]
%          .sub_vol.PS_noLog = power spectrum for ROI [linear scale]
%          .sub_vol.avg_PS_db = power spectrum for ROI [dB]
%          .sub_vol.attn_db = attenuation spectrum for ROI [dB]
%          .sub_vol.est_bsc_db (samp only) = estimated BSC [dB]

% 02/17/2021 (THL): Made attenuation and transmission compensation optional

%% Compute average power spectrum in ROI

% Sample
[samp.f,~,~,samp.sub_vol.PS_noLog] = getPowerSpec(hanning3D(samp.sub_vol.data), samp.fs, samp.N_fft);
samp.sub_vol.avg_PS_db = 10*log10( nanmean(samp.sub_vol.PS_noLog(:,:),2) );

% Reference
[ref.f,~,~,ref.sub_vol.PS_noLog]   = getPowerSpec(hanning3D(ref.sub_vol.data), ref.fs, ref.N_fft);
ref.sub_vol.avg_PS_db =  10*log10( nanmean(ref.sub_vol.PS_noLog(:,:),2) );

%% Compute the bandwidth limits for adaptive bandwidth
% Only do so if the flag is set in the samp structure
if isfield(samp,'adaptive_bw') && samp.adaptive_bw == true
    % Base the bandwidth on the reference spectrum
    ref_ps = ref.sub_vol.avg_PS_db;
    ref_ps = ref_ps - max(ref_ps(:));
    
%     keyboard;
    idx_1 = find(ref_ps >= -samp.bw_thresh,1,'first');
    idx_2 = find(ref_ps >= -samp.bw_thresh,1,'last');
    
    samp.min_bw = ref.f(idx_1);
    samp.max_bw = ref.f(idx_2);
end

%% Attenuation compensation

% Sample
if isfield(samp,'attn_params')
    d = abs( samp.surf_pos - roi.pos_z ) + (roi.len_z/2);
    [samp.sub_vol.attn_db, ~] = compAttn( samp.f, samp.attn_params, d );
else
    samp.sub_vol.attn_db = zeros(size(samp.sub_vol.avg_PS_db));
end

% Reference
if isfield(ref,'attn_params')
    d = abs( ref.surf_pos - roi.pos_z ) + (roi.len_z/2);
    [ref.sub_vol.attn_db, ~]  = compAttn( ref.f, ref.attn_params, d );
else
    ref.sub_vol.attn_db = zeros(size(samp.sub_vol.avg_PS_db));
end

%% Transmmission compensation
if ~isfield(samp,'transm_db')
    samp.transm_db = zeros(size(samp.sub_vol.avg_PS_db));
end
if ~isfield(ref,'transm_db')
    ref.transm_db = zeros(size(ref.sub_vol.avg_PS_db));
end

%% Sample BSC Estimation

% Make sure the vectors are all in the same direction
ref.sub_vol.avg_PS_db = reshape(ref.sub_vol.avg_PS_db, size(samp.sub_vol.avg_PS_db));
samp.sub_vol.attn_db  = reshape(samp.sub_vol.attn_db,  size(samp.sub_vol.avg_PS_db));
ref.sub_vol.attn_db   = reshape(ref.sub_vol.attn_db,   size(samp.sub_vol.avg_PS_db));
ref.theo_bsc_db       = reshape(ref.theo_bsc_db,       size(samp.sub_vol.avg_PS_db));
samp.transm_db        = reshape(samp.transm_db,        size(samp.sub_vol.avg_PS_db));
ref.transm_db         = reshape(ref.transm_db,         size(samp.sub_vol.avg_PS_db));
   
% Calculate
samp.sub_vol.est_bsc_db = (samp.sub_vol.avg_PS_db + samp.sub_vol.attn_db + samp.transm_db) - ...
                          (ref.sub_vol.avg_PS_db  + ref.sub_vol.attn_db  + ref.transm_db) + ...
                          ref.theo_bsc_db ;

end

