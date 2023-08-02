function [axial_vec, lateral_vec,c,elev_vec] = generate_image_axes(rf_data,sysParam,fs)
% Generate the axial and lateral image vectors based on the size of the RF
% frames and the acquisition parameters defined in sysParam. Using this
% function will ensure consistency during data processing
%
%   Inputs:
%       rf_data: 2D or 3D array of RF echo data
%       sysParam: struct containing acquisition parameters, returned by
%           placenta.read_visuasonics
%       fs: sampling frequency computed during IQ to RF conversion [Hz]
%
%   Output:
%       axial_vec: axial axes vector [mm]
%       lateral_vec: lateral axes vector [mm]
%       c: speed of sound defined during Vevo2100 acquisition [m/s]
%       elev_vec: elevational axes vector. Set to 1 if single frame.
%           Otherwise, generate based on the step size specified in
%           sysParam

% Get the depth offset, speed of sound,and image width from sysParam
c = str2double(sysParam.Sound_Speed_Tissue.value)/1000; % convert to [m/s]
del_mm = str2double(sysParam.B_Mode_Depth_Offset.value); % already in [mm]
im_width = str2double(sysParam.B_Mode_Width.value); % already in [mm]

axial_vec = [0:size(rf_data,1)-1]*c/(2*fs)*1000 + del_mm; % in [mm]
lateral_vec = linspace(0,im_width,size(rf_data,2));

if size(rf_data,3) == 1
    elev_vec = 1;
else
    % 3D data block. Get the elevational step size
    elev_step = str2double(sysParam.t3D_Step_Size.value);
    elev_dist = str2double(sysParam.t3D_Scan_Distance.value);

    elev_vec = 0:elev_step:elev_dist-elev_step;
end
