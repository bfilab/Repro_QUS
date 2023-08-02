function [samp] = convert_to_qus_struct(in_data)
% Convert the data structure created by placenta.iq2rf (which calls
% placenta.read_visualsonics) into a structure compatible with the QUS
% processing code
% 2023/1/12 - The current version of the codes assumes there's a single RF
%   frame (and its segmentation data) in the data struct
%
%   Inputs
%       in_data: struct return by placenta.iq2rf. Must contain
%           .del_mm
%           .fs
%           .rf_data
%           .sysParam
%
%   Ouputs
%       samp: struct with fields (same as the loadData function)
%             (1) 3D RF data
%             (2) 3D envelope data
%             (3) Vectors for dimensions of RF data, in meters
%             (4) Image resolution for each dimension, in meters/voxel
%             (5) Speed of sound (used to convert between time of flight and distance)
%             (6) Time of flight vector, seconds
%             (7) Axial sampling frequency, Hz
%             (8) Surface position, in meters

 % Generate the image axes
[axial_vec, lateral_vec,c] = placenta.generate_image_axes(in_data.rf_data,in_data.sysParam,in_data.fs*1e6);

% Get RF data
samp.data = in_data.rf_data;

% Get envelope data
samp.env = abs(hilbert(samp.data));

samp.x = lateral_vec/1000; % lateral axes vector in m
samp.y = 1; % single frame
samp.z = axial_vec/1000; % axial axes vector in m

% Set image resolution (m/voxel)
samp.dx = mean(diff(samp.x));
samp.dy = 1;
samp.dz = mean(diff(samp.z));

% Set sampling frequency
samp.c = c;
samp.t = dist2time(samp.c,samp.z);
samp.fs = in_data.fs*1e6; % fs in the RF data struct is in [MHz], needs to be in [Hz]
% samp.fs = 1/mean(diff(samp.t));

% Set the surface position
% This can be found from sec_struct, if it is a member of the in_data
% struct. Otherwise, return as empty member

if isfield(in_data,'seg_struct')
    samp_surf_seg = in_data.seg_struct.surf_roi.Position;
    samp_surf = min(samp_surf_seg(:,2)); % least distance from transducer to sample surface
    samp.surf_pos = samp_surf/1000; % in [m]
else
    samp.surf_pos = [];
end