function [ref_fid] = select_ref_data(ref_list, sysParam, calib_phantom, surf_dist)
% Function to select the most appropriate reference dataset based on the
% acquisition parameters and distance from transducer to sample surface
%
%   Inputs
%       ref_list: 2D array summarizing the available datasets. Same format
%           that is output from placenta.get_reference_data_list
%
%       sysParam: struct containing acquisition parameters. Format output
%           by placenta.read_visualsonics
%
%       calib_phantom: number indicating the calibration phantom (15 or 18)
%
%       surf_dist: distance from transducer to sample surface, in [mm]
%
%   Outputs
%       ref_fid: filename of reference dataset to use (does not include
%           full path)
%

% Get the name of the transducer and location of transmit focus from the
% acquisition parameters struct
trans = sysParam.Transducer_Name.value;
tx_pos = str2double(sysParam.B_Mode_Focal_Zones_Pos_Array1.value);
%tx_pos = str2double(sysParam.B_Mode_Focal_Zones_Pos1.value);

% Strip the 'LZ-' from the transducer name and convert to a double
trans = str2double(trans(4:end));

% Keep only rows for reference datasets collected with the specified
% transducer and Tx focus position for the chosen calibration phantom
ref_list = ref_list(ref_list(:,1) == trans,:);
ref_list = ref_list(ref_list(:,4) == tx_pos,:);
ref_list = ref_list(ref_list(:,2) == calib_phantom,:);

% Look for the reference dataset where the phantom surface distance is
% closest to the sample surface distance, constrained the the phantom
% surface distance must be less than or equal to the distance to the sample
surf_idx = find(ref_list(:,3) <= surf_dist,1,'last');

ref_list = ref_list(surf_idx,:);

% Recreate the filename for the reference dataset
ref_fid = sprintf('LZ%d_%dumPhantom_%dmmSurf_%dmmTx_RF_Data.mat', ...
    ref_list(1), ref_list(2), ref_list(3), ref_list(4));

