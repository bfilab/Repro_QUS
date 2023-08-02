function [sub_vol,roi_idx] = getROI(vol,roi)
%GETROI Extract 3D ROI from 3D volume.

% INPUTS:
%   vol.
%       data = 3D rf data (dim 1 = z, dim 2 = x, dim 3 = y)
%       x = position of each data voxel in x [meters]
%       y = position of each data voxel in y [meters]
%       z = position of each data voxel in z [meters]
%       dx = size of voxel [meters/voxel]
%       dy = size of voxel [meters/voxel]
%       dz = size of voxel [meters/voxel]
%   roi.
%       pos_x = initial position of roi in x [meters] (optional field)
%       pos_y = initial position of roi in y [meters]
%       pos_z = initial position of roi in z [meters]
%       len_x = length of roi in x [meters]           (optional field)
%       len_y = length of roi in x [meters]
%       len_z = length of roi in x [meters]

% if roi.pos_x,y,z are empty, default to 0 meters along dimension
% if roi.len_x,y,z are empty, default to full length of volume

% OUTPUTS:
%   sub_vol.
%           data = 3D rf data in ROI
%           x = position of each voxel in ROI in x [meters]
%           y = position of each voxel in ROI in y [meters]
%           z = position of each voxel in ROI in z [meters]
%   roi_idx.
%           pos_x = initial position of roi in x [indices]
%           pos_y = initial position of roi in y [indices]
%           pos_z = initial position of roi in z [indices]
%           len_x = length of roi in x [indices]
%           len_y = length of roi in x [indices]
%           len_z = length of roi in x [indices]
%           x_vec, y_vec, z_vec = vector of indices to sub into vol.data to
%                                 get sub_vol.data

%% Default values

% Default ROI position is initial point along dimension
if isempty(roi.pos_x)
   roi.pos_x = vol.x(1); 
end
if isempty(roi.pos_y)
   roi.pos_y = vol.y(1); 
end
if isempty(roi.pos_z)
   roi.pos_z = vol.z(1); 
end

% Default length of ROI is full length of volume
if isempty(roi.len_x)
   roi.len_x = size(vol.data,2)*vol.dx; 
end
if isempty(roi.len_y)
   roi.len_y = size(vol.data,3)*vol.dy; 
end
if isempty(roi.len_z)
   roi.len_z = size(vol.data,1)*vol.dz; 
end

%% Convert Lengths to Index Values
roi_idx.len_x = round(roi.len_x /vol.dx);
roi_idx.len_y = round(roi.len_y /vol.dy);
roi_idx.len_z = round(roi.len_z /vol.dz);

%% Find Position
[er_x, roi_idx.pos_x] = min(abs( vol.x - roi.pos_x ));
[er_y, roi_idx.pos_y] = min(abs( vol.y - roi.pos_y ));
[er_z, roi_idx.pos_z] = min(abs( vol.z - roi.pos_z));
if er_x>vol.dx || er_y>vol.dy || er_z>vol.dx
    error('Could not find desired ROI...')
end

%% Extract ROI

% Get ROI index vectors
roi_idx.x_vec = roi_idx.pos_x:(roi_idx.pos_x + roi_idx.len_x - 1);
roi_idx.y_vec = roi_idx.pos_y:(roi_idx.pos_y + roi_idx.len_y - 1);
roi_idx.z_vec = roi_idx.pos_z:(roi_idx.pos_z + roi_idx.len_z - 1);

% Trim ROI if necesssary
if sum(roi_idx.x_vec < 1) || ...                                            % check for negative indices in the x direction
        sum(roi_idx.x_vec > size(vol.data, 2)) || ...                       % check for indices > volume size in the x direction
        sum(roi_idx.y_vec < 1) || ...                                       % similarly for y direction
        sum(roi_idx.y_vec > size(vol.data, 3)) || ...
        sum(roi_idx.z_vec < 1) || ...                                       % similarly for z direction
        sum(roi_idx.z_vec > size(vol.data, 1))                            
    error('ROI was cut-off at: x=%i y=%i z=%i',roi.pos_x,roi.pos_y,roi.pos_z);
end
roi_idx.x_vec( roi_idx.x_vec < 1 | roi_idx.x_vec > size(vol.data, 2)) = [];
roi_idx.y_vec( roi_idx.y_vec < 1 | roi_idx.y_vec > size(vol.data, 3)) = [];
roi_idx.z_vec( roi_idx.z_vec < 1 | roi_idx.z_vec > size(vol.data, 1)) = [];

% Extract Sub-volume
sub_vol.x = vol.x(roi_idx.x_vec);
sub_vol.y = vol.y(roi_idx.y_vec);
sub_vol.z = vol.z(roi_idx.z_vec);
sub_vol.data = vol.data(roi_idx.z_vec, roi_idx.x_vec, roi_idx.y_vec);

end

