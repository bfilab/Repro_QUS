function [roi,roi_idx] = genROIpos(vol,roi)
%GENROIPOS Generate ROI positions.

% INPUTS:
%   vol.
%       x = position of each data voxel in x [meters]
%       y = position of each data voxel in y [meters]
%       z = position of each data voxel in z [meters]
%       dx = size of voxel [meters/voxel]
%       dy = size of voxel [meters/voxel]
%       dz = size of voxel [meters/voxel]
%   roi.
%       len_x = length of roi in x [meters]
%       len_y = length of roi in x [meters]
%       len_z = length of roi in x [meters]
%       o_x = overlap of ROI in x [%, e.g. 50]
%       o_y = overlap of ROI in y [%]
%       o_z = overlap of ROI in z [%]
%       init_x = starting point for roi in x [meters] (field optional)
%       init_y = starting point for roi in y [meters]
%       init_z = starting point for roi in z [meters]
%       end_x = ending point for  roi in x [meters]   (field optional)
%       end_y = ending point for  roi in y [meters]
%       end_z = ending point for  roi in z [meters] 

% OUTPUTS:
%   roi.         (appends to input struct)
%       step_x = distance between adjacent ROI [meters]
%       step_y = distance between adjacent ROI [meters]
%       step_z = distance between adjacent ROI [meters] 
%       pos_x = initial position of roi in x [meters]
%       pos_y = initial position of roi in y [meters]
%       pos_z = initial position of roi in z [meters]
%       init_x, init_y, init_z are added if not initially input
%       end_x, end_y, end_z are added if not initially input
%   roi_idx.
%       (these are the same as in roi, except in indices not meters)
%       init_x, init_y, init_z
%       end_x,  end_y,  end_z
%       len_x,  len_y,  len_z
%       step_x, step_y, step_z
%       pos_x,  pos_y,  pos_z

% 06/09/2020 (THL): Created
% 09/29/2020 (THL): Updated how roi_idx.pos is calculated. Previously, it
%                   included ROI that would be hanging off the edge of the
%                   image. Now it does not, so that it matches with
%                   roi.pos.
% 01/05/2021 (THL): Updated again so that we remove hanging ROI's but also
%                   are able to compute a single ROI.
% 01/14/2021 (THL): Changed roi.pos_x = roi_idx.pos_x*vol.dx to 
%                   roi.pos_x = vol.x(roi_idx.pos_x)
% 01/15/2021 (THL): Added a rounding when computing the step indices, so 
%                   roi_idx.step_x = roi_idx.len_x * (1-roi.o_x/100) to
%                   roi_idx.step_x = round(roi_idx.len_x * (1-roi.o_x/100));

%% Initialize values
% ROI will only be created between the roi.init and roi.end points.

% Default initial point is the start of the dimension
if ~isfield(roi,'init_x')
    roi.init_x = vol.x(1);
end
if ~isfield(roi,'init_y')
    roi.init_y = vol.y(1);
end
if ~isfield(roi,'init_z')
    roi.init_z = vol.z(1);
end

% Default ending point is the end of the dimension
if ~isfield(roi,'end_x')
    roi.end_x = vol.x(end);
end
if ~isfield(roi,'end_y')
    roi.end_y = vol.y(end);
end
if ~isfield(roi,'end_z')
    roi.end_z = vol.z(end);
end

[~,roi_idx.init_x] = min(abs(roi.init_x - vol.x));                          % Starting point for ROI generation
[~,roi_idx.init_y] = min(abs(roi.init_y - vol.y));
[~,roi_idx.init_z] = min(abs(roi.init_z - vol.z));

[~,roi_idx.end_x] = min(abs(roi.end_x - vol.x));                            % Ending point for ROI generation
[~,roi_idx.end_y] = min(abs(roi.end_y - vol.y));
[~,roi_idx.end_z] = min(abs(roi.end_z - vol.z));

%% Convert Lengths to Index Values
roi_idx.len_x = round(roi.len_x /vol.dx);
roi_idx.len_y = round(roi.len_y /vol.dy);
roi_idx.len_z = round(roi.len_z /vol.dz);
fprintf('ROI Length (X,Y,Z) = (%f, %f, %f) mm\n',...
    roi_idx.len_x*vol.dx*1e3,...
    roi_idx.len_y*vol.dy*1e3,...
    roi_idx.len_z*vol.dz*1e3);

%% Get Step Between Adjacent ROI
% roi.step_x indicates the space between the leftmost edge of the first ROI
% and the leftmost edge of the adjacent ROI in x. Similarly for roi.step_y
% and roi.step_z.

% Visual example,
% consider two overlapping rectangular ROI, adjacent in x:
% |----+---|----+
% |    +   |    +
% |----+---|----+
% *----*
% roi.step_x is the space between asterisks (*)
% Where | indicates the edge of the first ROI
% and + indicates the edge of the adjacent ROI

% Get overlap values
roi_idx.step_x = round(roi_idx.len_x * (1-roi.o_x/100));
roi_idx.step_y = round(roi_idx.len_y * (1-roi.o_y/100));
roi_idx.step_z = round(roi_idx.len_z * (1-roi.o_z/100));

roi.step_x = roi_idx.step_x .* vol.dx;                                      % [m]   
roi.step_y = roi_idx.step_y .* vol.dy;                                      % [m]
roi.step_z = roi_idx.step_z .* vol.dz;                                      % [m]

%% Generate ROI positions (starting position of each ROI)
roi_idx.pos_x   = roi_idx.init_x : roi_idx.step_x : roi_idx.end_x;
roi_idx.pos_y   = roi_idx.init_y : roi_idx.step_y : roi_idx.end_y;
roi_idx.pos_z   = roi_idx.init_z : roi_idx.step_z : roi_idx.end_z;

% Don't let any ROI hang off the edge
bad_idx_x = (roi_idx.pos_x + roi_idx.len_x - 1) > roi_idx.end_x;
roi_idx.pos_x(bad_idx_x) = [];
bad_idx_y = (roi_idx.pos_y + roi_idx.len_y - 1) > roi_idx.end_y;
roi_idx.pos_y(bad_idx_y) = [];
bad_idx_z = (roi_idx.pos_z + roi_idx.len_z - 1) > roi_idx.end_z;
roi_idx.pos_z(bad_idx_z) = [];

roi.pos_x = vol.x(roi_idx.pos_x);                % [m]
roi.pos_y = vol.y(roi_idx.pos_y);                % [m]
roi.pos_z = vol.z(roi_idx.pos_z);                % [m]    

end

