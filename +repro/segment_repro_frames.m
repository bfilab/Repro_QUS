function [seg_struct] = segment_repro_frames(in_fid,seg_frames,seg_in)
% Manual segmentation of mouse/rat repro frames
% Two segmentations are done:
%   1) contour of the repro
%   2) curve following the skin surface (for attenuation compensation)
%
%   Inputs:
%       in_fid: name of file containing RF echo data. Format should follow
%           the same is repro.iq2rf, which calls read_visualsonics
%   
%       seg_frames: vector specifying which frames to segment. If omitted
%           or empty, segment all frames
%
%   Outputs:
%       seg_struct: structure containing the segmented boundaries. The size
%           of this structure is the same size as the number of frames in
%           the data file, though the only segmentation data populated will
%           be for those frames specified in seg_frames. Has fields
%               p_roi - segmented boundary of repro, stored as Polygon
%                   object
%               surf_roi - segmented boundary of skin surface, stored as
%                   Polyline object
%

%%
% Load the data
load(in_fid);

% If seg_frames is not defined, set it to be the whole stack of frames
if ~exist('seg_frames','var') || isempty(seg_frames)
    seg_frames = 1:size(rf_data,3);
end

% Generate the images axes
% Use acquisition parameters defined in sysParam to maintain consistency
fs = fs*1e6; % Convert from [MHz] to [Hz]
[axial_vec,lateral_vec] = repro.generate_image_axes(rf_data,sysParam,fs);

% Initialize the structure that will contain the segmentation information
num_frames = size(rf_data,3);
seg_struct = struct([]);

% Add an entry for the last frame to ensure struct size matches the number
% of frames
seg_struct(num_frames).p_roi = [];
seg_struct(num_frames).surf_roi = [];

%% 
% Segment the specified frames
for f_idx=seg_frames
    this_rf = rf_data(:,:,f_idx); 
    
    % Display the B-mode image
    bmode_im = 20*log10(abs(hilbert(this_rf)));
    bmode_im = bmode_im - max(bmode_im(:));

    figure(1); clf;
    imagesc(lateral_vec,axial_vec,bmode_im,[-45, 0]);
    colormap gray;
    axis image;

    % Segment the repro first
    p_roi = drawpolygon(gca,'Color','r');

    % Segment the skin surface
    surf_roi = drawpolyline(gca,'Color','y');

    % To make things easier later, interpolate the points from the surface
    % segmentation over the number of A-lines
    s_points = surf_roi.Position;
    s_interp = interp1(s_points(:,1),s_points(:,2),lateral_vec,'linear','extrap');
    surf_roi.Position = [lateral_vec(:), s_interp(:)];

    seg_struct(f_idx).p_roi = p_roi;
    seg_struct(f_idx).surf_roi = surf_roi;

end

% Save the segmentation structure to the .mat file with the RF data
save(in_fid,'seg_struct','-append');