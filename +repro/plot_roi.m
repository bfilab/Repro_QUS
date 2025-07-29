function plot_roi(in_fig, roi_pts)
% Function to plot the ROI defined in roi_pts atop the figure specified by
% the handle in_fig.
% Assume the units of the figure are in mm.
%
% Inputs
%   roi_pts: list of points defining the corners of the ROI, in units of mm
%       The same vertex must appear in the first and last rows to define a 
%       closed shape.

% Ensure the starting point is repeated at the end of the array to close
% the polygon
if ~all(roi_pts(1,:) == roi_pts(end,:))
    roi_pts = cat(1,roi_pts,roi_pts(1,:));
end

% Plot each line segment
for ll=1:size(roi_pts,1)-1
  x1 = roi_pts(ll,1);
  x2 = roi_pts(ll+1,1);
  y1 = roi_pts(ll,2);
  y2 = roi_pts(ll+1,2);

  hold(in_fig,'on');
  plot(in_fig,[x1, x2], [y1, y2],'w','LineWidth',3);
end
    