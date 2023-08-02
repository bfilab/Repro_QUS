function [time] = dist2time(c,depth)
%DIST2TIME 
% INPUTS:
% c = speed of sound
% depth = depth

time = (depth/c)*2;

end

