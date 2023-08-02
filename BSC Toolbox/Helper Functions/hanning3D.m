function [V_gate] = hanning3D(V)
%Hanning3D Gate 3D RF data with a Hanning window. Assumes dimension 1 =
%axial dimension.
% INPUTS:
%   V: 3D RF data, where axial dimension = dimension 1
% OUTPUTS:
%   V_gate: Hanning windowed 3D RF data

% 4/29/2020 (THL): Created

% Create 3D Hanning window
han_win = hanning(size(V,1));
han_win = repmat(han_win,[1,size(V,2)]);
han_win = repmat(han_win,[1,1,size(V,3)]);

% Apply window
V_gate = V.*han_win;

end

