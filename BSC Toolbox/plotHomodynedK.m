function [hk_x,hk_y] = plotHomodynedK(alpha,s,omega,hi)
%PLOTHOMODYNEDK Return points to plot the homodyned K pdf for given
%parameters.
% INPUTS:

% 11/15/2021 (THL): Created

hk_x = hi.BinEdges(1:(end-1)) + hi.BinWidth/2;
hk_y = homok_func(hk_x, alpha, s, omega)*hi.BinWidth;

end

