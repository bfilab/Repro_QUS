function [J] = sph_bessel_1(x,n)
%SPH_BESSEL_1 Spherical Bessel function of the first kind.
% INPUTS:
%   x = input data on which to compute spherical Bessel function
%   n = order of Bessel function
% OUTPUTS:
%   J = output of spherical Bessel  function

% Reference: https://mathworld.wolfram.com/SphericalBesselFunctionoftheFirstKind.html

% 09/16/2020 (THL): Created

J = sqrt( pi./(2*x) ) .* besselj( n+0.5 , x);

end

