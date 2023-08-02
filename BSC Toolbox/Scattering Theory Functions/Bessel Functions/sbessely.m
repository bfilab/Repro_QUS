function [J] = sbessely(x,n)
%SBESSELY Spherical Bessel function of the second kind.
% https://mathworld.wolfram.com/SphericalBesselFunctionoftheSecondKind.html
% INPUTS:
%   x = data; column vector
%   n = order; row vector
% If size(x,2)>1, then sbessel will calculate the Bessel function for each
% x and each n.

six = size(x);
sin = size(n);

if any(six ~= sin)

    n  = n(:)';
    
    sin = size(n);
    nn = repmat( n,  [ six(1), 1] );
    
    
    J = zeros( six(1), sin(2), six(2) );
    
    for xi = 1:six(2)
        
        xx     = repmat( x(:, xi), [ 1, sin(2) ] );
        sqrtXi = sqrt( pi./(2 * xx) );
        
        J(:, :, xi) = sqrtXi .* bessely( nn+0.5, xx );
        
    end
    
else
    J = sqrt( pi./(2*x) ) .* bessely( n+0.5, x);
end

end

