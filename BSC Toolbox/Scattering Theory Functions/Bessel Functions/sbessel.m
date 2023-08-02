function [J] = sbessel(x,n)
% SBESSEL Spherical Bessel function of the first kind, order n.
% http://mathworld.wolfram.com/SphericalBesselFunctionoftheFirstKind.html

% Compute for a single order.
if length(n) == 1
    switch n
        case 0                                                              % Order 0
            J = sin(x) ./ x;
        case 1
            J = sin(x) ./ x.^2 - cos(x)./x;                                 % Order 1 (simplifies to negative derivative of sinc)
        case 2
            J = (3./x.^3 - 1./x) .* sin(x) - 3./x.^2 .* cos(x);             % Order 2
        otherwise
            J = sqrt( pi./(2*x) ) .* besselj( n+0.5 , x);
    end
    J( x==0 ) = 1;
    
% If you have multiple orders and multidimensional data.
% Compute the spherical Bessel function for all data and orders.
else
    
    six  = size(x);                                                         % Size of data
    sisn = size(n);                                                         % Number of orders
    
    % If x and n are not the same size:
    if any(six ~= sisn)
        
        n  = n(:)';                                                         % Get orders as a 1D vector
        
        sisn = size(n);                                                     % Number of orders
        nn   = repmat( n, [ six(1), 1] );                                   % Each column: one order to be applied to all data
        
        
        J = zeros( six(1), sisn(2), six(2) );
        
        for xi = 1:six(2)
            
            xx     = repmat( x(:, xi), [ 1, sisn(2) ] );                    % Each column: the data
            sqrtXi = sqrt( pi./(2 * xx) );
            
            J(:, :, xi) = sqrtXi .* besselj( nn+0.5, xx );                  % Each column: Bessel function of one order
            
        end
        
    % If n and x are of the same size, do the standard Bessel function.
    % This computes the n(i,j) order Bessel function for each x(i,j).
    else
        J = sqrt( pi./(2*x) ) .* besselj( n+0.5 , x);
        %  J = sqrt( pi./(2*x) ) .* bessely( n+0.5, x);
    end
    
end

end

