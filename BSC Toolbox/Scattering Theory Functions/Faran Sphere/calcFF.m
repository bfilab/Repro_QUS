function F = calcFF(f,poisson,a,c,c1,pwire,ph2o,Nume)
%CALCFF Calculate Faran sphere form factor (solid sphere). Derived from QUS
%processing source code.

% INPUTS:
%   f = frequency vector [Hz]
%   poisson = poisson for scatterer
%   a = scatterer radius (aka aeff)
%   c = speed of sound in surrounding medium [m/s]
%   c1 = speed of sound in scatterer [m/s]
%   pwire = scatterer density [g/cm^2]
%   ph2o = medium density [g/cm^2]
%   Nume = number of iterations

% OUTPUTS:
%   F = form factor function

% Reference: Jr., J. J. F. (1951). "Sound Scattering by Solid Cylinders and
% Spheres." The Journal of the Acoustical Society of America 23(4):
% 405-418.

%% Compute Wavenumber

% Calculate c2 (velocity of shear wave in the scatterer)
c2 = sqrt((.5*c1^2-poisson*c1^2)/(1-poisson));                              % [m/s]

% Calculate Wavenumber * Scatterer Radius (using c for surrounding fluid)
% k = w/c, where c = velocity of sound in fluid surrounding scatterer
k = (2*pi*f)/c;                                                             % [1/m]
ka = k*a;

% Calculate Wavenumber * Scatterer Radius (using c for scatterer)
% k1 = w/c1, where c1 = velocity compressional waves in scatterer
% k2 = w/c2, where c2 = velocity of shear waves in the scatterer
% x1 = k1*a;
% x2 = k2*a;
x1 = ka.*c/c1;
x2 = ka.*c/c2;

%% Calculate Spherical Bessel Functions of the First and Second Kind

% Calculate Sbessel min an plus one
in = -1:Nume+1;

% Calculate spherical Bessel of the first and second kind for ka
% (surrounding fluid)
J = sbessel(ka, in);
N = sbessely(ka, in);

% Calculate spherical Bessel of the first and second kind for ka
% (1. compressional waves in scatterer, 2. shear waves in the scatterer)
Jx1 = sbessel(x1, in);
Jx2 = sbessel(x2, in);

%%

xJp  = -0.5*J(:,2:end-1)   + bsxfun( @times, 0.5*ka, (J(:,1:end-2)   -   J(:, 3:end) ));
xNp  = -0.5*N(:,2:end-1)   + bsxfun( @times, 0.5*ka, (N(:,1:end-2)   -   N(:, 3:end) ));
x1Jp = -0.5*Jx1(:,2:end-1) + bsxfun( @times, 0.5*x1, (Jx1(:,1:end-2) - Jx1(:, 3:end) ));
x2Jp = -0.5*Jx2(:,2:end-1) + bsxfun( @times, 0.5*x2, (Jx2(:,1:end-2) - Jx2(:, 3:end) ));

del   = -J(:,2:end-1)./N(:,2:end-1);
alpha = -xJp./J(:,2:end-1);
beta  = -xNp./N(:,2:end-1);

alphax1 = -x1Jp./Jx1(:,2:end-1);
alphax2 = -x2Jp./Jx2(:,2:end-1);

iii   = repmat( in(2:end-1) .* in(2:end-1) + in(2:end-1), [ size(J, 1), 1 ] );

x22   = repmat(  x2.*x2/2, [ 1, size(alphax1, 2) ] );

Num   = alphax1./(alphax1+1)  -  iii./( alphax2+ iii -1 - x22 );
Denom = (2*alphax1 + iii - x22 )./( alphax1+1)-((iii) .* (alphax2+1))./(alphax2+iii-1-x22);

etc = (-x22).*Num./Denom;

phi = -etc .* ph2o/pwire;

eta_el = atan( del .* (alpha + phi)./(phi + beta));
eta    = atan((del .*  alpha)./beta );

expEta    = exp( -1i .* eta    );
expEta_el = exp( -1i .* eta_el );

ii   = ( 2.*in(2:end-1) + 1) .*(-1).^in(2:end-1);

sumEta    = nansum( bsxfun( @times, ii, sin( eta    ) ) .* expEta,    2 );
sumEtael  = nansum( bsxfun( @times, ii, sin( eta_el ) ) .* expEta_el, 2 );

sumEta   = 2.* abs(sumEta./(ka));
sumEtael = 2.* abs(sumEtael./(ka));

sume_el  = ((sumEtael)./ka.^2).^2;
sume_elN = sume_el/sume_el(2);

F = sume_elN;
F(1) = 1;

end
