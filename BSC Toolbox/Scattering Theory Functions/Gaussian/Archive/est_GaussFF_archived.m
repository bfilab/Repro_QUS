function [ac, aeff] = est_GaussFF(f, bsc_db, bandwidth, c)
%UNTITLED Estimate acoustic concentration and effective scatterer radius
%with a Gaussian form factor fit.
% INPUTS:
%   f = frequency [Hz] - Not MHz, the conversion is done in this function.
%   bsc_db = backscatter coefficient [dB]
%   bandwidth = bandwidth of transducer [Hz]
%   c = speed of sound [m/s]
% OUTPUTS:
%   ac = acoustic concentration [1/m^3]
%   aeff = effective scatterer radius [m]

%% Reshape to column vectors
f = reshape(f,[],1);
bsc_db = reshape(bsc_db,[],1);

%% Get bandwidth indices
[~,f1] = min(abs(bandwidth(1)-f));
[~,f2] = min(abs(bandwidth(2)-f));

%% Initial computations

% Subtract log frequency to get Gaussian scatter model
gs_db = bsc_db - 10*log10(f.^4);

% Get Gaussian scatter model in bandwidth
gs_db_bw = gs_db(f1:f2);
f_bw = f(f1:f2);

%% Get fit
[p,~,mu] = polyfit( (f_bw*1e-6).^2, gs_db_bw, 1 );                          % Fit in MHz
p(2) = p(2) - p(1)*mu(1)/mu(2);                                             % For computational purposes, fitting was done with scaled data (see help polyfit), so scale the coefficients back
p(1) = p(1)/mu(2);

% p(1) = slope, p(2) = intercept
M = p(1);                                                                   % dB/MHz
I = p(2);                                                                   % dB

% fit_curve = polyval(p,(f_bw*1e-6).^2);

%% Compute effective scatterer radius (aeff)
num = log( 10^(M/10) );
den = -0.827 * (2*pi/c)^2;
aeff = sqrt(num/den);                                                       % [um]
aeff = aeff*1e-6;                                                           % [m]

%% Compute Acoustic Concentration (ac)
% num = 10^(I/10);
den = (4/9)*(2*pi/c)^4*aeff^6;
% ac = num/den;
ac = I - 10*log10(den);                                                     % [dB]

end

