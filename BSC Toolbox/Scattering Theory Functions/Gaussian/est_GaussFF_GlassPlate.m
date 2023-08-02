function [ac,aeff,fit_curve,f_bw] = est_GaussFF_GlassPlate(f, bsc_db, bandwidth, c, roi_mid_z, roi_len_z,trans_type)
%UNTITLED Estimate acoustic concentration and effective scatterer radius
%with a Gaussian form factor fit.
% INPUTS:
%   f = frequency [Hz] - Not MHz, the conversion is done in this function.
%   bsc_db = backscatter coefficient [dB]
%   bandwidth = bandwidth of transducer [Hz]
%   c = speed of sound [m/s]
%   roi_mid_z = axial midpint of the ROI [m]
%   roi_len_z = axial length of ROI [m]
%   trans_type = string specifying transducer. Either '80' or 'GE'
% OUTPUTS:
%   ac = acoustic concentration [dB/m^3]
%   aeff = effective scatterer radius [meters] 
%   fit_curve = fitted Gaussian curve in bandwidth [dB]
%   f_bw = frequency vector corresponding to fit_curve [Hz]

% 10/26/2020 (THL): (1) returns fitted Gaussian curve

%% Reshape to column vectors
f = reshape(f,[],1);
bsc_db = reshape(bsc_db,[],1);

%% Get bandwidth indices
[~,f1] = min(abs(bandwidth(1)-f));
[~,f2] = min(abs(bandwidth(2)-f));

%% Initial computations

% Subtract log frequency to get Gaussian scatter model
gs_db = bsc_db - 10*log10(f.^4);                                            % dB/m/sr/MHz^4

% Get Gaussian scatter model in bandwidth
gs_db_bw = gs_db(f1:f2);
f_bw = f(f1:f2);

%% Get fit
[p,~,mu] = polyfit( (f_bw*1e-6).^2, gs_db_bw, 1 );                          % Fit in MHz
p(2) = p(2) - p(1)*mu(1)/mu(2);                                             % For computational purposes, fitting was done with scaled data (see help polyfit), so scale the coefficients back
p(1) = p(1)/mu(2);

% p(1) = slope, p(2) = intercept
M = p(1);                                                                   % dB/MHz^2
I = p(2);                                                          

fit_curve = polyval(p,(f_bw*1e-6).^2);
fit_curve = fit_curve + 10*log10(f_bw.^4);

%% Compute effective scatterer radius (aeff)
num = log( 10^(M/10) );                                                    
den = -0.827 * (2*pi/c)^2;                                                 
aeff = sqrt(num/den);                                                       % [um]
aeff = aeff*1e-6;                                                           % [m]

alt_aeff = real(2*sqrt(-p(1)/(4.34*(12.159+2.66*(3.18/roi_mid_z).^2))));

aeff = alt_aeff/2;
%% Compute Acoustic Concentration (ac)
% num = 10^(I/10);
den = (4/9)*(2*pi/c)^4*aeff^6;
% ac = num/den;
ac = I - 10*log10(den);                                                     % [dB]

alt_ac = p(2)-10*log10(185*roi_len_z*(0.5e-3/roi_mid_z).^2*(aeff)^6);

switch trans_type
    case 'GE'
        % Assume a constant F-number of 2 for GE data, so the q = (1/2F#) = 0.25
        alt_ac = p(2)-10*log10(185*roi_len_z*(0.25).^2*(aeff)^6);

    case '80'
        % Radius of 80MHz active element is about 1mm
        alt_ac = p(2)-10*log10(185*roi_len_z*(1e-3/roi_mid_z).^2*(aeff)^6);
end

ac = alt_ac;
