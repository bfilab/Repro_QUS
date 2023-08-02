function [bsc_faran,a_star] = fitFaranSphere(bsc,f,arange,param)
%FITFARANSPHERE Fit the Faran sphere form factor BSC.
% INPUTS:
%   bsc = backscatter coefficient [1/str/m]
%   f = frequency vector [Hz]
%   arange = range of radii for iterating fit [m]
%   param. = selected values for other parameters (kept constant)
%         c1 = speed of sound in scatterer [m/s]
%         c0 = speed of sound (background) [m/s]
%         n = average particle number density [m/s]
%         gamma = relative acoustic impedance between scatterer and surrounding
%                 medium
%         poisson = poisson for scatterer
%         pwire = scatterer density [g/cm^2]
%         ph2o = medium density [g/cm^2]
%         Nume = number iterations

% Iterate over radius
mse = nan(length(arange),1);
for i = 1:length(arange)
   
    % Compute Faran Sphere BSC
    ph = param;
    ph.a = arange(i);
    bsc_faran = computeBSC_mono(f,ph,ph.Nume);
    
    % Get Error
    mse(i) = mean((bsc - bsc_faran).^2);
    
    
    % Plot
    figure(1)
    subplot(2,1,1)
    plot(f*1e-6,10*log10(bsc),'LineWidth',2); hold on;
    plot(f*1e-6,10*log10(bsc_faran),'LineWidth',2); hold off;
    subplot(2,1,2)
    plot(arange*1e6,mse,'LineWidth',2);
    xlim([arange(1)*1e6,arange(end)*1e6]);
    
    pause;
    
    
end

% Get minimum error value
[~,i] = min(mse);
a_star = arange(i);

% Get fit BSC
ph.a = a_star;
bsc_faran = computeBSC_mono(f,ph,ph.Nume);

end

