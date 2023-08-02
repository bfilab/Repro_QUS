function [bsc_faran,param_faran] = fitFaranSphere_simple(bsc,f,arange,param)
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

% v = VideoWriter('test.avi');
% v.FrameRate = 15;
% open(v);

%% Convert to dB
bsc(bsc==0) = eps;
bsc_dB = 10*log10(bsc);

%% Remove frequency dependence
bsc2_dB = bsc_dB - 40*log10(f);

%% Fit the form factor by iterating over scatterer radius

% Preallocate variables
mse = nan(length(arange),1);                                                % mean squared error
all_shifts = nan(length(arange),1);                                         % shift to normalize magnitude
all_bsc_faran = nan(length(arange),length(f));                                % fitted BSC (linear scale)

% Iterate over radius
for i = 1:length(arange)
   
    % Select radius
    a = arange(i);
    
    % Compute Faran BSC
    param.a = a;
    bsc_faran = computeBSC_mono(f,param,param.Nume);
    bsc_faran(bsc_faran==0) = eps;
    
    % Convert to dB and subtract frequency dependence
    bsc_faran_dB = 10*log10(bsc_faran);
    bsc2_faran_dB = bsc_faran_dB - 40*log10(f);
    
    % Normalize
    shift_factor = mean(bsc2_dB) - mean(bsc2_faran_dB);
    bsc2_faran_dB_norm = bsc2_faran_dB + shift_factor;
    
    % Compute error
    mse(i) = mean((bsc2_faran_dB_norm - bsc2_dB).^2);
    
    % Store shift and bsc
    all_shifts(i) = shift_factor;
    all_bsc_faran(i,:) = 10.^((bsc2_faran_dB_norm + 40*log10(f))./10);
        
    % Plot
    %{

    figure(1);
    subplot(1,3,1)
    plot(f*1e-6,bsc_faran_dB,'r','LineWidth',2); axis tight;
    xlabel('MHz'); ylabel('dB');
    title('Theoretical Faran Sphere BSC');
    
    subplot(1,3,2)
    plot(f*1e-6,bsc2_dB,'LineWidth',2); hold on
    plot(f*1e-6,bsc2_faran_dB_norm,'r','LineWidth',2); hold off;
    xlabel('MHz'); ylabel('dB');
    title('Log-Freq Subtracted, Normalized, BSC');
    axis tight;
    
    subplot(1,3,3)
    plot(arange*1e6,mse,'LineWidth',2);
    title('Mean Squared Error');
    xlabel('Radius [um]');
    ylabel('Error');
    xlim([arange(1)*1e6,arange(end)*1e6]);

    sgtitle(['a = ', num2str(a*1e6), ' um']);
    
    set(gcf, 'Position', get(0, 'Screensize'));
    
    frame = getframe(gcf);
    writeVideo(v,frame);
    
    pause;
    %}
    
end

%close(v);

%% Find optimal a* and other parameters

% Radius
[~,i] = min(mse);
a_star = arange(i);
param_faran.a = a_star;

% Shift
shift_star = all_shifts(i);

% BSC (linear scale, not log-f subtracted)
bsc_faran = all_bsc_faran(i,:);

%% Find optimal acoustic concentration

% Calculate fluid filled sphere FF
poisson = param.poisson;
c1 = param.c1;
pwire = param.pwire;
ph2o = param.ph2o;
c0 = param.c0;
Nume = param.Nume;
F_dB = 10*log10(calcFF(f',poisson,a_star,c0,c1,pwire,ph2o,Nume));

% Subtract form factor and log frequency to get term for solving ac
alpha_dB = 10*log10(bsc_faran) - F_dB' - 40*log10(f);
alpha = 10.^(mean(alpha_dB)/10);

% Solve for ac
Vs = (4/3)*pi*a_star^3;
den = (2*pi/param.c0)^4 * Vs^2;
ac = alpha*(4*pi^2)/den;
param_faran.ac = ac;

end