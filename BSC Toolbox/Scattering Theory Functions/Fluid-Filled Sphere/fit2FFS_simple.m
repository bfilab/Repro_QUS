function [bsc_ffs,param_ffs] = fit2FFS_simple(bsc,f,arange)
%FIT2FFS Fit a BSC to a theoretical BSC using a fluid-filled sphere form
%factor.
% INPUTS:
%   bsc = backscatter coefficient to be fit [1/str/m], NOT in dB
%         THE BSC SHOULD ALREADY BE CROPPED TO ITS BANDWIDTH.
%   f = frequency vector for the input bsc [Hz]
%   arange = range in which to test scatterer radius [m]
% OUTPUTS:
%   bsc_ffs = fitted bsc with fluid-filled sphere form factor [1/str/m]
%   params_ffs.
%              a = scatterer radius [m]
%              ac = acoustic concentration
% 
% v = VideoWriter('test.avi');
% v.FrameRate = 15;
% open(v);

%% Convert to dB
bsc_dB = 10*log10(bsc);

%% Remove frequency dependence
bsc2_dB = bsc_dB - 40*log10(f);

%% Fit the form factor by iterating over scatterer radius

% Preallocate variables
mse = nan(length(arange),1);                                                % mean squared error
all_shifts = nan(length(arange),1);                                         % shift to normalize magnitude
all_bsc_ffs = nan(length(arange),length(f));                                % fitted BSC (linear scale)

% Iterate over radius
for i = 1:length(arange)
   
    % Select radius
    a = arange(i);
    
    % Setup fluid filled sphere parameters
    param.a = a;
    param.c0 = 1500; %[m/s]
    param.n = 1;     % dummy parameter
    param.gamma = 1; % dummy parameter
    
    % Compute fluid-filled sphere BSC
    [bsc_ffs] = computeBSC_mono_FFS(f,param);
    
    % Convert to dB and subtract frequency dependence
    bsc_ffs_dB = 10*log10(bsc_ffs);
    bsc2_ffs_dB = bsc_ffs_dB - 40*log10(f);
    
    % Normalize
    shift_factor = mean(bsc2_dB) - mean(bsc2_ffs_dB);
    bsc2_ffs_dB_norm = bsc2_ffs_dB + shift_factor;
    
    % Compute error
    mse(i) = mean((bsc2_ffs_dB_norm - bsc2_dB).^2);
    
    % Store shift and bsc
    all_shifts(i) = shift_factor;
    all_bsc_ffs(i,:) = 10.^((bsc2_ffs_dB_norm + 40*log10(f))./10);
        
    % Plot

 %{
    figure(1);
    subplot(1,3,1)
    plot(f*1e-6, bsc_dB - mean(bsc_dB),'LineWidth',2); hold on;
    plot(f*1e-6,bsc_ffs_dB - mean(bsc_ffs_dB),'r','LineWidth',2); hold off;
    axis tight;
    xlabel('MHz'); ylabel('dB/str/m');
    title('Theoretical Fluid-Filled Sphere BSC');
    
    subplot(1,3,2)
    plot(f*1e-6,bsc2_dB,'LineWidth',2); hold on
    plot(f*1e-6,bsc2_ffs_dB_norm,'r','LineWidth',2); hold off;
    xlabel('MHz'); ylabel('dB/str/m');
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
    %}
   
 
    
end

% close(v);

%% Find optimal a* and other parameters

% Radius
[~,i] = min(mse);
a_star = arange(i);
param_ffs.a = a_star;

% Shift
shift_star = all_shifts(i);

% BSC (linear scale, not log-f subtracted)
bsc_ffs = all_bsc_ffs(i,:);

%% Find optimal acoustic concentration

% Calculate fluid filled sphere FF
F_dB = 10*log10(calcFluidFilledSphereFF(f,a_star,param.c0));

% Subtract form factor and log frequency to get term for solving ac
alpha_dB = 10*log10(bsc_ffs) - F_dB - 40*log10(f);
alpha = 10.^(mean(alpha_dB)/10);

% Solve for ac
Vs = (4/3)*pi*a_star^3;
den = (2*pi/param.c0)^4 * Vs^2;
ac = alpha*(4*pi^2)/den;
param_ffs.ac = ac;

end

