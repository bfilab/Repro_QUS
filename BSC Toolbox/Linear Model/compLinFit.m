function [spectral_slope, intercept, midband_fit] = compLinFit(f,bsc_db,bw,cf)
%COMPLINFIT Compute linear fit and extract QUS parameters.
% INPUTS:
%   f = frequency for BSC [Hz]
%   bsc_db = backscatter coefficient [dB]
%   bw = bandwidth, 2-element vector [Hz]
%   cf = center frequency, one scalar [Hz]
% OUTPUTS:
%   spectral_slope = spectral slope [dB/MHz]
%   intercept = intercept [dB]
%   midband_fit = midband fit at center frequency [dB]

% 6/24/2020 (THL): Created

%% Ensure vectors are same size
f = reshape(f,[length(f),1]);
bsc_db = reshape(bsc_db,[length(bsc_db),1]);

%% Get frequency locations
[~,i1] = min( abs(f-bw(1)) );
[~,i2] = min( abs(f-bw(2)) );

%% Compute linear fit to BSC
p = polyfit(f(i1:i2)*1e-6, bsc_db(i1:i2),1);                                % Linear fit

%% Get parameters
spectral_slope = p(1);                                                      % [dB/MHz]
intercept = p(2);                                                           % [dB]
midband_fit = spectral_slope*(cf*1e-6) + intercept;                         % [dB]

end

% sub_f = f(i1:i2)/1e6;
% sub_bsc = bsc_db(i1:i2);
% fit_y = sub_f.*p(1) + p(2);
% figure(759); clf;
% plot(sub_f,sub_bsc,'k','LineWidth',2);
% hold on;
% plot(sub_f,fit_y,':r','Linewidth',2);
% xlim([55, 105]);
% xlabel('Frequency (MHz)');
% ylabel('W_{ROI}(f) (dB)');
% legend('Measured','Linear Fit','Location','best');
% set(gca,'FontSize',14);
% saveas(gcf,'example_linear_fit_3','fig');
% saveas(gcf,'example_linear_fit_3','png');


