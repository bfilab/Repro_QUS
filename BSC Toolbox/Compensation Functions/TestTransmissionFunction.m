%% Test Transmission Function

% Frequency
f = (0:0.01:70)*1e6;                                                        % [Hz]

%% Attenuation Parameters

% Attenuation Parameters (Saran)
%{
a = 5;                                                                      % [Nepers / m MHz^n]
n = 1.5;
attn_params.a = a * (1e-6)^n;                                               % [Nepers / m Hz^n]
attn_params.n = n;
%}

% Attenuation Parameters (TPX)

% These TPX values are obtained from Madsen, E. L., et al. (2011).
% "Properties of phantom tissuelike polymethylpentene in the frequency
% range 20-70 MHZ." Ultrasound Med Biol 37(8): 1327-1339.

% Note that in the manuscript, there are TWO sets of values obtained for
% TPX attenuation properties because they  tried two different estimation
% methods to get the TPX attenuation. The values here correspond to METHOD
% 2. Based on the results of the paper, the values of the two methods
% should be roughly equivalent, so it shouldn't matter too much. The
% equivalency of the two TPX attenuation values is demonstrated in Fig. 6
% of the manuscript.

a = 0.906;                                                                  % [dB / cm MHz^n]
n = 1.324;
attn_params.a = a * log(10)/20 * 1e2 * (1e-6)^n ;                           % [Nepers / m Hz^n]
attn_params.n = n;

%% Thickness of membrane
%{
% d_um = 25;                                                                % [um] saran
%}
% 128 um is the thickness of TPX used in Madsen 2011 for method 2
% (not sure about the thickness of TPX from our lab!)
d_um = 128;                                                                 % [um] TPX
d_m = d_um*1e-6;

% Transmission Parameters (Saran)
%{
% transm_params.Z1 = 1.5;                                                     % Z1 = acoustic impedance of first material [MRayl]
% transm_params.Z2 = 1.5;                                                     % Z2 = acoustic impedance of second material [MRayl]
% transm_params.Zs = 4.064;                                                   % Zs = acoustic impedance of surface membrane [MRayl]
% transm_params.c = 2400;                                                     % speed of sound in membrane [m/s]
%}

% Transmission Parameters (TPX)
% Again, these parameters are obtained from Madsen 2011. The acoustic
% impedances Z1 and Z2 are for water. The acoustic impedance Zs is for TPX.
% The speed of sound 2093 m/s corresponds to the speed of sound of TPX at
% 40 MHz. The speed of sound is frequency dependent with a curve as shown
% in Fig. 4 of Madsen 2011.
transm_params.Z1 = 1.5;                                                     % Z1 = acoustic impedance of first material [MRayl]
transm_params.Z2 = 1.5;                                                     % Z2 = acoustic impedance of second material [MRayl]
transm_params.Zs = 1.82;                                                    % Zs = acoustic impedance of surface membrane [MRayl]
transm_params.c = 2093;                                                     % speed of sound in membrane [m/s]

%% Compute Transmission

[T] = compTransm(f,transm_params,attn_params,d_m);

% This plot is equivalent to Figure 5 of Madsen 2011, which plots the
% transmission coefficient of TPX with a water-TPX-water interface (given
% in linear scale, not dB).
plot( f*1e-6 , abs(T{1}) );
xlabel('MHz'); ylabel('dB'); axis tight;
title('TPX Transmission Coefficient');

% With this, you can compensate for transmission in dB like so (note that I
% made the square negative!!! -which is different from what I sent you
% previously).

% samp.transm_db = 10*log10( ( abs(T{1}) .* abs(T{2}) ).^(-2) );
% samp.avg_PS_db_transmComp = samp.avg_PS_db + samp.transm_db;