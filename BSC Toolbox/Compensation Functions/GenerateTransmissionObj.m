%% Create Objects to Hold Transmission Parameters
% Created 5/5/2020 (THL)
% Modified 6/24/2020 (THL): Combined the two output structs into one.

%% Set Parameters

% Saran
% Reference: Wear, K. A., et al. (2005). "Interlaboratory comparison of ultrasonic backscatter coefficient measurements from 2 to 9 MHz." J Ultrasound Med 24(9): 1235-1250.
transm_params_saran.Z1 = 1.5;                                               % Z1 = acoustic impedance of first material [MRayl]
transm_params_saran.Z2 = 1.5;                                               % Z2 = acoustic impedance of second material [MRayl]
transm_params_saran.Zs = 4.0560;                                            % Zs = acoustic impedance of surface membrane [MRayl]
transm_params_saran.c = 2400;                                               % speed of sound in membrane [m/s]

n = 1.5;
attn_params_saran.a = 5 * (1e-6)^n;                                         % [Nepers / m Hz^n]
attn_params_saran.n = n;
transm_params_saran.attn = attn_params_saran;

% TPX
% Reference: Madsen, E. L., et al. (2011). "Properties of phantom tissuelike polymethylpentene in the frequency range 20-70 MHZ." Ultrasound Med Biol 37(8): 1327-1339.
transm_params_tpx.Z1 = 1.5;                                                 % Z1 = acoustic impedance of first material [MRayl]
transm_params_tpx.Z2 = 1.5;                                                 % Z2 = acoustic impedance of second material [MRayl]
transm_params_tpx.Zs = 1.82;                                                % Zs = acoustic impedance of surface membrane [MRayl]
transm_params_tpx.c = 2093;                                                 % speed of sound in membrane [m/s]

n = 1.324;
attn_params_tpx.a = 0.906 * log(10)/20 * 1e2 * (1e-6)^n;                    % [Nepers / m Hz^n]
attn_params_tpx.n = n;
transm_params_tpx.attn = attn_params_tpx;

% Description
description{1,1} = 'Z1';
description{1,2} = 'acoustic impedance of first material [MRayl]';
description{2,1} = 'Z2';
description{2,2} = 'acoustic impedance of second material [MRayl]';
description{3,1} = 'Zs';
description{3,2} = 'acoustic impedance of surface membrane [MRayl]';
description{4,1} = 'c';
description{4,2} = 'speed of sound in membrane [m/s]';
description{5,1} = 'attn.a';
description{5,2} = 'membrane attenuation [Nepers / m Hz^n]';
description{6,1} = 'attn.n';
description{6,2} = 'membrane attenuation, exponent for frequency';

%% Save

save('transm_attn_params_saran','transm_params_saran','description');
save('transm_attn_params_tpx','transm_params_tpx','description');