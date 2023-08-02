%% Create Objects to Hold Attenuation Parameters
% 05/05/2020 (THL): Created
% 06/29/2020 (THL): Added cartilage
% 01/27/2021 (THL): Added lung

%% Set parameters

% 15 um phantom
% Reference: 
% \\rri-usa.org\labs\BioMed\Data\Technical Documents\Phantoms\15 mu and 18 mu TPX\Attenuation\Attenuatoin2.xlsx
% \\rri-usa.org\labs\BioMed\Data\Technical Documents\Phantoms\MatlabRefObjects 
n = 1.842;
attn_params_15.a = 0.0254 * 1e2 * (1e-6)^n;                                 % [dB / m Hz^n]
attn_params_15.n = n;
attn_params_15.c = 0 * 1e2;                                                 % [dB / m ]

% 18 um phantom
% Reference: 
% \\rri-usa.org\labs\BioMed\Data\Technical Documents\Phantoms\15 mu and 18 mu TPX\Attenuation\Attenuatoin2.xlsx
% \\rri-usa.org\labs\BioMed\Data\Technical Documents\Phantoms\MatlabRefObjects 
n = 1.862;
attn_params_18.a = 0.0279 * 1e2 * (1e-6)^n;                                 % [dB / m Hz^n]
attn_params_18.n = n;
attn_params_18.c = 0 * 1e2;                                                 % [dB / m] 

% 60 um phantom
% Reference:
% \\rri-usa.org\labs\BioMed\Data\Technical Documents\Phantoms\MatlabRefObjects 
n = 1;
attn_params_60.a = 0.55 * 1e2 * (1e-6)^n;                                   % [dB / m Hz^n]
attn_params_60.n = n;
attn_params_60.c = -1 * 1e2;                                                % [dB / m]   

% prostate
% Reference: Rohrbach, D., et al. (2018). "High-Frequency Quantitative Ultrasound for Imaging Prostate Cancer Using a Novel Micro-Ultrasound Scanner." Ultrasound in Medicine & Biology 44(7): 1341-1354.
n = 1;
attn_params_pros.a = 0.5 * 1e2 * (1e-6)^n;                                  % [dB / m Hz^n]
attn_params_pros.n = n;
attn_params_pros.c = 0 * 1e2;

% cartilage
% Reference: Rohrbach, D., et al. (2017). "Regular chondrocyte spacing is a potential cause for coherent ultrasound backscatter in human articular cartilage." J Acoust Soc Am 141(5): 3105.
% Using the extracellular matrix value.
n = 1;
attn_params_cart.a = (10.65/45) * 1e3 * (1e-6)^n;                           % [dB / m Hz^n]
attn_params_cart.n = n;
attn_params_cart.c = 0;

% lung
n = 1;
attn_params_lung.a = 0.5 * 1e2 * (1e-6)^n;                                  % [dB / m Hz^n]
attn_params_lung.n = n;
attn_params_lung.c = 0;

% Description
description{1,1} = 'a*f^n + c';
description{2,1} = 'a';
description{2,2} = '[dB / m Hz^n]';
description{3,1} = 'n';
description{3,2} = 'exponent for frequency';
description{4,1} = 'c';
description{4,2} = '[dB / m]';

%% Save
save('attn_params_15um','attn_params_15','description');
save('attn_params_18um','attn_params_18','description');
save('attn_params_60um','attn_params_60','description');
save('attn_params_prostate','attn_params_pros','description');
save('attn_params_cartilage','attn_params_cart','description');
save('attn_params_lung','attn_params_lung','description');