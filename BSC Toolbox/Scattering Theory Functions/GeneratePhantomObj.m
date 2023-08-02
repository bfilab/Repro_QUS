%% Create Objects to Hold Phantom Scattering Parameters
% Created 5/5/2020 (THL)
% Modified 6/24/2020 (THL): changed average particle number density to
%                           particles/m^3

% Parameters were referenced from:
% \\rri-usa.org\labs\BioMed\Data\Technical Documents\Phantoms\15 mu and 18 mu TPX\Bead Specs.xlsx

%% Set Parameters

% 15 um phantom parameters
ph_15.a = (15.9/2)*1e-6;                                                    % Scatterer radius [m]
ph_15.c1 = 5640;                                                            % speed of sound in scatterer [m/s]
ph_15.c0 = 1500;                                                            % speed of sound in medium [m/s]
ph_15.n = 2.926e12;                                                         % average particle number density [particles/m^3]
ph_15.gamma = 8.4;                                                          % relative acoustic impedance between scatterer and surrounding medium
ph_15.poisson = 0.2;                                                        % poisson
ph_15.pwire = 2.5;                                                          % scatterer density [g/cm^2]
ph_15.ph2o = 1;                                                             % medium density [g/cm^2]
ph_15.Nume = 20;

% 18 um phantom parameters
ph_18.a = (18.2/2)*1e-6;                                                    % Scatterer radius [m]
ph_18.c1 = 5640;                                                            % speed of sound in scatterer [m/s]
ph_18.c0 = 1500;                                                            % speed of sound in medium [m/s]
ph_18.n = 1.56e12;                                                          % average particle number density [particles/m^3]
ph_18.gamma = 8.4;                                                          % relative acoustic impedance between scatterer and surrounding medium
ph_18.poisson = 0.2;                                                        % poisson
ph_18.pwire = 2.5;                                                          % scatterer density [g/cm^2]
ph_18.ph2o = 1;                                                             % medium density [g/cm^2]
ph_18.Nume = 20;

% 60 um phantom parameters
ph_60.a = (60/2)*1e-6;                                                      % Scatterer radius [m]
ph_60.c1 = 5640;                                                            % speed of sound in scatterer [m/s]
ph_60.c0 = 1500;                                                            % speed of sound in medium [m/s]
ph_60.n = 4.681e10;                                                         % average particle number density [particles/m^3]
ph_60.gamma = 8.4;                                                          % relative acoustic impedance between scatterer and surrounding medium
ph_60.poisson = 0.2;                                                        % poisson
ph_60.pwire = 2.5;                                                          % scatterer density [g/cm^2]
ph_60.ph2o = 1;                                                             % medium density [g/cm^2]
ph_60.Nume = 20;

% Description
description{1,1} = 'a';
description{1,2} = 'Scatterer radius [m]';
description{2,1} = 'c1';
description{2,2} = 'speed of sound in scatterer [m/s]';
description{3,1} = 'c0';
description{3,2} = 'speed of sound in medium [m/s]';
description{4,1} = 'n';
description{4,2} = 'average particle number density [particles/m^3]';
description{5,1} = 'gamma';
description{5,2} = 'relative acoustic impedance between scatterer and surrounding medium';
description{6,1} = 'poisson';
description{6,2} = 'poisson';
description{7,1} = 'pwire';
description{7,2} = 'scatterer density [g/cm^3]';
description{8,1} = 'ph2o';
description{8,2} = 'medium density [g/cm^3]';
description{9,1} = 'Nume';
description{9,2} = 'number of iterations to compute scattering theory';

%% Save
save('ph_15','ph_15','description');
save('ph_18','ph_18','description');
save('ph_60','ph_60','description');