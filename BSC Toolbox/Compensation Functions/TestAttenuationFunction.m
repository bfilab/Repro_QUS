%% Inputs

% Frequency
f = (10:0.1:30)*1e6;                                                        % [Hz]

% Attenuation Parameters
% a*f^n + c
a = 0.0254;                                                                 % [dB / cm MHz^n]
n = 1.842;
c = 0;                                                                      % [dB / cm]

attn_params.a = a * 1e2 * (1e-6)^n;                                         % [dB / m Hz^n]
attn_params.n = n;
attn_params.c = c * 1e2;                                                    % [dB / m]

% Attenuation distance
d_cm = 0.086624996;
d_m = d_cm*1e-2;

%% Calculate

% Prior to unit conversions (dB)
A1 = 4*(a.*(f*1e-6).^n + c)*(d_cm);

% After unit conversions (dB)
[A2, attn_lin] = compAttn(f,attn_params,d_m);

% Compute in linear scale, then convert to dB
A3 = 10*log10(attn_lin);

plot(A1,'LineWidth',3); hold on;
plot(A2,'r--','LineWidth',3);
plot(A3,'g.','MarkerSize',10); hold off;
axis tight;
