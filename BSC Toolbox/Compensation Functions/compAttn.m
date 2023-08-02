function [attn_db, attn_lin] = compAttn(f,attn_params,d)
%COMPATTN Generate attenuation compensation.
% The attenuation compensation is calculated with a*f^n + b [dB].
% To compensate a power spectrum in dB: PS_db + attn_db
% To compensate a power spectrum in linear scale: PS_lin*attn_lin

% INPUTS:
% f = frequency [Hz]
% attn_params.
%           a = [dB/Hz^n/m]
%           n
%           c = [dB/m]
% d = distance over which attenuation occurs [m]

% OUTPUTS
% attn_db = attenuation in dB
% attn_lin = attenuation in linear scale

%% Unpack Inputs
a = attn_params.a;
n = attn_params.n;
c = attn_params.c;

%% Compute Attenuation in dB
attn_db = 4 * (a .* f.^n + c) * d;

%% Compute Attenuation in Linear Scale
attn_lin = exp(4*(log(10)/40)*attn_db);

end

