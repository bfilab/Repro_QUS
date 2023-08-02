function [T] = compTransm(f,transm_params,attn_params,d)
%COMPTRANSM Compensate transmission through the phantom membrane (both
% directions through the membrane). To compensate a phantom power spectrum,
% do: PS_phantom * (abs(T{1}*abs(T{2}))^2.

% Reference: Wear, K. A., et al. (2005). "Interlaboratory comparison of
% ultrasonic backscatter coefficient measurements from 2 to 9 MHz." J
% Ultrasound Med 24(9): 1235-1250.

% INPUTS
%   f = frequency [Hz]
% transm_params.
%          Z1 = acoustic impedance of first material [MRayl]
%          Z2 = acoustic impedance of second material [MRayl]
%          Zs = acoustic impedance of surface membrane [MRayl]
%          c = speed of sound in membrane [m/s]
% attn_params. = attenuation parameters in membrane (a*f^n)
%           a = [Nepers/Hz^n/m]
%           n
%   d = thickness of membrane [m]

%% Unpack Inputs

Z1 = transm_params.Z1;
Z2 = transm_params.Z2;
Zs = transm_params.Zs;
c = transm_params.c;

%% Compute Parameters

% Calculate attenuation
alf = attn_params.a.*f.^(attn_params.n);                                    % [nepers / meter]

% Compute k
k   = 2.*pi.*f./ c - 1i .* alf;

% Multiply by thickness
kl  = k.*d;

%% Compute transmission coefficient
T = cell(2,1);
for i = 1:2
    
    if i == 1
        Zi = Z1;
        Zf = Z2;
    else
        Zi = Z2;
        Zf = Z1;
    end
    
    % Compute transmission coefficient
    N = (Zi+Zf).*cos(kl) + 1i.*( Zs + Zi*Zf/Zs ).*sin(kl);
    T{i} = 2*Zf ./ N;
    
end

end


