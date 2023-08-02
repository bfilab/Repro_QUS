function [hk_kl] = computeKLDiv_HK(env,model)
%COMPUTEKLDIV_HK Compute KL divergence for Homodyned-K fit.
% INPUTS:
%   env = envelope values
%   model
%        .alpha = scatterer clustering parameter
%        .s = (ratio of coherent to diffuse signal)*sqrt(omega);   
%        .omega


% 03/01/2021 (THL): Created

%% Get distribution of envelope

[hi.Values,hi.BinEdges] = histcounts(env(:),...
                                     100,...
                                     'Normalization','Probability');
hi.BinWidth = mean(diff(hi.BinEdges));

%% Compute KL divergence of HK fit

alpha = model.hk.alpha;
s     = model.hk.s;
omega = model.hk.omega;

hk_x = hi.BinEdges(1:(end-1)) + hi.BinWidth/2;
hk_y = homok_func(hk_x, alpha, s, omega)*hi.BinWidth;
hk_kl = KLdiv(hi.Values,hk_y);

end

