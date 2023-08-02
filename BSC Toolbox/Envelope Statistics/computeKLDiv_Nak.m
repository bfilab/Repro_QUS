function [kl_nak] = computeKLDiv_Nak(env,nak_fit)
%COMPUTEKLDIV_NAK Compute KL divergence for fitted Nakagami model.
% INPUTS:
%   env = envelope values
%   nak_fit = fitdist(env(:),'nakagami')
% OUTPUTS:
%   kl_nak = kl divergence for fitted Nakagami

% 03/01/2021 (THL): Created

% Envelope statistics
[hi.Values,hi.BinEdges] = histcounts(env(:),...
                                     100,...
                                     'Normalization','Probability');
hi.BinWidth = mean(diff(hi.BinEdges));

% Get KL divergence for fitted Nakagami
nak_x = hi.BinEdges(1:(end-1))+hi.BinWidth/2;
nak_y = pdf(nak_fit, nak_x)*hi.BinWidth;
kl_nak = KLdiv(hi.Values,nak_y);

end

