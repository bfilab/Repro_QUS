function [nak_param,scale_factor,nak_fit] = computeNakagami(env)
%COMPUTENAKAGAMI Compute Nakagami fit and parameters for envelope.
% INPUTS:
%   env = envelope data
% OUTPUTS:
%   nak_param = nakagami parameter
%   scale_factor = scaling factor
%   nak_fit = fitted distribution object

% 06/22/2020 (THL): Created
% 11/19/2020 (THL): Added Nakagami object as output

% Compute Nakagami fit
nak_fit = fitdist(env(:),'nakagami');                                      
mu = nak_fit.Params(1);
alpha = nak_fit.Params(2);

% Store results
nak_param = mu;
scale_factor = alpha;

end

