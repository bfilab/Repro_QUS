function [alpha,kappa,s,omega] = computeHomodynedK(env)
%COMPUTEHOMEDYNEDK Compute Homodyned K fit and parameters for envelope,
%using the XU estimator.

% Generally following notation in Destrempes 2013, the parameters are:
% alpha = scatterer clustering parameter
% kappa = structure parameter
% mu = mean intensity
% hk = ratio of coherent to diffuse signal
% epsilon = coherent signal
% sigma = related to incoherent signal
% omega = related to sigma

% References:
% [1] Destrempes, F., et al. (2013). "Estimation Method of the Homodyned
% K-Distribution Based on the Mean Intensity and Two Log-Moments." SIAM J
% Imaging Sci 6(3): 1499-1530.
% [2} Cristea, A., et al. (2020). "Quantitative assessment of media
% concentration using the Homodyned K distribution." Ultrasonics 101:
% 105986.

% INPUTS:
%   env: envelope of RF data
% OUTPUTS:
%   alpha: scatterer clustering parameter
%   kappa: structure parameter

% 11/19/2020 (THL): Added s and omega outputs for plotting fit.

%% Compute HK fit with XU Estimator
[gamma, alpha, epsilon, sigma] = HK_estim_XU(env(:));           

%% Compute parameters
kappa = epsilon^2/(2*sigma^2*alpha);                               
mu = mean(env(:).^2);                                           
hk = epsilon/(sigma*sqrt(alpha));                                 
omega = mu/(hk^2+2);                                       
s = hk*sqrt(omega);                                            

end

