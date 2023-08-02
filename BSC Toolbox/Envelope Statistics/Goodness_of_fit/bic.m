function bic = bic(PD)
%BIC Bayesian information criterion for evaluating the goodness of fit of a
% distribution.
% INPUTS:
%   PD = distribution object (e.g. NakagamiDistribution)
% OUTPUTS:
%   bic = Bayesian information criterion

NLL = PD.NLogL;
n = numel(PD.InputData.data);
k = PD.NumParameters;
bic =-2*(-NLL)+k*log(n); 

end

