function [contrast,CNR,gCNR] = computeContrastMetrics(data1,data2,edges)
%COMPUTECONTRAST Compute the contrast metrics between two different sets of
%data (group 1 and group 2).
% INPUTS:
%   data1: data points corresponding to group 1
%   data2: data points corresponding to group 2
%   edges: edges for histogram
% OUTPUTS:
%   contrast: 3x1 vector, (1) u1/u2, (2) = u2/u1, (3) = max contrast
%   CNR: contrast to noise ratio; ratio of means normalized by standard
%        deviation
%   gCNR: measure of overlap in PDFs

% 10/14/2020 (THL): Created
% 10/15/2020 (THL): Added outputs for contrast, now a 3-elem vector.

% Reference: Morris, D. C., D. Y. Chan, T. H. Lye, H. Chen, M. L. Palmeri,
% T. J. Polascik, W.-C. Foo, J. Huang, J. Mamou and K. R. Nightingale
% (2020). "Multiparametric Ultrasound for Targeting Prostate Cancer:
% Combining ARFI, SWEI, QUS and B-Mode." Ultrasound in Medicine & Biology.

%% Convert to row vector
data1 = data1(:);
data2 = data2(:);

%% set minimum to 0
min_val = min([data1;data2]);
data1 = data1 - min_val;
data2 = data2 - min_val;

%% Comptue Statistics
mu1 = mean(data1);
mu2 = mean(data2);
omega1 = std(data1);
omega2 = std(data2);

%% Contrast
contrast(1) = mu1/mu2;
contrast(2) = mu2/mu1;
contrast(3) = max([contrast(1),contrast(2)]);

%% CNR
CNR = abs(mu1-mu2)/sqrt(omega1^2 + omega2^2);

%% gCNR

% Get edges
if isempty(edges)
    numBins = 300;
    edges = linspace(min([data1;data2]),max([data1;data2]),numBins);
end

% Estimate PDFS
[N1,~] = histcounts(data1,edges,'Normalization', 'Probability');
[N2,~] = histcounts(data2,edges,'Normalization', 'Probability');

% Take minimum and compute gCNR
ovl = sum(min([N1; N2]));
gCNR = 1-ovl;

end

