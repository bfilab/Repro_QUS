function [map] = convertUnits(map,param_name)
%CONVERTUNITS Convert units from QUS results
% INPUTS:
%   map = QUS values
%   param_name = parameter name for the QUS values
% OUTPUTS:
%   map = map with converted units

% 02/25/2021 (THL): Created

if strcmp(param_name,'Nak Scale Factor')
    map = log10(map*1e3);
elseif strcmp(param_name,'Effective Scatterer Size')
    map = 2*map*1e6;
elseif strcmp(param_name,'HK Scatterer Clustering Param')
    map = log10(map);
end

% Check imaginary
if sum(imag(map(:)))~=0
    fprintf('%s is imaginary. Only retaining real part.\n',param_name);
    map = real(map);
end

end

