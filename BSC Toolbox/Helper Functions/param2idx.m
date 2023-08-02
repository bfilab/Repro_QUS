function [idx] = param2idx(params,param_name)
%PARAM2IDX Given a parameter name, give the index of the cell in "params"
%that has a matching string.
% INPUTS:
% params = cell array of parameter names
% param_name = string, parameter name to search for
% OUTPUTS:
% idx = index such that params{idx} = param_name

% 6/22/2020 (THL): Created

out = cellfun(@(x) strcmp(x, param_name), params );
idx = find(out==1);

end

