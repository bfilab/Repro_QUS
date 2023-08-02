function [data_cleaned] = cleanValues(data,varargin)
%CLEANVALUES Remove NaNs, Infs, and imaginary numbers from the data.
% INPUTS:
%   data         = 1D or 2D data
%   varargin{1}  = feature_name, for printing warnings
%   varargin{2}  = dim, dimension of features
%                  Each feature has a 1D vector of data
%                  Clean the data so that each feature still has the same
%                  number of data points, without removing features.
%                  e.g. if we have 4 features, each with 100 samples
%                       and a 4x100 matrix, then dim = 1.
% OUTPUTS:
%   data_cleaned = data with NaNs and Infs removed

% 10/15/2020 (THL): Created
% 10/21/2020 (THL): Corrected a bug where imaginary numbers were not be
%                   detected properly (because isreal and isimag don't work
%                   per element).

%% Print warnings
if ~isempty(varargin)
    
    % Get feature name
    feature_name = varargin{1};
    
    if sum(isinf(data(:)))
        warning('%s had Inf values',feature_name);
    end
    if sum(isnan(data(:)))
        warning('%s had NaN values',feature_name);
    end
    if sum(double(imag(data(:))~=0))
        warning('%s had imaginary numbers',feature_name);
    end
    
end

%% Clean data
if isvector(data)
    % vector data
    data(isinf(data)) = [];
    data(isnan(data)) = [];
    data(imag(data)~=0) = [];
    data_cleaned = real(data);
else
    
    % check dimensions of matrix data
    if length(varargin)<2||isempty(varargin{2})
        error('Must input dim to clean matrix data.');
    end
    feat_dim = varargin{2};
    if feat_dim~=1&&feat_dim~=2
        error('dim must be 1 or 2');
    end

    % get indices of bad data points
    all_bad_idx = [];
    for i = 1:size(data,feat_dim)
        if feat_dim == 1, data_line = data(i,:);
        else,             data_line = data(:,i);
        end
        inf_idx = find(isinf(data_line));
        nan_idx = find(isnan(data_line));
        imag_idx = find(imag(data_line)~=0);
        all_bad_idx = unique([all_bad_idx;inf_idx;nan_idx;imag_idx]);
    end
    
    % clear bad data points for all features
    if feat_dim==1, data(:,all_bad_idx) = [];
    else          , data(all_bad_idx,:) = [];
    end
    data_cleaned = real(data);
    
end

end

