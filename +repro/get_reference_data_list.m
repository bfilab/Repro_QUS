function [ref_list,ref_header] = get_reference_data_list(ref_dir)
% Return an array containing details of the reference datasets
% The array will indicate (by column)
%   1) Transducer number (250 or 550)
%   2) Calibration phantom (15 or 18)
%   3) Surface location of phantom (in mm)
%   4) Tx focus location (in mm)
%
%   Inputs
%       ref_dir: folder containing the reference datasets
%  
%   Outputs:
%       ref_list: array (as specified above) with details of the reference
%           datasets
%
%       ref_header: cell array with column names
%

ref_header = {'Transducer','Phantom','Surf. Location','Tx Depth'};

% Get a list of all the reference datasets in the folder
search_str = sprintf('%s\\*RF_Data.mat',ref_dir);
all_files = dir(search_str);

ref_list = zeros(length(all_files),4); % initialize the array

for r_idx=1:length(ref_list)
    this_fid = all_files(r_idx).name;

    % Parse the filename for the acquisition details
    [this_trans, this_phan, this_surfpos, this_tx] = placenta.parse_ref_filename(this_fid);

    ref_list(r_idx,:) = [this_trans, this_phan, this_surfpos, this_tx];

end

% Sort the array in order of: transducer, phantom, surface location, Tx
% depth
ref_list = sortrows(ref_list,[1,2,3,4]);