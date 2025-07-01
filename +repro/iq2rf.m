function [save_fid] = iq2rf(id,num_frames)
% Generates RF data from IQ and XML files. This function does the IQ to RF
% remixing for all .bmode/.xml files in the current working directory

% Get a list of all files in the current folder
all_bmode = dir('*.bmode');

 % all_bmode = dir('*.3d.bmode');
for b_idx=1:length(all_bmode)
    this_file = all_bmode(b_idx).name;
    
    % Need only the file name with no extensions (i.e., remove the
    % .iq.3d.bmode)
    str_tok = strsplit(this_file,'.');
    this_file = str_tok{1};

   
    % Convert IQ to RF
    [rf_data,~,~,fs,params,sysParam] = repro.read_visualsonics(this_file,[],num_frames);
    
    del_mm = params.depth_offset; % store separately, just because 

    % Save the data to a .mat file
    save_fid = sprintf('%s_RF',id);
    save(save_fid,'rf_data','fs','del_mm','params','sysParam','-v7.3','-nocompression');
end
