%% Rename the .bmode or .xml file

% f_suffix = '.iq.3d.bmode';
% search_str = '*.bmode';

f_suffix = '.iq.xml';
search_str = '*.xml';

% f_prefix = 'LZ550';

f_prefix = 'LZ250';

bmode_list = dir(search_str);

for b_idx=1:length(bmode_list)
    this_fid = bmode_list(b_idx).name;
    % Get the scatterer size, surface position, and Tx focus position
    scat_size = regexpi(this_fid,'\d{1,2}um','match','once');
    mm_strs = regexpi(this_fid,'\d{1,2}mm','match');

    % For LZ-550 data
%     surf_pos = mm_strs{1}; % surface location is first token
%     tx_pos = mm_strs{2}; % Tx focus is second token

    % For LZ-250 data
    tx_pos = mm_strs{1}; % Tx focus is first token
    surf_pos = mm_strs{2}; % surface location is second token

    % Create the new file name
    save_str = sprintf('%s_%sPhantom_%sSurf_%sTx%s', ...
        f_prefix,scat_size,surf_pos,tx_pos,f_suffix);

    movefile(this_fid,save_str);

end