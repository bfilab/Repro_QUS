function compile_qus_results(qus_fid,fid_cell_array,id_cell_array,qus_frame)
% Compile all the QUS results into a single table and write to a .csv file

%%
% Create a cell array containing the file names (with full paths) to the
% .mat files containing the RF datasets. Ideally, this will be
% generated programmatically
%{
fid_cell_array = {'F:\OneDrive - med.cornell.edu\Documents\Photoacoustic_repro\Preeclampsia Model\RTG11\Day 14\RTG_4_808_2020-05-28-10-49-53_RF_Data.mat'; 
    'F:\OneDrive - med.cornell.edu\Documents\Photoacoustic_repro\Preeclampsia Model\RTG11\Day 16\RTG_4_808_2020-05-30-10-19-54_RF_Data.mat'; 
    'F:\OneDrive - med.cornell.edu\Documents\Photoacoustic_repro\Preeclampsia Model\RTG11\Day 18\RTG_5_808_2020-06-01-11-15-02_RF_Data.mat'};
%}

% Name of the .mat file containing the QUS results
%qus_fid = 'qus_results_2023-01-12.mat';

% Frame number in the dataset that was processed
%qus_frame = 1; % when processing the data for the R21 proposal, I had chosen frame 50
                % I suggest sticking with frame 1 when the data are 2D

% Cell array containing the ID for each dataset
% I suggest following the syntax {Animal}_{GestationDay}_{Frame}
% Whatever format is chosen, be consistent! Otherwise, post-processing is
% going to become a headache
%id_cell_array = {'RTG11_14_1'; 'RTG11_16_1';'RTG_18_1'};

% Name of the .csv the results will be written to
% Note: this will overwrite any .csv with the same name!
% Include the full folder path if it should not be written to the current
% working directory
%date_str = datestr(datetime,'yyyy-mm-dd');
date_str = char(datetime(datetime,'Format','yyyy-MM-dd'));
csv_save_fid = [date_str,'.csv'];
% csv_save_fid = '2023-07-25.csv';

% Loop over all the datasets and combine the QUS results into a single
% table
full_res_table = table;

for f_idx=1:length(fid_cell_array)
    % for q = 1:length(qus_fid)      
    %     try 
    %         this_res_table = repro.get_qus_results(fid_cell_array{f_idx},qus_fid(q),qus_frame,id_cell_array{f_idx});
    %     catch
    %         continue;
    %     end
    % 
    %     disp(fid_cell_array{f_idx});
    %     break;
    % end
    this_res_table = repro.get_qus_results(fid_cell_array{f_idx},qus_fid,qus_frame,id_cell_array{f_idx});
    full_res_table = [full_res_table; this_res_table];
end

writetable(full_res_table,csv_save_fid);