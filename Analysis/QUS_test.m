%For processing QUS data from 3D bmode scans of reproductive tract tissues
%6/4/2024
%Updated 7/3/2025
ProcessIQ = false;
Segment = false;
ProcessQUS = false; 
OverlayQUS = true;
Compile = false;

addpath(genpath('D:\Andrew\Repro_QUS\'));
data_dir = 'D:\Andrew\Prolapse Model';
cd(data_dir);
directories = {dir('M955*'),dir('M949*')};
%date_str = char(datetime(datetime,"Format","uuuu-MM-dd"));
date_str = '2025-07-07';
anatomy = readtable('anatomy_qus.csv');

for d = 1:length(directories)
    samples = directories{d};
    for i = 1:length(samples)
        %Go to sample folder
        fid = samples(i).name;
        cd(fid);
        anatomy_table = anatomy(anatomy.ID == str2num(fid(2:end)),:);
        anatomy_array = table2array(anatomy_table(1,3:6));
    
        %Convert IQ data file into .mat file with RF data
        if ProcessIQ
            mat_file = repro.iq2rf(fid);
        else
            mat_file = strcat(fid,'_RF');
        end
          
        %Draw an ROI around the tissue, then a line on the surface
        if Segment
            for a = anatomy_array
                repro.segment_repro_frames(mat_file,a);
            end
        end
        
        %Compute the QUS parameters in the entire image
        
        for j = 1:length(anatomy_array)
            organ = "";
            frame = anatomy_array(j);
            switch j
                case 1
                    organ = 'vagina';
                case 2
                    organ = 'external';
                case 3
                    organ = 'internal';
                case 4
                    organ = 'uterus';
            end
            qus_file = strcat(date_str,'_',organ,'_',num2str(frame));

            if ProcessQUS
                repro.qus_processing(mat_file,frame,qus_file);
            end

            if OverlayQUS
                repro.plot_bmode_qus_overlay(qus_file,frame);
            end

            if Compile
                fid_cell = {strcat(data_dir,'\',fid,'\',mat_file,'.mat')};
                id_cell = {fid};
                repro.compile_qus_results(qus_file,fid_cell,id_cell,frame,organ)
            end
        end
    
        cd(data_dir);
    end
end

%For all of the samples
%Compile qus results into a table for data analysis
%List of paths to each RF data file
% fid_cell_array = {'E:\QPA\data2\G3_RF.mat';'E:\QPA\data\G0_RF.mat'};
%List of data IDs to lab
% id_cell_array = {'M955';'M830'};
% repro.compile_qus_results(qus_file,fid_cell_array,id_cell_array,frame_num);