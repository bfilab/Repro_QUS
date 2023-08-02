function qus_processing_select_multiple_roi(in_fid,qus_frames)
% % QUS processing of RF echo data acquired at Tulane using Vevo 2100
% The basis of this code is Process_SUNY_GE_Data. However, this version
% does not perform batch processing starting at a top-level folder. I can't
% assume the folder structure with data from Tulane will follow a specific
% format. 
%
% Therefore, this function operates on the data specified by in_fid. 
%
%   Inputs
%       in_fid: filename of RF echo data. If not in the working directory,
%           the full path must be specified
%       qus_frames: (optional) specifies with RF echo frames should be
%           processed. If omitted or empty, process all echo frames in the
%           RF dat array

%% QUS Processing
% 2023/1/11 (CH):  Modified from Process_SUNY_GE_Data for use with echo
%   data acquired at Tulane with Vevo 2011

% All GE data are assumed to be a single frame stored in each folder. The
% calibration phantom method is used for BSC estimation and the calibration
% phantom comprises 60um scatterers.

close all; clc;

save_folder = cd;

% in_fid = 'RTG_4_808_2020-05-28-10-49-53_RF_Data.mat';
% in_fid = 'RTG_4_808_2020-05-30-10-19-54_RF_Data.mat';
% in_fid = 'RTG_5_808_2020-06-01-11-15-02_RF_Data.mat';
% qus_frames = 50;

if ~exist('qus_frames','var') || isempty(qus_frames)
    qus_frames = Inf; % default flag to indicate processing all frames later
end

% Set some QUS processing parameters
do_adaptive_bw = false; % should the processing bandwidth be adaptive?
adaptive_bw_thresh = 6; % -dB threshold for selecting adaptive bandwidth

% Specify which calibration phantom to use
calib_phantom = 18; % 15 or 18 (15 has artifacts, best not to use for now)

% Default file name for QUS results contains the processing date
date_str = datestr(datetime,'yyyy-mm-dd');
qus_save_str = sprintf('qus_results_%s',date_str);

% cd('D:\GE_Data_for_Processing');
% cd('E:\Documents\SUNY_LN\GE_Data_for_Processing');

%% 1. Add path to BSC toolbox
% edit this as needed.
% toolbox_path = ['C:\Users\tlye\OneDrive - Riverside Research\',...
%     'Documents\MATLAB\BSC Toolbox'];
toolbox_path = 'C:\Users\acm4005\Box\WCM_Tulane_Shared_Folder\Placenta_QUS\BSC Toolbox';
% toolbox_path = 'C:\Cameron_Matlab_Files\Riverside\BSC Toolbox';

addpath(genpath(toolbox_path));

%% 2. Inputs

% Set the folder containing all of the calibration datasets
%ref_dir = 'F:\OneDrive - med.cornell.edu\Documents\Photoacoustic_Placenta\Calibration_Data';
ref_dir = 'C:\Users\acm4005\Box\WCM_Tulane_Shared_Folder\Phantom Data\Calibration_Data';

% Get a list of the reference datasets
[ref_data_list, ref_data_header] = placenta.get_reference_data_list(ref_dir);

% Load the RF data to be processed
samp_data = load(in_fid);
num_rf_frames = size(samp_data.rf_data,3); % number of RF frames in dataset

% If specific frames were specified for processing, analyze all of them 
if isinf(qus_frames)
    qus_frames = 1:num_rf_frames;
end

if length(qus_frames) > 1
    error('Can only process single frames for now');
end

% Ensure all frames indices for processing fall within the range of
% available frames
qus_frames(qus_frames < 1) = 1; % negative indices not allowed
qus_frames(qus_frames > num_rf_frames) = num_rf_frames; % Can't index beyond available frames
qus_frames = unique(qus_frames); % no redundant QUS processing of frames

% Generate the image axes vectors for later use
[axial_vec,lateral_vec] = placenta.generate_image_axes(samp_data.rf_data,samp_data.sysParam,samp_data.fs*1e6);


%% 3. ROI Info
% set this as desired
% X typically assumed to be the 2nd dimension of the MATLAB matrix
% Y typically assumed to be the 3rd dimension of the MATLAB matrix
% Z typically assumed to be the 1st dimension of the MATLAB matrix

clear roi

% ROI Length (in meters)
% 1mm is approximately 10 wavelengths and 10 A-lines
roi.len_x = 1.75e-3;
roi.len_y = 1;
roi.len_z = 1e-3;

% Note: roi.len y is actually just the interval between planes
%       (since we are processing only one plane for this example)

% ROI Overlap (in percent)
roi.o_x = 75;
roi.o_y = 0;
roi.o_z = 75;

% Portion of the image to process
% All need to be in [m], hence the division by 1000 of the axes vectors
roi.init_z = axial_vec(1)/1000; % Will be adjusted based on the sample surface location
roi.end_z = axial_vec(end)/1000;
roi.init_x = lateral_vec(1)/1000;
roi.end_x = lateral_vec(end)/1000;

% Set these to only process a sub-region of the full image volume.
% If you just want to process the entire image volume, do not set these
% fields (see ROI Functions -> genROIpos.m)
% roi.init_z = 5*1e-3;
% roi.end_z = 8*1e-3;
% roi.init_x = 11*1e-3;
% roi.end_x = 18*1e-3;
% 
% roi.init_z = [];
% roi.end_z = [];
% roi.init_x = [];
% roi.end_z = [];

% roi.end_z = 20e-3; % beyond this point is likely to be low SNR. 
% 
% roi.init_z = 11e-3;
% roi.end_z = 15e-3;
% roi.init_x = 10.5e-3;
% roi.end_x = 12.5e-3;

% % For plotting example BSC curves
% roi.init_x = 12e-3;
% roi.end_x = 15e-3;
% roi.init_z = 10e-3;
% roi.end_z = 15e-3;


%% 4. Spectrum Info
% Padding value for FFT
N_fft = 2048;
% Bandwidth parameters (Hz)
min_bw = 9e6;
max_bw = 19e6;

%% 5. Reference and Compensation Parameters
% This section loads the structs holding the reference acoustical
% properties and compensation parameters. You can go to the directories
% listed here and check the Generate____Obj.m scripts to see how the
% parameters are set.

% Theoretical BSC params for reference
ph_dir = fullfile(toolbox_path,'Scattering Theory Functions');
if calib_phantom == 15
    ph_fname = 'ph_15.mat';
    attn_fname_ref = 'attn_params_15um.mat';
elseif calib_phantom == 18
    ph_fname = 'ph_18.mat';
    attn_fname_ref = 'attn_params_18um.mat';
else
    error('Invalid calibration phantom');
end

% Attenuation params
attn_dir_ref = fullfile(toolbox_path,'Compensation Functions');

attn_dir_samp = fullfile(toolbox_path,'Compensation Functions');
attn_fname_samp = 'attn_params_lung.mat';
% This says lung, but I just used the same parameters for liver and lung.

% Transmission params
transm_dir_ref = fullfile(toolbox_path,'Compensation Functions');
transm_fname_ref = 'transm_attn_params_tpx.mat';

%% 6. Specify desired parameters to save
% These can be anything you want.

params = {'HK Structure Param',...
    'HK Scatterer Clustering Param',...
    'Nak Shape Param',...
    'Nak Scale Factor',...
    'Spectral Slope',...
    'Intercept',...
    'Midband Fit', ...
    'Effective Scatterer Size', ...
    'Acoustic Concentration', ...
    'XCoord', ...
    'YCoord', ...
    'ZCoord', ...
    'BSC_dB', ...
    'env', ...
    'HK_s', ...
    'HK_omega', ...
    'Nak Fit', ...
    'Gauss Fit', ...
    'Freq Band', ...
    'BW Limits', ...
    'SNR'};

% Set which parameters to save in a results table for later recall
results_params = {'HK Structure Param',...
    'HK Scatterer Clustering Param',...
    'Nak Shape Param',...
    'Nak Scale Factor',...
    'Spectral Slope',...
    'Intercept',...
    'Midband Fit', ...
    'Effective Scatterer Size', ...
    'Acoustic Concentration'};

%% 7. Setup reference

% 2023/1/11 (CH) For now, processing only a single frame. Choose the best
% reference dataset based on the single frame. Later, I will need to figure
% out a way to generalize this section
%
% Using the sample data acquisition parameters and the point of the tissue
% surface closest to the transducer, select the appropriate reference
% dataset
samp_surf_seg = samp_data.seg_struct(qus_frames).surf_roi.Position;
samp_surf = min(samp_surf_seg(:,2)); % least distance from transducer to sample surface
%force_samp_surf = 9; % Set the sample surface to 9mm to surface calibration data where surface is at 8mm
ref_fname = placenta.select_ref_data(ref_data_list,samp_data.sysParam,calib_phantom,samp_surf);

% Set the start of the processing region to the sample surface location
% roi.init_z = samp_surf/1000; % needs to be in [m]
% roi.init_z = 13e-3; % start somewhere in the tissue for debugging

% "loadData" formats your data into a struct used for this code
% For your data, you need to add your own case.
% You need to fill out all the fields of the samp struct, as shown for the
% other cases already filled out in the function.
ref         = loadData(ref_dir,ref_fname,'Tulane_Vevo2100');
ref.N_fft   = N_fft;

% Load attenuation data and theoretical reference BSC
load(fullfile(attn_dir_ref,attn_fname_ref));
load(fullfile(ph_dir,ph_fname));

if calib_phantom == 15
    ref.attn_params = attn_params_15;
    ref.ph = ph_15;
elseif calib_phantom == 18
    ref.attn_params = attn_params_18;
    ref.ph = ph_18;
end

% Theoretical reference BSC
ref.f = ref.fs.*(0:(ref.N_fft/2-1))/ref.N_fft;
ref.theo_bsc_db  = 10*log10( computeBSC_mono(ref.f,ref.ph,ref.ph.Nume) );
% 
% Load transmission data
load(fullfile(transm_dir_ref,transm_fname_ref));
ref.transm_params = transm_params_tpx;
% ref.transm_params.d = 25*1e-6;                                              % tpx thickness in meters
ref.transm_params.d = 52e-6; % TPX thickness in [m]

% Compute transmission compensation
if isfield(ref,'transm_params')
    T = compTransm(ref.f, ...
        ref.transm_params, ...
        ref.transm_params.attn,...
        ref.transm_params.d);
    ref.transm_db = 10*log10( ( abs(T{1}) .* abs(T{2}) ).^(-2) );
else
    ref.transm_db = zeros(size(ref.theo_bsc_db));
end

%% 8. Setup sample

this_qus_frame = qus_frames(1);
samp = samp_data;
samp.rf_data = samp.rf_data(:,:,this_qus_frame);
samp.seg_struct = samp.seg_struct(this_qus_frame);
% Convert to struct compatible with this QUS code
samp = placenta.convert_to_qus_struct(samp);

% close all; 

% Plot the B-mode image that will be used for ROI selection
bmode_fig = figure(1); clf;
this_bmode = 20*log10(abs(hilbert(samp.data(:,:,1))));
this_bmode = this_bmode - max(this_bmode(:));
imagesc(lateral_vec,axial_vec,this_bmode,[-50, 0]);
set(gca,'FontSize',16);
xlabel('Lateral (mm)')
ylabel('Axial (mm)')
axis image
colormap gray

while true
    
    samp.adaptive_bw = do_adaptive_bw;
    samp.bw_thresh = adaptive_bw_thresh;

    % Set sample frequency parameters
    samp.N_fft  = N_fft;
    samp.min_bw = min_bw;
    samp.max_bw = max_bw;

    % Attenuation
    % Assume it's 0.5 dB/MHz/cm, n=0, c=0
    load(fullfile(attn_dir_samp,attn_fname_samp));
    samp.attn_params = attn_params_lung;
%     Try using attenuatinon value for the phantom
%     load(fullfile(attn_dir_ref,'attn_params_15um.mat'));
%     samp.attn_params = attn_params_15;

    % Have the user select the ROI location

    figure(bmode_fig);

    roi_pt = drawpoint(gca);
    
    roi_z = roi_pt.Position(2)*1e-3;
    roi_x = roi_pt.Position(1)*1e-3;
    
    delete(roi_pt)
    if exist('roi_rect')
        %delete(roi_rect);
    end
    
    % ROI Length (in meters)
%     roi.len_x = 2e-3;
%     roi.len_y = 1;
%     roi.len_z = 2e-3;
    
    roi.init_x = roi_x - roi.len_x/2;
    roi.end_x = roi_x + roi.len_x/2;
    roi.init_z = roi_z - roi.len_z/2;
    roi.end_z = roi_z + roi.len_z/2; % Increase size slightly to make sure one ROI fits
    
    % Note: roi.len y is actually just the interval between planes
    %       (since we are processing only one plane for this example)
    
    % ROI Overlap (in percent)
    roi.o_x = 75;
    roi.o_y = 0;
    roi.o_z = 75;

    %% 9. Generate ROI
    [all_roi,all_roi_idx] = genROIpos(samp,roi);

    % Plot the ROI on the B-mode image
    figure(bmode_fig);
    hold on;
    roi_rect = rectangle('Position',[roi.init_x,roi.init_z,roi.len_x,roi.len_z]*1000,'EdgeColor','y','LineWidth',2);

    %% 10. Compute QUS parameters for the select ROI

    % All the results will be stored in the cell array, result. Each cell
    % in the cell array will store a matrix/cell array holding the QUS
    % parameters (or other desired parameters) in each ROI, i.e. the
    % parameter maps.
    result = cell(size(params));

    % Process each ROI
    for i = 1:length(all_roi.pos_x)
        for j = 1 % Just do one plane for this example
            for k = 1:length(all_roi.pos_z)


                %% Extract ROI

                % Get ROI coordinates
                roi.pos_x = all_roi.pos_x(i);
                roi.pos_y = all_roi.pos_y(j);
                roi.pos_z = all_roi.pos_z(k);
                
%                     roi_z_mid = range([roi.init_z,roi.end_z])/2 + roi.init_z;
                
                % Sample
                [samp.sub_vol,samp.roi_idx] = getROI(samp,roi);
                result{param2idx(params,'XCoord')}(k,i,j) = mean(samp.sub_vol.x);
                result{param2idx(params,'YCoord')}(k,i,j) = mean(samp.sub_vol.y);
                result{param2idx(params,'ZCoord')}(k,i,j) = mean(samp.sub_vol.z);

                % Get Sample Surface
    %             surf_vec_idx = samp.surf_mat(samp.roi_idx.x_vec,j);
    %             surf_pos_idx = round(mean(surf_vec_idx));
    %             samp.surf_pos = samp.z(surf_pos_idx);
%                 samp.surf_pos = samp.z(1);
                
                % Reference
                ref_roi = roi;
                % These are set to empty vectors so that the full lateral
                % extent of the phantom volume is used, see getROI.m.
%                 ref_roi.pos_x = [];
                ref_roi.pos_y = [];
%                 ref_roi.len_x = [];
                ref_roi.len_y = [];
                [ref.sub_vol,ref.roi_idx] = getROI(ref,ref_roi);

                %% BSC Estimation
                [samp,ref] = computeBSC_RefPhantom(samp,ref,roi);
                result{param2idx(params,'BSC_dB')}{k,i,j} = samp.sub_vol.est_bsc_db;

                % Compute the SNR
%                 result{param2idx(params,'SNR')}{k,i,j} = compute_ge_snr(samp.f,samp.sub_vol.avg_PS_db);
                result{param2idx(params,'SNR')}{k,i,j} = 0; % Don't worry about computing SNR for now
                
                  %% Power Spectra
%                     s_max = max(samp.sub_vol.avg_PS_db,[],'all');
%                     figure;
%                     plot(samp.f/1e6,samp.sub_vol.avg_PS_db)
%                     axis([5 30 s_max-20 s_max])
%                     ylabel('Power (dB)')
%                     xlabel('Frequency (MHz)')
%                     title('Sample Power Spectrum')
    
                %% Linear Model
%                 [SS, I0, MF] = compLinFit(samp.f,...
%                     samp.sub_vol.est_bsc_db,...
%                     [samp.min_bw, samp.max_bw],...
%                     (samp.min_bw+samp.max_bw)/2 );
%                 result{param2idx(params,'Spectral Slope')}(k,i,j) = SS;
%                 result{param2idx(params,'Intercept')}(k,i,j) = I0;
%                 result{param2idx(params,'Midband Fit')}(k,i,j) = MF;
% 
% %                     %% Gaussian Form Factor
% %                     [ac, aeff,gauss_fit,fit_freq] = est_GaussFF_GlassPlate(samp.f, ...
% %                         samp.sub_vol.est_bsc_db, ...
% %                         [samp.min_bw, samp.max_bw], ...
% %                         samp.c, roi_z_mid, roi.len_z, 'GE');
% %                     result{param2idx(params,'Effective Scatterer Size')}(k,i,j) = aeff;
% %                     result{param2idx(params,'Acoustic Concentration')}(k,i,j) = ac;
% %                     result{param2idx(params,'Gauss Fit')}{k,i,j} = gauss_fit;
% 

%                     %% Gaussian Form Factor (ref phantom method)
                [ac, aeff,gauss_fit,fit_freq] = est_GaussFF(samp.f, ...
                    samp.sub_vol.est_bsc_db, ...
                    [samp.min_bw, samp.max_bw], ...
                    samp.c);
                result{param2idx(params,'Effective Scatterer Size')}(k,i,j) = aeff;
                result{param2idx(params,'Acoustic Concentration')}(k,i,j) = ac;
                result{param2idx(params,'Gauss Fit')}{k,i,j} = gauss_fit;
                result{param2idx(params,'Freq Band')}{k,i,j} = fit_freq;
                
                if isfield(samp,'adaptive_bw') && samp.adaptive_bw == true
                    result{param2idx(params,'BW Limits')}{k,i,j} = [samp.min_bw, samp.max_bw];
                end

                % Get bandwidth indices
                [~,tf1] = min(abs(samp.min_bw - samp.f));
                [~,tf2] = min(abs(samp.max_bw - samp.f));
                sub_f = samp.f(tf1:tf2);
                sub_est_bsc = samp.sub_vol.est_bsc_db(tf1:tf2);

                gs_fig = figure(860+01); %clf;
                %plot(sub_f./1e6,sub_est_bsc,'k','LineWidth',2);
                plot(sub_f./1e6,sub_est_bsc,'LineWidth',2);
                hold on;
                %plot(sub_f./1e6,gauss_fit,':c','LineWidth',2);
                plot(sub_f./1e6,gauss_fit,':','LineWidth',2);
                xlabel('Frequency (MHz)');
                ylabel('BSC (dB/(m\bulletsr)');
                %ylabel('W_{ROI}(f) (dB)');
                legend('Measured','Gaussian FF','Location','best');
                set(gca,'FontSize',14)
                xlim([sub_f(1)-0.2e6, sub_f(end)+0.2e6]./1e6)

%                 lin_fit = sub_f/1e6.*SS + I0;
%                 lin_fig = figure(861+10); clf;
%                 plot(sub_f./1e6,sub_est_bsc,'k','LineWidth',2);
%                 hold on;
%                 plot(sub_f./1e6,lin_fit,':r','LineWidth',2);
%                 xlabel('Frequency (MHz)');
%                 ylabel('W_{ROI}(f) (dB)');
%                 legend('Measured','Linear Fit','Location','best');
%                 set(gca,'FontSize',14)
%                 xlim([sub_f(1)-0.2e6, sub_f(end)+0.2e6]./1e6)
%                 
% 
%                 
%                 %% Envelope Statistics
% 
%                 % Get envelope ROI
%                 samp.env_sub_vol = samp.env( samp.roi_idx.z_vec, ...
%                     samp.roi_idx.x_vec, ...
%                     samp.roi_idx.y_vec);
%                 result{param2idx(params,'env')}{k,i,j} = samp.env_sub_vol;
% 
%                 % Compute Nakagami parameters
%                 [nak_param,scale_factor,nak_fit] = computeNakagami(samp.env_sub_vol);
%                 result{param2idx(params,'Nak Shape Param')}(k,i,j) = nak_param;
%                 result{param2idx(params,'Nak Scale Factor')}(k,i,j) = scale_factor;
%                 result{param2idx(params,'Nak Fit')}(k,i,j) = nak_fit;
% 
%                 % Compute Homodyned-K parameters
%                 [hk_scat_clust_param, hk_struc_param, s, omega] = computeHomodynedK(samp.env_sub_vol);
%                 result{param2idx(params,'HK Scatterer Clustering Param')}(k,i,j) = hk_scat_clust_param;
%                 result{param2idx(params,'HK Structure Param')}(k,i,j) = hk_struc_param;
%                 result{param2idx(params,'HK_s')}(k,i,j) = s;
%                 result{param2idx(params,'HK_omega')}(k,i,j) = omega;
%                 
%                 a = 1;
%                 
% %                 % Plot the histogram with the pdf fits
%                     hist_fig = figure(900+10); clf;
%                     env_hist = histogram(samp.env_sub_vol(:),'NumBins',70,'Normalization','probability');
%                     [hk_x,hk_y] = plotHomodynedK(hk_scat_clust_param,s,omega,env_hist);
%                     [nak_x,nak_y] = plotNakagami(nak_fit,env_hist);
%                     hold on
%                     plot(hk_x,hk_y,'--r','LineWidth',3);
%                     hold on
%                     plot(nak_x,nak_y,'-.b','LineWidth',3);
%                     %legend('Envelope Histogram','H-k Fit');
%                     legend('Envelope Histogram','H-k Fit','Nakagami Fit');
%                     xlabel('Envelope Value (mV)')
%                     ylabel('Probability');
%                     set(gca,'FontSize',14);
%                     
%                     a = 1;
%                     
%                     {'ESD', aeff; 'EAC', ac; 'H-k alpha', hk_scat_clust_param; 
%                         'H-k k', hk_struc_param; 'Nak alpha', nak_param; 
%                         'Nak Omega', scale_factor; 'ROI X',roi.pos_x*1000; 
%                         'ROI Z',roi.pos_z*1000}
                    
            end
        end
    end

%     % Create the results vector and save the QUS parameter
%     % estimates for this ROI to a .mat file
%     results_vec = zeros(length(results_params),1);
%     for p_idx=1:length(results_params)
%         results_vec(p_idx) = result{param2idx(results_params,results_params{p_idx})}(k,i,j);
%     end
%     results_table = table(results_vec,'RowNames',results_params);
%     roi_results = result(k,i,j);
% 
%     %save('_roi_results','results_table','roi_results');
%     save('lz_roi_results','samp');

    
end

%%
% saveas(bmode_fig,'_bmode_with_roi','fig');
% saveas(bmode_fig,'_bmode_with_roi','png');
% saveas(gs_fig,'_gs_fit_fig','fig')
% saveas(gs_fig,'_gs_fit_fig','png');
% saveas(lin_fig,'_lin_fit_fig','fig');
% saveas(lin_fig,'_lin_fit_fig','png');
% saveas(hist_fig,'_hist_fig_fig','fig');
% saveas(hist_fig,'_hist_fig_fig','png');
