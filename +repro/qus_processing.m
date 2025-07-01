function qus_processing(in_fid,qus_frames,qus_save_str)
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

[save_folder,~,~] = fileparts(in_fid);  % the results in the same folder as the RF data
if isempty(save_folder)
    % Assume data should be saved in current working directory
    save_folder = cd;
end


% in_fid = 'RTG_4_808_2020-05-28-10-49-53_RF_Data.mat';
% in_fid = 'RTG_4_808_2020-05-30-10-19-54_RF_Data.mat';
% in_fid = 'RTG_5_808_2020-06-01-11-15-02_RF_Data.mat';
% in_fid = 'RTG_4_808_2020-05-30-10-19-54_RF_Data.mat';

% 
if ~exist('qus_frames','var') || isempty(qus_frames)
    qus_frames = Inf; % default flag to indicate processing all frames later
end

% Set some QUS processing parameters
do_adaptive_bw = false; % should the processing bandwidth be adaptive?
adaptive_bw_thresh = 6; % -dB threshold for selecting adaptive bandwidth

% Specify which calibration phantom to use
calib_phantom = 15; % 15 or 18 (15 has artifacts, best not to use for now)

% Default file name for QUS results contains the processing date
%date_str = datestr(datetime,'yyyy-mm-dd');
%qus_save_str = sprintf('qus_results_%s',date_str);

% cd('D:\GE_Data_for_Processing');
% cd('E:\Documents\SUNY_LN\GE_Data_for_Processing');

%% 1. Add path to BSC toolbox
% edit this as needed.
% toolbox_path = ['C:\Users\tlye\OneDrive - Riverside Research\',...
%     'Documents\MATLAB\BSC Toolbox'];
toolbox_path = 'C:\Users\lshi6\Box\QUS\QUS_prolapse\Repro_QUS\BSC Toolbox';
%toolbox_path = 'E:\Cameron_Matlab_Files\Riverside\BSC Toolbox';
% toolbox_path = 'C:\Cameron_Matlab_Files\Riverside\BSC Toolbox';

addpath(genpath(toolbox_path));

%% 2. Inputs

% Set the folder containing all of the calibration datasets
%ref_dir = 'F:\OneDrive - med.cornell.edu\Documents\Photoacoustic_repro\Calibration_Data';
ref_dir = 'C:\Users\lshi6\Box\QUS\QUS_prolapse\Calibration_Data';

% Get a list of the reference datasets
[ref_data_list, ref_data_header] = repro.get_reference_data_list(ref_dir);

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
[axial_vec,lateral_vec] = repro.generate_image_axes(samp_data.rf_data,samp_data.sysParam,samp_data.fs*1e6);


%% 3. ROI Info
% set this as desired
% X typically assumed to be the 2nd dimension of the MATLAB matrix
% Y typically assumed to be the 3rd dimension of the MATLAB matrix
% Z typically assumed to be the 1st dimension of the MATLAB matrix

clear roi

% ROI Length (in meters)
% 1mm is approximately 10 wavelengths and 10 A-lines
%roi.len_x = 1.75e-3;
roi.len_x = 1.0e-3;
%roi.len_x = 5.25e-3;
roi.len_y = 1;
roi.len_z = 1.0e-3;
%roi.len_z = 1.25e-3;

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
min_bw = 10e6;
max_bw = 35e6;
%max_bw = 29e6;

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
    'SNR', ...
    'ROI Z Dist'};

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
force_samp_surf = 4; % Set the sample surface to 4mm to surface calibration data where surface is at 4mm
samp_surf = force_samp_surf;
ref_fname = repro.select_ref_data(ref_data_list,samp_data.sysParam,calib_phantom,samp_surf);


% Set the start of the processing region to the sample surface location
roi.init_z = samp_surf/1000; % needs to be in [m]
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
% This is where the loop starts for batch processing

top_dir = cd;
for f_idx=1:length(qus_frames)
    
    close all; 
    
    this_qus_frame = qus_frames(f_idx);
    
    samp = samp_data;
    samp.rf_data = samp.rf_data(:,:,this_qus_frame);
    samp.seg_struct = samp.seg_struct(this_qus_frame);

    % Convert to struct compatible with this QUS code
    samp = repro.convert_to_qus_struct(samp);
    
    samp.adaptive_bw = do_adaptive_bw;
    samp.bw_thresh = adaptive_bw_thresh;

    % Set sample frequency parameters
    samp.N_fft  = N_fft;
    samp.min_bw = min_bw;
    samp.max_bw = max_bw;

    % Attenuation
    % Assume it's 0.5 dB/MHz/cm, n=1, c=0
    load(fullfile(attn_dir_samp,attn_fname_samp));
    samp.attn_params = attn_params_lung;
    %Try using attenuatinon value for the phantom
%     load(fullfile(attn_dir_ref,'attn_params_15um.mat'));
%     samp.attn_params = attn_params_15;

    % Load detected surface of volume
    % load('samp_surf');
    % samp.surf_mat = surf_mat;

    %% 9. Generate ROI
    [all_roi,all_roi_idx] = genROIpos(samp,roi);

    plotROIonBmode(samp, all_roi);
    saveas(gcf,'bmode_roi','fig');

    %% 10. Compute QUS parameters for each ROI

    tic

    f = waitbar(0,'Processing ROI');
    iter_num = 1;
    total_iter = length(all_roi.pos_x)*...
        1*...
        length(all_roi.pos_z);


    % All the results will be stored in the cell array, result. Each cell
    % in the cell array will store a matrix/cell array holding the QUS
    % parameters (or other desired parameters) in each ROI, i.e. the
    % parameter maps.
    result = cell(size(params));

    % Process each ROI
    for i = 1:length(all_roi.pos_x)
        for j = 1 % Just do one plane for this example
            for k = 1:length(all_roi.pos_z)

                % Update waitbar
                waitbar(iter_num/total_iter,f,'Processing ROI');
                iter_num = iter_num+1;

                %% Extract ROI

                % Get ROI coordinates
                roi.pos_x = all_roi.pos_x(i);
                roi.pos_y = all_roi.pos_y(j);
                roi.pos_z = all_roi.pos_z(k);
                
%                     roi_z_mid = range([roi.init_z,roi.end_z])/2 + roi.init_z;

                % Based on the lateral position of the ROI, adjust the
                % defined z-position of the sample using the skin
                % segmentation curve
                roi_mid_x = roi.pos_x + roi.len_x/2; % mid-lateral point of ROI
                roi_mid_x = roi_mid_x*1000; % convert to mm for compatibility with segmentation data
                surf_seg = samp_data.seg_struct(this_qus_frame).surf_roi.Position;

                this_samp_surf = interp1(surf_seg(:,1),surf_seg(:,2),roi_mid_x);
                this_samp_surf = this_samp_surf/1000; % back to [m] 
                samp.surf_pos = this_samp_surf;

                % Also compute the distance from the sample surface to
                % z-axis midpoint of the ROI. Will save this as a variable.
                roi_mid_z = roi.pos_z + roi.len_z/2;
                roi_z_dist = roi_mid_z - samp.surf_pos;

                result{param2idx(params,'ROI Z Dist')}(k,i,j) = roi_z_dist;
                
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
                ref_roi.pos_x = [];
                ref_roi.pos_y = [];
                ref_roi.len_x = [];
                ref_roi.len_y = [];
                [ref.sub_vol,ref.roi_idx] = getROI(ref,ref_roi);


                %% BSC Estimation
                [samp,ref] = computeBSC_RefPhantom(samp,ref,roi);
                result{param2idx(params,'BSC_dB')}{k,i,j} = samp.sub_vol.est_bsc_db;

                % Compute the SNR
%                 result{param2idx(params,'SNR')}{k,i,j} = compute_ge_snr(samp.f,samp.sub_vol.avg_PS_db);
                result{param2idx(params,'SNR')}{k,i,j} = 0; % Don't worry about computing SNR for now
                
                %% Linear Model
                [SS, I0, MF] = compLinFit(samp.f,...
                    samp.sub_vol.est_bsc_db,...
                    [samp.min_bw, samp.max_bw],...
                    (samp.min_bw+samp.max_bw)/2 );
                result{param2idx(params,'Spectral Slope')}(k,i,j) = SS;
                result{param2idx(params,'Intercept')}(k,i,j) = I0;
                result{param2idx(params,'Midband Fit')}(k,i,j) = MF;

%                     %% Gaussian Form Factor
%                     [ac, aeff,gauss_fit,fit_freq] = est_GaussFF_GlassPlate(samp.f, ...
%                         samp.sub_vol.est_bsc_db, ...
%                         [samp.min_bw, samp.max_bw], ...
%                         samp.c, roi_z_mid, roi.len_z, 'GE');
%                     result{param2idx(params,'Effective Scatterer Size')}(k,i,j) = aeff;
%                     result{param2idx(params,'Acoustic Concentration')}(k,i,j) = ac;
%                     result{param2idx(params,'Gauss Fit')}{k,i,j} = gauss_fit;


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
%                     [~,tf1] = min(abs(samp.min_bw - samp.f));
%                     [~,tf2] = min(abs(samp.max_bw - samp.f));
%                     sub_f = samp.f(tf1:tf2);
%                     sub_est_bsc = samp.sub_vol.est_bsc_db(tf1:tf2);
% 
%                     figure(860); clf;
%                     plot(sub_f./1e6,sub_est_bsc,'k','LineWidth',2);
%                     hold on;
%                     plot(sub_f./1e6,gauss_fit,':r','LineWidth',2);
%                     xlabel('Frequency (MHz)');
%                     ylabel('W_{ROI}(f) (dB)');
%                     legend('Measured','Gaussian FF','Location','best');
%                     set(gca,'FontSize',14)
%                     xlim([sub_f(1)-0.2e6, sub_f(end)+0.2e6]./1e6)
                

                
                %% Envelope Statistics

                % Get envelope ROI
                samp.env_sub_vol = samp.env( samp.roi_idx.z_vec, ...
                    samp.roi_idx.x_vec, ...
                    samp.roi_idx.y_vec);
                result{param2idx(params,'env')}{k,i,j} = samp.env_sub_vol;

                % Compute Nakagami parameters
                [nak_param,scale_factor,nak_fit] = computeNakagami(samp.env_sub_vol);
                result{param2idx(params,'Nak Shape Param')}(k,i,j) = nak_param;
                result{param2idx(params,'Nak Scale Factor')}(k,i,j) = scale_factor;
                result{param2idx(params,'Nak Fit')}(k,i,j) = nak_fit;

                % Compute Homodyned-K parameters
                [hk_scat_clust_param, hk_struc_param, s, omega] = computeHomodynedK(samp.env_sub_vol);
                result{param2idx(params,'HK Scatterer Clustering Param')}(k,i,j) = hk_scat_clust_param;
                result{param2idx(params,'HK Structure Param')}(k,i,j) = hk_struc_param;
                result{param2idx(params,'HK_s')}(k,i,j) = s;
                result{param2idx(params,'HK_omega')}(k,i,j) = omega;
                
                a = 1;
                
%                 % Plot the histogram with the pdf fits
%                     figure(900); clf;
%                     env_hist = histogram(samp.env_sub_vol(:),'NumBins',70,'Normalization','probability');
%                     [hk_x,hk_y] = plotHomodynedK(hk_scat_clust_param,s,omega,env_hist);
%                     [nak_x,nak_y] = plotNakagami(nak_fit,env_hist);
%                     hold on
%                     plot(hk_x,hk_y,'--r','LineWidth',3);
%                     hold on
%                     plot(nak_x,nak_y,'-.b','LineWidth',3);
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
%                     
                a = 1;
            end
        end
    end

    close(f);

%         save('.\qus_results.mat','result','params')

    %% Plot the parameter maps

    params2plot = {'HK Structure Param',...
        'HK Scatterer Clustering Param',...
        'Nak Shape Param',...
        'Nak Scale Factor',...
        'Spectral Slope',...
        'Intercept',...
        'Midband Fit', ...
        'Effective Scatterer Size', ...
        'Acoustic Concentration'};

    title_labels = {['HK Structure Parameter (',char(954),')'],...
        ['HK Scatterer Clustering Parameter (',char(945),' (log))'],...
        'Nak Shape Param (m)',... % note: m is not meters; just the symbol m
        ['Nak Scale Factor (',char(937),' (log(mV^2))'],...
        'Spectral Slope (dB/sr/m/MHz)',...
        'Intercept (dB/sr/m)',...
        'Midband Fit (dB/sr/m)',...
        'Effective Scatterer Diameter (um)',...
        'Acoustic Concentration (dB/m^3)',...
        };

    figure(2);
    for i = 1:length(params2plot)

        % Get parameter name
        param_fname = params2plot{i};
        param_title = title_labels{i};

        % Get coordinates
        x = result{param2idx(params,'XCoord')};
        x = x(1,:)*1e3;
        z = result{param2idx(params,'ZCoord')};
        z = z(:,1)*1e3;

        subplot(3,3,i)
        param_map = result{param2idx(params,params2plot{i})};
        param_map = convertUnits(param_map,param_fname);
        imagesc(x,z,param_map);
        title(param_title);
        colorbar;
        axis image;

    end
    
    %save_results_str = sprintf('%s\\%s',save_folder,qus_save_str);
    save_results_str = qus_save_str;
    save(save_results_str,'result','params','fit_freq', ...
    'min_bw','max_bw','all_roi','all_roi_idx','do_adaptive_bw', ...
    'adaptive_bw_thresh','-v7.3');
    
        %save_fig_str = sprintf('%s\\all_param_maps_corrected',save_folder);
        save_fig_str = 'all_param_maps_corrected'
        saveas(gcf,save_fig_str,'fig');

    
        cd(top_dir);
end
