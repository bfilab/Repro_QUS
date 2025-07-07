function plot_bmode_qus_overlay(qus_results_fid,frame_num)
% Plot the B-mode image with QUS parameter color overlay. The parameter
% map covers only the repro
%
%   Inputs
%       in_fid: filename (with full path) of data to display
%
%       frame_num: frame in the RF data array to display
%
%       qus_params_to_plot: cell array specifying with QUS parameters to
%           overlay. The total number of figures generated is equal to the
%           length of this cell array
%% 
%%

clc; close all;

%%ADDED BY ANDREW
%Boolean to indicate whether to add lines to images to indicate the bounds
%of the included data around the transducer focus
AddLines = false; 

top_dir = cd;
date_str = char(datetime(datetime,'Format','yyyy-MM-dd'));
%date_str = datestr(datetime,'yyyy-mm-dd');
%date_str = '2023-07-21'

copy_dir = [cd,'\',qus_results_fid,'_bmode_overlay'];
mkdir(copy_dir);




% in_fid = 'RTG_5_808_2020-06-01-11-15-02_RF_Data.mat';
% in_fid = 'RTG_4_808_2020-05-28-10-49-53_RF_Data.mat';
%day = 18

%in_fid = dir('*RF.mat');
in_fid = dir('LZ550_18umPhantom_6mmSurf_7mmTx_RF_Data*');
in_fid = in_fid(1).name;

%qus_results_fid = '.\qus_results_2023-03-22.mat';
%frame_num = 1;

cd(top_dir); % move to the main folder

% ge_folders = get_ge_folders_to_analyze;

qus_params = {'HK Structure Param',...
    'HK Scatterer Clustering Param',...
    'Nak Shape Param',...
    'Nak Scale Factor',...
    'Spectral Slope',...
    'Intercept',...
    'Midband Fit', ...
    'Effective Scatterer Size', ...
    'Acoustic Concentration'};

% qus_params_to_plot = {'Effective Scatterer Size','Intercept','HK Scatterer Clustering Param'};
qus_params_to_plot = qus_params;
% qus_clims = {[40,60],[10,40],[0 3]};
% qus_clims = {[45,70],[90, 107],[10,40],[0 3]};
% qus_str = {'ESD','Int','HK-alpha'};
% qus_str = {'ESD','EAC','Int','HK-alpha'};
qus_str = {'HK-alpha','log10(HK-k)','Naka-m','log10(Naka-Omega)','SS','I0','MBF','ESD','EAC'};
%qus_clims = {[0, 1],[0, 5],[0, 1.2],[5, 8],[0, 1.5],[0, 20],[12, 30],[12, 30],[115, 140]};
qus_clims = {[0, 0.1],[0, 2],[0.6, 1.1],[5, 7],[0, 1],[0, 25],[8, 18],[15, 20],[100, 140]};

% Red/Green colormap
qus_cmap = {[1 0 0; 0 1 0], [1 0 0; 0 1 0]};

% Pink/Yellow colormap for colorblindness
qus_cmap = {[254, 254, 98]./255, [254, 254, 98]./255};

qus_mean = zeros(length(qus_params),1);
qus_std = zeros(length(qus_params),1);
qus_roi_N = zeros(length(qus_params),1);
qus_roi_zero = zeros(length(qus_params),1); % keep track of the number of omitted QUS ROIs

qus_row_count = 0;

% keyboard; % codes needs tested and debugged

% Load the RF echo data and generate the image axes
load(in_fid);
[axial_vec,lateral_vec,c,~] = repro.generate_image_axes(rf_data,sysParam,fs*1e6);
% Keep only the specified frame
rf_data = rf_data(:,:,frame_num);

this_bmode = 20*log10(abs(hilbert(rf_data)));
this_bmode = this_bmode - max(this_bmode(:));

% Load repro segmentation data
roi_poly = seg_struct(frame_num).p_roi.Position; 

% For compatiblity with the rest of the code, convert the units of the
% segmentation points from [mm] to [m]
roi_poly = roi_poly./1000;


%%ADDED BY ANDREW
%Iterate through different QUS fids until the right one is found.
% for q = 1:length(qus_results_fid)      
%     try 
%         load(qus_results_fid(q));    
%     catch
%         continue;
%     end
%     disp(qus_results_fid(q));
%     break;
% end

load(qus_results_fid); 

% Load the QUS results
%load(qus_results_fid);

%%
% Find all the QUS processing ROIs whose (z,x) midpoint falls within
% the node ROI
num_z = length(all_roi.pos_z);
num_x = length(all_roi.pos_x);

% List of the indices of the ROIs whose midpoints are within the node
% ROI
keep_z_inds = [];
keep_x_inds = [];
keep_count = 0;




%[xm,zm] = meshgrid(all_roi_x,all_roi_z);

%[in_roi_idx] = inpolygon(xm(:),zm(:),roi_poly(:,1),roi_poly(:,2));


% all_roi_z1 = all_roi.pos_z(:);
% all_roi_z2 = all_roi.pos_z(:) + all_roi.len_z;
% all_roi_x1 = all_roi.pos_x(:);
% all_roi_x2 = all_roi.pos_x(:) + all_roi.len_x;

all_roi_z = all_roi.pos_z(:)+ all_roi.len_z/2;
all_roi_x = all_roi.pos_x(:)+ all_roi.len_x/2;
[xm,zm] = meshgrid(all_roi_x,all_roi_z);

%[xm,zm] = meshgrid(all_roi_x,all_roi_z);

[in_roi_idx] = inpolygon(xm(:),zm(:),roi_poly(:,1),roi_poly(:,2));

%ADDED BY ANDREW
%Pos_z and Pos_x at the midpoint?
% Making sure all 4 corners of the QUS ROIs are within the repro ROI
% all_roi_z2 = all_roi_z + all_roi.len_z;
% all_roi_x2 = all_roi_x + all_roi.len_x;
% [xm2,zm2] = meshgrid(all_roi_x2,all_roi_z2);
% [xm3,zm3] = meshgrid(all_roi_x,all_roi_z2);
% [xm4,zm4] = meshgrid(all_roi_x2,all_roi_z);
% in_roi_idx12 = inpolygon(xm2(:),zm2(:),roi_poly(:,1),roi_poly(:,2));
% in_roi_idx21 = inpolygon(xm3(:),zm3(:),roi_poly(:,1),roi_poly(:,2));
% in_roi_idx22 = inpolygon(xm4(:),zm4(:),roi_poly(:,1),roi_poly(:,2));
% 
% in_roi_idx = in_roi_idx | in_roi_idx12 | in_roi_idx21 | in_roi_idx22;

in_roi_idx = find(in_roi_idx);
% in_roi_idx12 = find(in_roi_idx12);
% in_roi_idx21 = find(in_roi_idx21);
% in_roi_idx22 = find(in_roi_idx22);

% Check for any ROI that has the ESD as an imaginary number. For these,
% all QUS parameters should be ignored (or set to zero, which
% effectively ignores it below)
esd_map = result{param2idx(params,'Effective Scatterer Size')};
esd_map = real(esd_map);
zeros_map = (esd_map <= 0) | (isnan(esd_map));


% For reference, plot the repro segmentation on top of the B-mode image
figure(5); clf;
set(gcf,'Units','Normalized','Position',[0.2, 0.25, 0.7, 0.5]);
imagesc(lateral_vec,axial_vec,this_bmode,[-60,0]); 
colormap gray; axis image;
xlabel('Lateral (mm)');
ylabel('Axial (mm)');
ylim([axial_vec(1),axial_vec(end)])
hold on;
repro.plot_roi(gca,roi_poly*1000);
set(gca,'FontSize',16);

%%ADDED BY ANDREW
%Add horizontal lines to indicate limits of inclusion
if AddLines
    plot(gca,[0, max(lateral_vec)], [12,12],'y','LineWidth',3);
    plot(gca,[0, max(lateral_vec)], [13,13],'g','LineWidth',3);
    plot(gca,[0, max(lateral_vec)], [14,14],'b','LineWidth',3);
    plot(gca,[0, max(lateral_vec)], [15,15],'m','LineWidth',3);
    plot(gca,[0, max(lateral_vec)], [16,16],'r','LineWidth',3);
end

save_str = sprintf('%s\\bmode_with_segmentation',copy_dir);
saveas(gcf,save_str,'png');
saveas(gcf,save_str,'fig');

% For each specified QUS parameter, compute the mean and std. And keep
% track of the number of ROIs used in the calculation
this_qus_cell = cell(1,length(qus_params_to_plot));

all_res_table = table;

for p_count=1:length(qus_params_to_plot)
    
%     scatter(xm(in_roi_idx).*1000,zm(in_roi_idx).*1000,'.r');

    figure(5+p_count); clf;
    set(gcf,'Units','Normalized','Position',[0.2, 0.25, 0.7, 0.5]);
    imagesc(lateral_vec,axial_vec,this_bmode,[-60,0]); 
    colormap gray; axis image;
    xlabel('Lateral (mm)');
    ylabel('Axial (mm)');
    xlabel('');
    ylim([axial_vec(1),axial_vec(end)])

    title_str = qus_str{p_count};
    title(title_str);

    len_z = all_roi.len_z*1000;
    len_x = all_roi.len_x*1000;
    color_list = lines(4);
    overlap_x = all_roi.o_x; 
    overlap_z = all_roi.o_z;
    
    
    this_p_idx = param2idx(params,qus_params_to_plot{p_count});
    this_p_map = result{param2idx(params,qus_params_to_plot{p_count})}; % extract the QUS parameter map
    this_p_map = real(this_p_map);
    
    this_p_map(zeros_map) = 0;
    
    this_alpha_map = zeros(size(this_p_map));
    %this_alpha_map(in_roi_idx) = 0.4;
    this_alpha_map(in_roi_idx) = 0.5;
    
    this_alpha_map(this_p_map == 0) = 0;
    this_alpha_map(isnan(this_p_map)) = 0;
    
    if strcmp(qus_params_to_plot{p_count},'Effective Scatterer Size')
        this_p_map = this_p_map * 1e6;
    elseif strcmp(qus_params_to_plot{p_count},'Nak Scale Factor')
        this_p_map = log10(this_p_map);
    elseif strcmp(qus_params_to_plot{p_count},'HK Scatterer Clustering Param')
        this_p_map = log10(this_p_map);
    elseif strcmp(qus_params_to_plot{p_count},'Burr_b')
        this_p_map = log10(this_p_map);
    end
    
    figure(5+p_count);
    sub_ax = gca;
    set(gca,'FontSize',14);
    colorbar('Location','southoutside','Visible','off');
    sub_pos = get(gca,'Position');
    new_ax = axes;
    imagesc(new_ax,(all_roi_x)*1000,(all_roi_z)*1000,...
        this_p_map,'AlphaData',this_alpha_map);
    axis(new_ax,'image');
    rg_cmap = interp1([0;1],[1 0 0; 0 1 0],linspace(0,1,256));
%         this_cmamp = interp1([0;1],qus_cmap{p_count},linspace(0,1,256));
    colormap(new_ax,isolum(256,'invert',false));
    new_ax.Visible = false;
    this_cbar = colorbar('Location','southoutside');
    set(gca,'FontSize',14);
    caxis(qus_clims{p_count});

    repro.plot_roi(new_ax,roi_poly*1000);

    %%ADDED BY ANDREW
    %Add horizontal lines to indicate limits of inclusion
    if AddLines
        plot(gca,[0, max(lateral_vec)], [12,12],'y','LineWidth',3);
        plot(gca,[0, max(lateral_vec)], [13,13],'g','LineWidth',3);
        plot(gca,[0, max(lateral_vec)], [14,14],'b','LineWidth',3);
        plot(gca,[0, max(lateral_vec)], [15,15],'m','LineWidth',3);
        plot(gca,[0, max(lateral_vec)], [16,16],'r','LineWidth',3);
    end
    
    new_ax.XLim = sub_ax.XLim;
    new_ax.YLim = sub_ax.YLim;
    new_ax.Position = sub_ax.Position;
    xlabel(new_ax,'');
    
    this_p_map = this_p_map(in_roi_idx);

    % Add the results for each ROI to the results table
    all_res_table.(qus_params_to_plot{p_count}) = this_p_map(:);
    
    % Any pixel with value 0 or NaN should be removed
    this_p_map = this_p_map(:);
    this_nan = isnan(this_p_map);
    
    this_remove_N = 0;
    this_remove_N = this_remove_N + sum(this_nan);
    this_p_map(this_nan) = [];
    
    this_zero = find(this_p_map == 0);
    this_remove_N = this_remove_N + length(this_zero);
    this_p_map(this_zero) = [];

    % Remove any Inf values (particular issue for log10(naka omega)
    this_p_map = this_p_map(~isinf(this_p_map));
    
    qus_mean(this_p_idx) = mean(this_p_map,'omitnan');
    qus_std(this_p_idx) = std(this_p_map);
    qus_roi_N(this_p_idx) = numel(this_p_map);
    qus_roi_zero(this_p_idx) = this_remove_N;
    
    this_qus_cell{p_count} = qus_mean(this_p_idx);
    
    save_str = sprintf('%s\\bmode_with_qus_overlay_%s',copy_dir,qus_str{p_count});
    saveas(gcf,save_str,'png');
    saveas(gcf,save_str,'fig');
% 
%     copy_fid = sprintf('%s\\%d_%d_%s',copy_dir,this_mrn,this_dataset,qus_str{p_count});
%     copyfile([save_str,'.png'],[copy_fid,'.png']);
%     copyfile([save_str,'.fig'],[copy_fid,'.fig']);
    
end



qus_row_count = qus_row_count + 1;
qus_result_cell(qus_row_count,:) = this_qus_cell;

% all_res_table.Day = ones(size(all_res_table,1),1)*day;
 end
% save('.\qus_results_stats_2023-1-13','qus_mean','qus_std','qus_roi_N','qus_roi_zero', ...
%     'qus_params_to_plot','all_res_table');

%     save('.\qus_result_stats.mat','qus_mean','qus_std','qus_roi_N','qus_roi_zero', ...
%         'params','qus_params','this_mrn','this_dataset','this_node_sus');

    


% header_cell = [{'MRN','Dataset','Score'},qus_params];

% full_qus_cell = [header_cell; qus_result_cell];