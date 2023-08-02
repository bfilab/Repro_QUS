%%
if exist('del_mm','var')
    clear del_mm
end

in_fid = uigetfile('*.mat');

in_mat = matfile(in_fid);

% in_mat.DepthOffset = 1;

del_mm = in_mat.DepthOffset;


%%
del_mm = 1;
save(in_fid,'del_mm','-append');

%%

for f_idx=1:size(in_mat,'RfData',3)
    in_data.rf_data = in_mat.RfData(:,:,f_idx);
    in_data.fs = in_mat.fs_int;
    in_data.del_mm = del_mm;
    
    if f_idx == 1
        [~,~,~,~,spec_roi] = utils.plot_spectrum_roi_arb(in_data);
    else
        utils.plot_spectrum_roi_arb(in_data,spec_roi);
    end

    input('');
end


%%

in_data.rf_data = in_mat.RfData(:,:,1);
in_data.fs = in_mat.fs_int;
in_data.del_mm = in_mat.DepthOffset;
% in_data.del_mm = 2;

[~,~,~,~,spec_roi] = utils.plot_spectrum_roi_arb(in_data,roi20);

%%

in_data.rf_data = RfData(:,:,1);
in_data.fs = fs_int;
in_data.del_mm = del_mm;

[~,~,~,~,spec_roi] = utils.plot_spectrum_roi_arb(in_data);