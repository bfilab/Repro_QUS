%%
eye_fid = dir('*.eye');
eye_fid = eye_fid(1).name;

eye_sig = eye.EyeData;
eye_sig.readEye(eye_fid);

fs = eye_sig.settings.fs*1e6;
c = 1460;
del_mm = eye_sig.settings.del;

pass_band = [0.3,3]*1e6;
[b,a] = cheby2(5,80,pass_band/(fs/2));

rf_data = eye_sig.data;
% rf_filt = rf_data;
rf_filt = filtfilt(b,a,rf_data);

axial_vec = [0:size(rf_data,1)-1]*c/(2*fs)*1000+del_mm;
dx = 100e-3;
lateral_vec = [0:size(rf_data,2)-1]*dx;

figure(40); clf;
plot(axial_vec,rf_data(:,75));
figure(41); clf;
plot(axial_vec,rf_filt(:,75));


bmode_im = 20*log10(abs(hilbert(rf_filt)));
bmode_im = bmode_im - max(bmode_im(:));
figure(42); clf;
imagesc(lateral_vec,axial_vec,bmode_im); colormap gray;

in_data = [];
in_data.rf_data = rf_filt;
in_data.fs = fs;