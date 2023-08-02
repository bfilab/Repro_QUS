function [final_avg_spec_db,final_avg_spec_noLog,f] = getPowerSpecWelch(data,fs,N_fft)
% GETPOWERSPECWELCH Use Welch's method for computing the power spectrum.
% INPUTS:
%   data: 3D matrix of rf data to be processed (dimension 1 should be axial)
%   fs: sampling frequency of data
%   N_fft: N for FFT
% OUTPUTS:
%   final_avg_spec_db: average spectrum for given rf data, processed with
%                      Welch's method
%   final_avg_spec_noLog: same as above, but in linear scale
%   f: frequency vector in Hz for spectrum

% Get size of data and appropriate window size
dim_axial = size(data,1);
dim_axial_win = floor(dim_axial/2);

% Compute spectrum in each window
all_avg_spec_noLog = zeros(N_fft/2,3);
for i = 1:3
    
    % Get window in depth
    win_idx1 = 1 + (i-1)*(floor(dim_axial_win/2));
    win_idx2 = win_idx1 + dim_axial_win - 1;
    
    % Window the data
    data_win = data(win_idx1:win_idx2,:,:);
    han_win = hanning(size(data_win,1));
    han_win = repmat(han_win,[1,size(data_win,2)]);
    han_win = repmat(han_win,[1,1,size(data_win,3)]);
    data_han_win = data_win.*han_win;
    
    % Compute the average spectra
    [f,~,~,spec_noLog] = getPowerSpec(data_han_win,fs,N_fft);
    avg_spec_noLog = mean(spec_noLog,2);
    avg_spec_noLog = mean(avg_spec_noLog,3);
    
    % Add to compilation
    all_avg_spec_noLog(:,i) = avg_spec_noLog;
    
end

% Average to get final spectrum
final_avg_spec_noLog = mean(all_avg_spec_noLog,2);
final_avg_spec_db = 10*log10(final_avg_spec_noLog);

end

