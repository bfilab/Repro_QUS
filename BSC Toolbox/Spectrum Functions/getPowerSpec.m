function [f,PS,Phase,PS_noLog] = getPowerSpec(y,Fs,N)
%GETPOWERSPEC Takes in 3D matrix of rf data, where time is along
%dimension 1, and outputs the spectra of each rf signal, where frequency is
%along dimension 1. Inputs and outputs are both 3D.
% INPUTS:
%	y = 3D matrix of RF signals, where time is along dimension 1.
%	Fs = sampling frequency [1/Hz]
%	N = N-point DFT
% OUTPUTS:
%   f = frequency [Hz]
%   PS = power spectrum in dB
%   Phase = phase
%   PS_noLog = power spectrum in linear scale

Y     = fft(y,N,1);
Phase = angle(Y(1:(N/2)+1,:, :));

% power spectral density
Pyy = abs(Y).^2;
f = Fs.*(0:(N/2-1))/N;

PS_noLog = Pyy(1:(N/2),:,:); 
PS       = 10*log10(PS_noLog);

end

