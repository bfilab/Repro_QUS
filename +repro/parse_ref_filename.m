function [trans,ph,surf_pos,tx_pos] = parse_ref_filename(in_fid)
% Parse the filename of the reference data to extract the transducer
% number, calibration phantom, phantom surface location, and Tx focus
% location
%
%   Inputs
%       in_fid: file name of the reference data

% The reference datasets filenames follow the format
% LZTTT_PPumPhantom_SmmSurf_TmmTx, where
%   TTT: transducer number (250 or 550) - return as double
%   PP : calibration phantom (15 or 18) - return as double
%   S  : phantom surface position in [mm] - return as double
%   T  : Tx focus location in [mm] - return as double

trans = regexpi(in_fid,'(?<=LZ)\d{2,4}(?=_)','match','once');
ph = regexpi(in_fid,'(?<=_)\d{1,3}(?=umPhantom)','match','once');
surf_pos = regexpi(in_fid,'(?<=_)\d{1,2}(?=mmSurf)','match','once');
tx_pos = regexpi(in_fid,'(?<=_)\d{1,2}(?=mmTx)','match','once');

trans = str2double(trans);
ph = str2double(ph);
surf_pos = str2double(surf_pos);
tx_pos = str2double(tx_pos);