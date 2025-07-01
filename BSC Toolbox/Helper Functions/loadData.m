function [samp] = loadData(samp_dir,samp_fname,type)
%LOADDATA Load RF imaging data into a struct used with the BSC toolbox.

%The struct is filled with the following data:
% (1) 3D RF data
% (2) 3D envelope data
% (3) Vectors for dimensions of RF data, in meters
% (4) Image resolution for each dimension, in meters/voxel
% (5) Speed of sound (used to convert between time of flight and distance)
% (6) Time of flight vector, seconds
% (7) Axial sampling frequency, Hz
% (8) Surface position, in meters

% Because data can be stored in different formats, the different cases of
% the switch-case correspond to different datasets.
% (1) prostate_invivo - P###_invivoData.mat, from Duke University; data
%                       collected with the ER7B and 12L4 transducer
% (2) 2Dphantom - 2D phantom data that was collected with the ER7B and 12L4
%                 transducers, used in conjunction with the prostate_invivo
%                 data.
% (3) beamform - beamformed RF data from lung project using Mark's code.
% (4) beamform2 - beamformed RF data from lung project using MUST.
% (5) etc...

% INPUTS:
%   samp_dir = directory containing imaging data .mat file
%   samp_fname = filename of .mat file containing imaging data
%   type = format of data (string)
% OUTPUTS:
%   samp =  struct containing imaging data and associated parameters

% 06/24/2020 (THL): Created
% 08/25/2020 (THL): Added option for beamformed RF data
% 08/31/2020 (THL): Added option for beamformed RF data using MUST.
% 09/10/2020 (THL): Removed previous beamformed RF options; added option
%                   for MUST IQ data.
% 09/21/2020 (THL): Added SUNY option (not complete, missing surface
%                   location and all planes; note this is the data
%                   formatted into the attenuation estimation format)
% 11/17/2020 (THL): Added IRM option for lung project.
% 02/17/2021 (THL): Added LiverRF for liver project.

%% Load data
samp_l = load(fullfile(samp_dir,samp_fname));


%% Fill struct fields
switch type
    
    case 'prostate_invivo'
        
        % Convert to double
        samp_l.rfData = double(samp_l.rfData);
        samp_l.lat = double(samp_l.lat);
        samp_l.axial = double(samp_l.axial);
        
        % Get rf data
        samp.data = samp_l.rfData;
        
        % Get envelope data
        samp.env = abs(hilbert(samp.data));
        
        % Set dimension vectors (in m)
        samp.x = samp_l.lat*1e-3;
        samp.y = [1:length(samp_l.eleAngle)]*1e-3;
        samp.z = (samp_l.axial(1):mean(diff(samp_l.axial)):samp_l.axial(end))*1e-3;
        
        % Set image resolution (m/voxel)
        samp.dx = mean(diff(samp.x));
        samp.dy = mean(diff(samp.y));
        samp.dz = mean(diff(samp.z));
        
        % Set sampling frequency & time vector
        samp.c = 1500;
        samp.t = dist2time(samp.c,samp.z);
        samp.fs = 1/(mean(diff(samp.t)));
        
        % Surface position
        samp.surf_pos = samp.z(1);
        
    case '2Dphantom'
        
        % Get rf data
        samp.data = samp_l.data;
        
        % Get envelope data
        samp.env = abs(hilbert(samp.data));
        
        % Set dimension vectors (in m)
        samp.x = samp_l.lateral_mm*1e-3;
        samp.y = 1;
        samp.z = samp_l.axial_mm(1):mean(diff(samp_l.axial_mm)):samp_l.axial_mm(end);
        samp.z = samp.z*1e-3;
        
        % Set image resolution (m/voxel)
        samp.dx = mean(diff(samp.x));
        samp.dy = 1;
        samp.dz = mean(diff(samp.z));
        
        % Set sampling frequency
        samp.c = samp_l.c;
        samp.t = dist2time(samp.c,samp.z);
        samp.fs = 1/(mean(diff(samp.t)));
        
        % Surface position
        samp.surf_pos = samp.z(1);
        
    case '2DphantomCH'
        % Used when the calibration data file was created with the script 
        % create_QUS_refPhantom_calibration_file.m
        
        % Get 3D RF data
        samp.data = samp_l.rf_data;
        
        % Get envelope data
        samp.env = abs(hilbert(samp.data));
        
        % Set dimension vectors (in m)
        samp.x = samp_l.lateral_vec/1000; % convert from mm to m
        del_y = mean(diff(samp.x)); % assume elevational step size is same as lateral
        samp.y = [0:size(samp.data,3)-1].*del_y; 
        samp.z = samp_l.axial_vec/1000; % convert from mm to m
        
        % Set image resolution (m/voxel)
        samp.dx = mean(diff(samp.x));
        samp.dy = 1; %mean(diff(samp.y));
        samp.dz = mean(diff(samp.z));
        
        if length(samp.y) == 1
            samp.y = 1;
            samp.dy = 1;
        end
        
        % Set sampling frequency
        samp.c = samp_l.c;
        samp.t = dist2time(samp.c,samp.z);
        samp.fs = 1/(mean(diff(samp.t)));
        
        % Surface position
        samp.surf_pos = samp_l.surf_location/1000;
        
        
    case 'MUSTIQ'
        
        % Setup
        n = 4;
        param.fs = 62.5*1e6;
        param.pitch = 0.3*1e-3;
        param.fc = 7.8e6;
        param.c = 1540;
        param.fnumber = 1;
        N = size(samp_l.full_BFIQ(:,:,1),1)*n;
        
        % Get rf data
        samp.data = iq2rf(samp_l.full_BFIQ(:,:,1),param.fs/n,param.fc,N);
        
        % Get envelope data
        samp.env = abs(hilbert(samp.data));
        
        % Set dimension vectors (in m)
        samp.x = samp_l.x;
        samp.y = 1;
        samp.z = samp_l.z;
        
        % Set image resolution (m/voxel)
        samp.dx = mean(diff(samp.x));
        samp.dy = 1;
        samp.dz = mean(diff(samp.z));
        
        % Set sampling frequency
        samp.c = param.c;
        samp.t = dist2time(samp.c,samp.z);
        % samp.fs = 1/mean(diff(samp.t));
        samp.fs = param.fs;
        
        % Surface position
        samp.surf_pos = samp.z(1);
        
    case 'MUSTRF'
        
        % Setup
        param.fs = 62.5*1e6;
        param.pitch = 0.3*1e-3;
        param.fc = 7.8e6;
        param.c = 1540;
        param.fnumber = 1;
        
        % Get rf data
        samp.data = samp_l.full_BFRF;
        
        % Get envelope data
        samp.env = abs(hilbert(samp.data));
        
        % Set dimension vectors (in m)
        samp.x = samp_l.x;
        samp.y = 1:size(samp.data,3);
        samp.z = samp_l.z;
        
        % Set image resolution (m/voxel)
        samp.dx = mean(diff(samp.x));
        samp.dy = 1;
        samp.dz = mean(diff(samp.z));
        
        % Set sampling frequency
        samp.c = param.c;
        samp.t = dist2time(samp.c,samp.z);
        samp.fs = 1/mean(diff(samp.t));
        % samp.fs = param.fs;
        
        % Surface position
        samp.surf_pos = samp.z(1);
        
    case 'SUNY'
        
        % Get rf data
        samp.data = samp_l.data(:,:,1);
        
        % Get envelope data
        samp.env = abs(hilbert(samp.data));
        
        % Set dimension vectors (in m)
        samp.x = samp_l.lateral_mm*1e-3;
        samp.y = 1;
        samp.z = samp_l.axial_mm*1e-3;
        
        % Set image resolution (m/voxel)
        samp.dx = mean(diff(samp.x));
        samp.dy = 1;
        samp.dz = mean(diff(samp.z));
        
        % Set sampling frequency
        samp.c = samp_l.c;
        samp.t = dist2time(samp.c,samp.z);
        samp.fs = 1/mean(diff(samp.t));
        
        % Surface position
        samp.surf_pos = samp.z(1);
      
        
    case 'IRM'
        
        % Setup
        param.fs = 62.5*1e6;
        param.pitch = 0.3*1e-3;
        param.fc = 7.8e6;
        param.c = 1540;
        param.fnumber = 1;
        
        dz = (1/(2*param.fs))*param.c;
       
        % Get rf data
        samp.data = samp_l.Rf;
        
        % Get envelope data
        samp.env = abs(hilbert(samp.data));
        
        % Set dimension vectors (in m)
        samp.x = param.pitch*[0:(size(samp.data,2)-1)];
        samp.y = param.pitch*[0:(size(samp.data,3)-1)];
        samp.z = dz*[0:(size(samp.data,1)-1)];
        
        % Set image resolution (m/voxel)
        samp.dx = mean(diff(samp.x));
        samp.dy = mean(diff(samp.y));
        samp.dz = mean(diff(samp.z));
        
        % Set sampling frequency
        samp.c = param.c;
        samp.t = dist2time(samp.c,samp.z);
        samp.fs = 1/mean(diff(samp.t));
        % samp.fs = param.fs;
        
        % Surface position
        samp.surf_pos = samp.z(1);
      
    case 'PigEye'
        
        c = 1540;
        
        % Get RF data
        samp.data = samp_l.full_RF;
        
        % Get envelope data
        samp.env = abs(hilbert(samp.data));
        
        % Set axes vectors (in m)
        samp.x = samp_l.lateral_vec/1000; % convert from mm to m
        samp.z = samp_l.axial_vec/1000; % convert from mm to m
        samp.y = [0:(size(samp.data,3)-1)]*samp_l.delta_y/1000;
        
        % Set image resolution (m/voxel)
        samp.dx = mean(diff(samp.x));
        samp.dy = mean(diff(samp.y));
        samp.dz = mean(diff(samp.z));
        
        % Set sampling frequency
        samp.c = c;
        samp.t = dist2time(samp.c,samp.z);
        samp.fs = 1/mean(diff(samp.t));
        
        samp.fs = samp_l.fs;
        
        if length(samp.y) == 1
            samp.y = 1;
            samp.dy = 1;
        end
        
    case 'LiverRF'
        
        % Setup
        c = 1540;
        param = samp_l.param;
        ss = param.stepsize*1e-3;
        fs = samp_l.fs_int;
        dt = 1/fs;
        dz = (dt/2)*c;
        
        % Get rf data
        samp.data = samp_l.full_RF;
        
        % Get envelope data
        samp.env = abs(hilbert(samp.data));
        
        % Set dimension vectors (in m)
        samp.x = linspace(0,param.BmodeWidth*1e-3,size(samp.data,2));
        samp.y = 0:ss:((size(samp.data,3)-1)*ss);
        samp.z = 0:dz:((size(samp.data,1)-1)*dz);
        samp.z = samp.z + param.BmodeDepthOffset*1e-3;
        
        % Set image resolution (m/voxel)
        samp.dx = mean(diff(samp.x));
        samp.dy = mean(diff(samp.y));
        samp.dz = mean(diff(samp.z));
        
        % Set sampling frequency
        samp.c = c;
        samp.t = dist2time(samp.c,samp.z);
        samp.fs = 1/mean(diff(samp.t));
        
        % Surface position
        % samp.surf_pos = samp.z(1);
        
    case 'SUNY_GE'
        
        % For this case, samp_1 will contain the file name of the GE
        % .bin file. Use other functions to load the echo data and set the
        % parameters
        % 
        % An optional parameter is surf_pos, which indicates the location
        % fo the target surface from the transducer. If non-existent or
        % empty, assume the transducer is in contant with the surface.
        
        ge_struct = readIQRFDataFrame(fullfile(samp_dir,samp_l.ge_fname),'RF');
        c = 1540;
        fs = ge_struct.samplingRate;
        num_pts = ge_struct.nSamplesPerBeam;
        
        delta_x = 0.1/1000; % line spacing in m, assuming 50um pitch
        
        % Get RF data
        samp.data = ge_struct.data{1};
        
        % Get envelope data
        samp.env = abs(hilbert(samp.data));
        
        samp.x = [0:size(samp.data,2)-1]*delta_x; % lateral axes vector in mm
        samp.y = 1; % single frame
        samp.z = [0:num_pts-1].*c/(2*fs); % axial axes vector in m
        
        % Set image resolution (m/voxel)
        samp.dx = mean(diff(samp.x));
        samp.dy = 1;
        samp.dz = mean(diff(samp.z));
        
        % Set sampling frequency
        samp.c = c;
        samp.t = dist2time(samp.c,samp.z);
        samp.fs = 1/mean(diff(samp.t));
        
        % Set the surface position
        if isfield(samp_l,'surf_pos') && ~isempty(samp_l.surf_pos)
            samp.surf_pos = samp_l.surf_pos;
        else
            samp.surf_pos = samp.z(1);
        end
        
%         keyboard;

    case 'Tulane_Vevo2100'

        % Generate the image axes
        [axial_vec, lateral_vec,c,elev_vec] = repro.generate_image_axes(samp_l.rf_data,samp_l.sysParam,samp_l.fs*1e6);

        % Get RF data
        samp.data = samp_l.rf_data;
        
        % Get envelope data
        samp.env = abs(hilbert(samp.data));
        
        samp.x = lateral_vec/1000; % lateral axes vector in m
        samp.y = elev_vec;
        samp.z = axial_vec/1000; % axial axes vector in m
        
        % Set image resolution (m/voxel)
        samp.dx = mean(diff(samp.x));
        samp.dz = mean(diff(samp.z));
        % Only set an elevational step size if more than 1 frame is present
        if length(samp.y) > 1
            samp.dy = mean(diff(samp.y));
        else
            samp.dy = 1;
        end

        % Set sampling frequency
        samp.c = c;
        samp.t = dist2time(samp.c,samp.z);
        samp.fs = samp_l.fs*1e6; % fs in the RF data struct is in [MHz], needs to be in [Hz]
        % samp.fs = 1/mean(diff(samp.t));
        
        % Set the surface position
        % This is specified in the filename for reference data
        % For sample data, this is a field that will have to be updated
        % later based on the segmentation
        surf_pos = regexpi(samp_fname,'(?<=_)\d{1,2}(?=mmSurf)','match','once');
        if ~isempty(surf_pos)
            surf_pos = str2double(surf_pos);
            samp.surf_pos = surf_pos/1000; % convert to [m]
        else
            samp.surf_pos = [];
        end
        
    otherwise
        
        disp('Unknown data format.');
        
end

end

