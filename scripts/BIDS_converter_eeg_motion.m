%% Convert source data to BIDS
clc; clear; close all;

%% Load required paths and libraries

mainpath = 'C:\Users\koes4731\Desktop\eeglab2023.0\'; % eeglab folder
addpath('C:\Users\koes4731\Desktop\fieldtrip-20230215\'); % add Fieldtrip
path = 'C:\Users\koes4731\Desktop\Thesis\MediaPipe\'; % raw data path
outpath = 'C:\Users\koes4731\Desktop\Thesis\data_hoops\'; % BIDS output path
files = dir(fullfile(path, '\*.xdf')); % listing datasets

ft_defaults; % Fieldtrip defaults

%% Set general BIDS conversion configuration

cfg = [];
cfg.InstitutionName = 'University of Oldenburg';
cfg.dataset_description.Name = 'EEG_Basketball';
cfg.dataset_description.BIDSVersion = '1.9';
cfg.method = 'convert'; % The original data is in a BIDS-compliant format and can be copied
cfg.bidsroot = './data/bids';  % write to the present working directory

%% Loop over datasets

for sub = 1: 1%length(files)

    participant = extractBefore(files(sub).name, '.xdf');  % get subject name
    data = load_xdf([path, files(sub).name]); % Saving the data in a variable
    out_subfold = [outpath, participant, '\\'];
    load([out_subfold, 'events_', participant,'.mat']); % Loading events file


    %% Extract EEG and motion data from XDF file
    for i = 1:length(data)
        currentName = data{1, i}.info.name;

        if contains(currentName, 'Pose', 'IgnoreCase', true)
            mp = data{1, i}; % Pose (motion) data
        elseif contains(currentName, 'Android_EEG', 'IgnoreCase', true)
            eeg = data{1, i}; % EEG data
        end
    end

    %% EEG Data Conversion

    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab; % Open EEGLAB

    % Import EEG, add channel locations, and add events
    EEG = pop_loadxdf([path, files(sub).name], 'streamtype', 'EEG', 'exclude_markerstreams', {});
    EEG = pop_chanedit(EEG, 'lookup', [mainpath, 'plugins/dipfit/standard_BEM/elec/standard_1005.elc']); % Add channel info
    EEG.event = events; % Add events
    EEG = eeg_checkset(EEG, 'eventconsistency'); % Check event consistency
    eeglab redraw; % Update GUI

    % Channel labels and electrode info
    eeg_channel_labels = {EEG.chanlocs.labels};
    elec_file = [mainpath, 'plugins/dipfit/standard_BEM/elec/standard_1005.elc'];
    elec = ft_read_sens(elec_file); % Load electrode file
    cfg.elec = elec;
    
    % Define Fieldtrip EEG structure (eeglab2fieldtrip)
    data_eeg = eeglab2fieldtrip( EEG, 'raw', 'none');
    eeg = data_eeg;

    % % Define Fieldtrip EEG structure
    % ft_data_eeg = [];
    % ft_data_eeg.trial{1} = eeg.time_series; % EEG data (channels x timepoints)
    % ft_data_eeg.time{1} = eeg.time_stamps; % Time points for data
    % ft_data_eeg.label = eeg_channel_labels; % 42 EEG channel labels
    % ft_data_eeg.elec = elec; % Electrode structure for the 42 channels

    % BIDS-specific settings for EEG
    cfg.sub = extractAfter(participant, 'sub_');
    cfg.datatype = 'eeg';
    ft_checkdata(eeg);
    cfg.task = 'Freethrow';

    cfg.eeg.PowerLineFrequency = 50; % Power line frequency (50Hz for EU)
    cfg.eeg.EEGReference = 'TP9 TP10';  % Reference electrodes
    cfg.eeg.InstitutionName = 'University of Oldenburg';
    cfg.eeg.InstitutionAddress = 'Ammerlaender Heerstr. 114-118, 26129 Oldenburg, Germany';
    cfg.eeg.ManufacturersModelName = 'mbraintrain';

    % Generate dataset_description.json
    cfg.dataset_description.Name = 'Readiness Potential in Basketball';
    cfg.dataset_description.BIDSVersion = 'v1.9';
    cfg.dataset_description.Authors = {'Miguel Contreras-Altamirano', 'Stefan Debener'};
    cfg.dataset_description.License = 'Creative Commons';
    cfg.dataset_description.DatasetType = 'raw';

    % Generate README file
    README = sprintf('The experiment included 27 participants. \n- Miguel Contreras-Altamirano (December, 2024)');

    % Call data2bids for EEG
    data2bids(cfg, data_eeg);


    %% Generate EEG coordsystem.json
    
% % Ensure the 'eeg' subfolder exists before trying to write the file
% eeg_folder = fullfile(out_subfold, 'eeg');
% if ~exist(eeg_folder, 'dir')
%     mkdir(eeg_folder);  % Create the 'eeg' folder if it doesn't exist
% end
% 
% % Generate the coordsystem.json file for EEG data
% coordsystem_json = fullfile(eeg_folder, [participant '_coordsystem.json']);
% fid = fopen(coordsystem_json, 'w');
% 
% % Ensure the file opened successfully
% if fid == -1
%     error('Could not open file for writing: %s', coordsystem_json);
% end
% 
% % Write the required keys into the coordsystem.json
% fprintf(fid, '{\n');
% fprintf(fid, '    "EEGCoordinateSystem": "CTF",\n'); % Replace "CTF" with the correct system if applicable
% fprintf(fid, '    "EEGCoordinateUnits": "mm",\n');  % Units for the EEG coordinates
% fprintf(fid, '    "Manufacturer": "mbraintrain",\n'); % Manufacturer of the EEG device
% fprintf(fid, '    "ManufacturersModelName": "mbraintrain",\n'); % Model name of the EEG system
% fprintf(fid, '    "CapManufacturersModelName": "mbraintrain cap",\n'); % Model of the EEG cap
% fprintf(fid, '    "EEGCoordinateSystemDescription": "Standard 1005",\n'); % Description of the EEG system
% fprintf(fid, '    "EEGCoordinateSystemReference": "https://doi.org/10.1371/journal.pbio.1002295"\n'); % Reference link for the coordinate system
% fprintf(fid, '}');
% fclose(fid);


    %% Motion Data (MediaPipe) Conversion

    % Prepare motion data structure
    mocap.trial{1} = [mp.time_series];
    mocap.time{1} = mp.time_stamps;
    mocap.label = {
        'nose_X', 'nose_Y', 'nose_Z', ...
        'left eye (inner)_X', 'left eye (inner)_Y', 'left eye (inner)_Z', ...
        'left eye_X', 'left eye_Y', 'left eye_Z', ...
        'left eye (outer)_X', 'left eye (outer)_Y', 'left eye (outer)_Z', ...
        'right eye (inner)_X', 'right eye (inner)_Y', 'right eye (inner)_Z', ...
        'right eye_X', 'right eye_Y', 'right eye_Z', ...
        'right eye (outer)_X', 'right eye (outer)_Y', 'right eye (outer)_Z', ...
        'left ear_X', 'left ear_Y', 'left ear_Z', ...
        'right ear_X', 'right ear_Y', 'right ear_Z', ...
        'mouth (left)_X', 'mouth (left)_Y', 'mouth (left)_Z', ...
        'mouth (right)_X', 'mouth (right)_Y', 'mouth (right)_Z', ...
        'left shoulder_X', 'left shoulder_Y', 'left shoulder_Z', ...
        'right shoulder_X', 'right shoulder_Y', 'right shoulder_Z', ...
        'left elbow_X', 'left elbow_Y', 'left elbow_Z', ...
        'right elbow_X', 'right elbow_Y', 'right elbow_Z', ...
        'left wrist_X', 'left wrist_Y', 'left wrist_Z', ...
        'right wrist_X', 'right wrist_Y', 'right wrist_Z', ...
        'left pinky_X', 'left pinky_Y', 'left pinky_Z', ...
        'right pinky_X', 'right pinky_Y', 'right pinky_Z', ...
        'left index_X', 'left index_Y', 'left index_Z', ...
        'right index_X', 'right index_Y', 'right index_Z', ...
        'left thumb_X', 'left thumb_Y', 'left thumb_Z', ...
        'right thumb_X', 'right thumb_Y', 'right thumb_Z', ...
        'left hip_X', 'left hip_Y', 'left hip_Z', ...
        'right hip_X', 'right hip_Y', 'right hip_Z', ...
        'left knee_X', 'left knee_Y', 'left knee_Z', ...
        'right knee_X', 'right knee_Y', 'right knee_Z', ...
        'left ankle_X', 'left ankle_Y', 'left ankle_Z', ...
        'right ankle_X', 'right ankle_Y', 'right ankle_Z', ...
        'left heel_X', 'left heel_Y', 'left heel_Z', ...
        'right heel_X', 'right heel_Y', 'right heel_Z', ...
        'left foot index_X', 'left foot index_Y', 'left foot index_Z', ...
        'right foot index_X', 'right foot index_Y', 'right foot index_Z'
        };

    % BIDS motion data settings
    cfg.tracksys = 'MediaPipe';
    cfg.motion.TrackingSystemName = 'MediaPipe';
    cfg.motion.samplingrate = sampling_rate_mp;

    % specify channel details, this overrides the details in the original data structure
    cfg.channels = [];
    cfg.channels.name = mocap.label;
    cfg.channels.tracked_point = mocap.label;
    cfg.channels.type = cellstr(repmat('POS',length(mocap.label),1));
    cfg.channels.units = cellstr(repmat('m',length(mocap.label),1));

    mocap = ft_datatype_raw(mocap);

    cfg.datatype = 'motion';
    data2bids(cfg, mocap);

end

