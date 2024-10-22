clc, clear, close all;

%% Settings data

mainpath = 'C:\Users\micua\Desktop\eeglab2023.0\'; % eeglab folder
path = 'C:\Users\micua\OneDrive - Benemérita Universidad Autónoma de Puebla\NCP_Basketball\MediaPipe\';
outpath = 'C:\\Users\\micua\\OneDrive - Benemérita Universidad Autónoma de Puebla\\Oldenburg_University\\Thesis\\data_hoops\\';
files = dir( fullfile( path,'\*.xdf')); % listing data sets

channels_RP_bottom = {'C3','Cz','C4'};  % Channels for bottom row (central)
channels_RP_top = {'FC1','Fz','FC2'};   % Channels for top row (frontal)

%% Process for both sets of channels

% Initialize cell arrays to hold matrices for both sets of channels
R_Squared_Matrices_RP_bottom = cell(1, length(channels_RP_bottom));
R_Squared_Matrices_RP_top = cell(1, length(channels_RP_top));

% This variable will hold feature names (assuming they are the same across all channels)
feature_names = {};

% Process Readiness Potential for both sets of channels
for ch = 1:length(channels_RP_bottom)
    [R_Squared_Matrices_RP_bottom{ch}, feature_names] = process_channels(files, channels_RP_bottom{ch}, outpath, 'RP');
end

for ch = 1:length(channels_RP_top)
    [R_Squared_Matrices_RP_top{ch}, ~] = process_channels(files, channels_RP_top{ch}, outpath, 'RP');  % We already got the feature names
end

%% Visualization of features variance across participants

figure('units','normalized','outerposition', [0 0 1 1]);

% Plotting for top row (Frontal channels: FC1, Fz, FC2)
for ch = 1:length(channels_RP_top)
    subplot(2, 3, ch);
    imagesc(R_Squared_Matrices_RP_top{ch});
    colorbar;
    title(['Explained Variance of ', channels_RP_top{ch}, ' [R^2]'], 'FontSize', 12);
    subtitle('Trials prior movement / [Readiness Potential]', 'FontSize', 12);
    xlabel('Participants', 'FontWeight','bold', 'FontSize', 11.5);
    ylabel('Features', 'FontWeight','bold', 'FontSize', 11.5);
    axis tight;
    xticks(1:length(files) + 2);  % Number of participants + 2 for padding/mean
    xticklabels([arrayfun(@num2str, 1:length(files), 'UniformOutput', false), ' ', 'Mean']);
    yticks(1:length(feature_names) + 1);  % Number of features + mean
    yticklabels(feature_names);  % Use feature names for the y-axis labels
    colormap('jet');
    clim([0, max(R_Squared_Matrices_RP_top{ch}(:))]);
end

% Plotting for bottom row (Central channels: C3, Cz, C4)
for ch = 1:length(channels_RP_bottom)
    subplot(2, 3, ch+3);
    imagesc(R_Squared_Matrices_RP_bottom{ch});
    colorbar;
    title(['Explained Variance of ', channels_RP_bottom{ch}, ' [R^2]'], 'FontSize', 12);
    subtitle('Trials prior movement / [Readiness Potential]', 'FontSize', 12);
    xlabel('Participants', 'FontWeight','bold', 'FontSize', 11.5);
    ylabel('Features', 'FontWeight','bold', 'FontSize', 11.5);
    axis tight;
    xticks(1:length(files) + 2);  % Number of participants + 2 for padding/mean
    xticklabels([arrayfun(@num2str, 1:length(files), 'UniformOutput', false), ' ', 'Mean']);
    yticks(1:length(feature_names) + 1);  % Number of features + mean
    yticklabels(feature_names);  % Use feature names for the y-axis labels
    colormap('jet');
    clim([0, max(R_Squared_Matrices_RP_bottom{ch}(:))]);
end

% Save the figure
ROI = 'fronto-central';

% saveas(gcf, [outpath, 'Feat_variance_', ROI, '.jpg']);
saveas(gcf, [outpath, '\\group_analysis\\','Feat_variance_', ROI, '.jpg']);

%% Local function definition 

function [R_Squared_Matrix, feature_names] = process_channels(files, channel, outpath, condition)
    numParticipants = length(files); % number of participants
    numFeatures = 30; % 30 features

    % Initialize a matrix to hold R_Squared values for all features across all participants
    R_Squared_Matrix = zeros(numFeatures, numParticipants);
    feature_names = {};  % Initialize feature names

    for sub = 1:length(files)
        participant = extractBefore(files(sub).name, '.xdf');
        out_subfold = [outpath, participant, '\\'];
        load([out_subfold, 'Biserial_corr_', condition, '_', channel, '_', participant, '.mat']);

        % Fill in the matrix with R_Squared values for this participant
        R_Squared_Matrix(:, sub) = Biserial_corr_RP.R_Squared;

        % Save the feature names if they haven't been assigned yet
        if isempty(feature_names)
            feature_names = Biserial_corr_RP.Features;
        end
    end

    % Adding extra layer for visual separation
    R_Squared_Matrix(numFeatures+1 , ( size(R_Squared_Matrix, 2) + 1)) = NaN;

    % Calculate means
    R_Squared_Matrix(numFeatures+2, :) = mean(R_Squared_Matrix, 1, 'omitmissing');
    R_Squared_Matrix(:, size(R_Squared_Matrix, 2) + 1) = mean(R_Squared_Matrix, 2, 'omitmissing');
    feature_names{numFeatures+2} = 'Mean';
end
