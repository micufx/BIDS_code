clc, clear, close all;

%% MediaPipe - Pose Landmark Detection (PLD)

% MediaPipe timestamps are in units of microseconds. However, the LSD send
% the data in seconds.

% Units:
% x and y : Landmark coordinates normalized between 0.0 and 1.0 by
% the image width ( x ) and height ( y ).
% z : The landmark depth, with the depth at the midpoint of the hips as
% the origin. The smaller the value, the closer the landmark is to the
% camera. The magnitude of z uses roughly the same scale as x .

% The pose landmarker model tracks 33 body landmark locations, representing
% the approximate location of body parts. The data is arranged by 99 rows,
% representing x,y and z coordinates of each body part in secuencial order,
% the columns represent the time series.

% The coordinate system is oriented with the positive Y-axis pointing
% downwards.

% Body points:
% 0 - nose
% 1 - left eye (inner)
% 2 - left eye
% 3 - left eye (outer)
% 4 - right eye (inner)
% 5 - right eye
% 6 - right eye (outer)
% 7 - left ear
% 8 - right ear
% 9 - mouth (left)
% 10 - mouth (right)
% 11 - left shoulder
% 12 - right shoulder
% 13 - left elbow
% 14 - right elbow
% 15 - left wrist
% 16 - right wrist
% 17 - left pinky
% 18 - right pinky
% 19 - left index
% 20 - right index
% 21 - left thumb
% 22 - right thumb
% 23 - left hip
% 24 - right hip
% 25 - left knee
% 26 - right knee
% 27 - left ankle
% 28 - right ankle
% 29 - left heel
% 30 - right heel
% 31 - left foot index
% 32 - right foot index

% Miguel Contreras-Altamirano, 2024


%% Loading xdf files

mainpath = 'C:\Users\micua\Desktop\eeglab2023.0\'; % eeglab folder
path = 'C:\Users\micua\OneDrive - Benemérita Universidad Autónoma de Puebla\NCP_Basketball\MediaPipe\'; % raw data
outpath = 'C:\\Users\\micua\\OneDrive - Benemérita Universidad Autónoma de Puebla\\Oldenburg_University\\Thesis\\data_hoops\\';
files = dir( fullfile( path,'\*.xdf')); % listing data sets


%% Selecting participant

for sub=1 : 1 %length(files)

    participant = extractBefore(files(sub).name, '.xdf');
    out_subfold = [outpath, participant, '\\'];
    data = load_xdf([path, files(sub).name]); % Saving the data in a variable


    % Finding the number struct of MP and EEG
    for i = 1:length(data)
        currentName = data{1, i}.info.name;

        % Check if the current data is MP
        if contains(currentName, 'Pose', 'IgnoreCase', true)
            mp = i;
        end

        % Check if the current data is EEG
        if contains(currentName, 'Android_EEG', 'IgnoreCase', true)
            eeg = i;
        end
    end

    % Assigning variables

    timeseries_mp = data{1, mp}.time_series; % MP time-series
    timeseries_eeg = data{1, eeg}.time_series; % EEG time-series

    timestamps_mp = data{1, mp}.time_stamps;  % MP timestamps
    timestamps_eeg = data{1, eeg} .time_stamps; % EEG timestamps


    % Inspecting data characteristics

    % EEG sampling rate
    srate_eeg = 1 / mean(diff(timestamps_eeg)); % calculates differences between adjacent elements of X and takes the mean divided by 1.
    timestamp_range_eeg = range(timestamps_eeg); %  returns the difference between the maximum and minimum values of sample data in X.
    disp(['Sampling rate EEG: ', num2str(srate_eeg)]);

    % Length of the recording
    length_minutes_eeg = timestamp_range_eeg / 60; % in minutes
    duration_seconds_eeg = str2double(data{1, eeg}.info.last_timestamp) - str2double(data{1, eeg}.info.first_timestamp);  % in seconds

    % Double check based on data info (more accurate)
    sampling_rate_eeg = str2double(data{1, eeg}.info.sample_count) / duration_seconds_eeg; % in seconds


    % PLD sampling rate
    srate_mp = 1 / mean(diff(timestamps_mp)); % calculates differences between adjacent elements of X and takes the mean divided by 1.
    timestamp_range_mp = range(timestamps_mp); %  returns the difference between the maximum and minimum values of sample data in X.
    disp(['Sampling rate PLD: ', num2str(srate_mp)]);

    % Length of the recording
    length_minutes_mp = timestamp_range_mp / 60; % in minutes
    duration_seconds_mp = str2double(data{1, mp}.info.last_timestamp) - str2double(data{1, mp}.info.first_timestamp); %  in seconds

    % Double check based on data info (more accurate)
    sampling_rate_mp = str2double(data{1, mp}.info.sample_count) / duration_seconds_mp; % in seconds


    % Display the results
    disp(['Sampling rate EEG: ', num2str(sampling_rate_eeg), ' samples per second']);
    disp(['Sampling rate PLD: ', num2str(sampling_rate_mp), ' samples per second']);
    disp(['PLD recording: ', num2str(length_minutes_mp), ' minutes']);
    disp(['EEG recording: ', num2str(length_minutes_eeg), ' minutes']);


    % % Interpolation
    % timeseries_mp = interp1(timestamps_mp, timeseries_mp', timestamps_eeg, 'linear')';
    % timestamps_mp = timestamps_eeg;

    % % Limit EEG timestamps to the range of MP timestamps
    % valid_indices = timestamps_eeg >= min(timestamps_mp) & timestamps_eeg <= max(timestamps_mp);
    % timestamps_eeg_valid = timestamps_eeg(valid_indices);
    %
    % % Interpolation only for the valid range
    % timeseries_mp = interp1(timestamps_mp, timeseries_mp', timestamps_eeg_valid, 'linear')';
    % timestamps_mp = timestamps_eeg;

    % Enable extrapolation in interp1
    timeseries_mp = interp1(timestamps_mp, timeseries_mp', timestamps_eeg, 'linear', 'extrap')';
    timestamps_mp = timestamps_eeg;


    %% Defining onset events based on Y coordinates interaction (right eye - right wrist)

    % Determine the number of frames and landmarks
    numFrames = size(timeseries_mp, 2); % Number of frames
    numLandmarks = size(timeseries_mp, 1) / 3; % Number of landmarks (33 according to MediaPipe)

    % Define the indices for the right thumb and right shoulder
    EyeIndices = 6; % Since body landmarks start from 0, you always add up +1 to the index of the desired body part
    WristIndices = 17; %

    % Initialize onset marker variables
    onsetTimes = [];  % Initialize an empty array to store onset times
    onsetFrames = [];  % Initialize an empty array to store onset frames

    % Extract landmarks for the first frame
    landmarks = timeseries_mp(:, 1);

    % Initialize a flag to indicate onset detected
    onsetDetected = false;

    % Define the cooldown duration (e.g., 5000 ms for 5 seconds). This is to
    % avoid detecting movements after the first eye/wrist interaction, in case
    % that participants scratch or something

    cooldownDuration =  10; % 0.5;  for interpolated data
    cooldownEndTime = 0;

    % Loop through the frames starting from the second frame
    % Here I start after 2000 frames because the close eyes recording at
    % the beginning.

    for i = 2000:numFrames %

        % Extract landmarks for the current frame
        landmarks = timeseries_mp(:, i);

        % Extract the Y-coordinates of the right index and right shoulder
        EyeY = landmarks(EyeIndices * 3 - 1);
        WristY = landmarks(WristIndices * 3 - 1);

        % It checks if the Y-coordinate of the eye is greater than or equal to
        % the Y-coordinate of the wrist. This means that the onset will be
        % detected when the eye is at or below the level of the wrist.

        if   EyeY >= WristY

            % Check if an onset hasn't been detected yet and it's not in cooldown
            if ~onsetDetected && timestamps_mp(i) >= cooldownEndTime
                % Onset event detected
                onsetFrames(end+1) = i;  % Store the onset frame
                onsetTimes(end+1) = timestamps_mp(i);  % Store the onset time
                onsetDetected = true;  % Set the flag to indicate onset detected

                % Set the cooldown end time
                cooldownEndTime = timestamps_mp(i) + cooldownDuration;
            end

        else
            % Reset the onset detected flag
            onsetDetected = false;
        end

        % If an onset has been detected and at least 5000 [ms] (adjust as needed) have passed,
        % reset the flag to allow detection of the next onset (next shot)
        if onsetDetected && timestamps_mp(i) >= cooldownEndTime
            onsetDetected = false;
        end
    end


    % Adjusting cooling period time in case of not detection

    if length(onsetTimes) < 120  % If there are not enough trials detected, adjust the cooling period due to participant's speed

        clear onsetDetected
        clear onsetFrames
        clear onsetTimes


        % Determine the number of frames and landmarks
        numFrames = size(timeseries_mp, 2); % Number of frames
        numLandmarks = size(timeseries_mp, 1) / 3; % Number of landmarks (33 according to MediaPipe)

        % Define the indices for the right thumb and right shoulder
        EyeIndices = 6; % Since body landmarks start from 0, you always add up +1 to the index of the desired body part
        WristIndices = 17; %

        % Initialize onset marker variables
        onsetTimes = [];  % Initialize an empty array to store onset times
        onsetFrames = [];  % Initialize an empty array to store onset frames

        % Extract landmarks for the first frame
        landmarks = timeseries_mp(:, 1);

        % Initialize a flag to indicate onset detected
        onsetDetected = false;

        % Define the cooldown duration (e.g., 5000 ms for 5 seconds)
        cooldownDuration =  9; % 0.5;  for interpolated data
        cooldownEndTime = 0;

        % Loop through the frames starting from the second frame
        % Here I start after 2000 frames because the close eyes recording at
        % the beginning.

        for i = 2000:numFrames % in steps of 468 for interpolated data

            % Extract landmarks for the current frame
            landmarks = timeseries_mp(:, i);

            % Extract the Y-coordinates of the right index and right shoulder
            EyeY = landmarks(EyeIndices * 3 - 1);
            WristY = landmarks(WristIndices * 3 - 1);

            % It checks if the Y-coordinate of the eye is greater than or equal to
            % the Y-coordinate of the wrist. This means that the onset will be
            % detected when the eye is at or below the level of the wrist.

            if   EyeY >= WristY

                % Check if an onset hasn't been detected yet and it's not in cooldown
                if ~onsetDetected && timestamps_mp(i) >= cooldownEndTime
                    % Onset event detected
                    onsetFrames(end+1) = i;  % Store the onset frame
                    onsetTimes(end+1) = timestamps_mp(i);  % Store the onset time
                    onsetDetected = true;  % Set the flag to indicate onset detected

                    % Set the cooldown end time
                    cooldownEndTime = timestamps_mp(i) + cooldownDuration;
                end

            else
                % Reset the onset detected flag
                onsetDetected = false;
            end

            % If an onset has been detected and at least 5000 [ms] (adjust as needed) have passed,
            % reset the flag to allow detection of the next onset (next shot)
            if onsetDetected && timestamps_mp(i) >= cooldownEndTime
                onsetDetected = false;
            end
        end
    end


    %% Displaying events

    % Display the onset frames and times
    if ~isempty(onsetFrames)
        %disp(['Total blocks: ', num2str(blockCounter - 1)])
        disp(['Total shots: ', num2str(length(onsetFrames))])
    else
        disp('No onsets detected!');
    end


    % Marker events
    % Average body movement

    % Averaging all 33 channel points to have a common representation
    ave_channels = mean(timeseries_mp, 1);

    % Time series plot for average body data
    figure('units','normalized','outerposition', [0 0 1 1]);
    %subplot 211
    plot(timestamps_mp, ave_channels);
    hold on;
    for i = 1:length(onsetTimes)
        plot([onsetTimes(i), onsetTimes(i)], get(gca, 'YLim'), 'r--');  % Vertical red line at the onset time
        onsetFrameLabel = num2str(onsetFrames(i)); % Label with the onset frame number
        text(onsetTimes(i), max(ave_channels), onsetFrameLabel, 'Color', 'k',...
            'FontSize', 7, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'Rotation', 90, 'FontWeight','bold');
        axis tight; % Adjust Y-axis limits to fit the data (automatic adjustment)
    end
    xlabel('Time [s]', 'FontSize', 11); % Customize the plot
    ylabel('Coordinates [x, y, z]', 'FontSize', 11);
    title('Basketball onset detection [PLD event markers]', 'FontSize', 12);
    subtitle(['Sub. [', num2str(sub), ']',' / Average body coordinates'], 'FontSize', 11.5);
    legend('Average body', 'Onset times', 'Location', 'southwest');  % Show the legend
    grid on;  % Display the plot
    hold off;

    % Wrist/eye movement

    % Extract the Y-coordinate of the right wrist
    WristY = timeseries_mp(WristIndices * 3 - 1, :);
    % Extract the Y-coordinate of the right eye
    EyeY = timeseries_mp(EyeIndices * 3 - 1, :);

    % Time series plot for right wrist and right eye data
    %subplot  212
    plot(timestamps_mp, WristY, 'k', 'LineWidth', 0.8);
    hold on;
    plot(timestamps_mp, EyeY, 'Color', "#0072BD", 'LineWidth', 0.8);

    for i = 1:length(onsetTimes)
        plot([onsetTimes(i), onsetTimes(i)], get(gca, 'YLim'), 'r--');  % Vertical red line at the onset time
        onsetFrameLabel = num2str(onsetFrames(i)); % Label with the onset frame number
        %text(onsetTimes(i), max([WristY, EyeY]), onsetFrameLabel, 'Color', 'k',...
        %    'FontSize', 7, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'Rotation', 90, 'FontWeight','bold');
        axis tight; % Adjust Y-axis limits to fit the data (automatic adjustment)
    end

    xlabel('Time [s]', 'FontSize', 18);
    ylabel('Y-Coordinates', 'FontSize', 18);
    title(['Sub. [', num2str(sub), ']',' / Wrist-eye coordinates'], 'FontSize', 20);
    legend('Wrist', 'Eye', 'Onset times', 'FontSize', 17,'Location', 'northeast');  % Provide custom legends
    grid on;
    hold off;


    % If after all there are still false trials detected, proceed to manual
    % trial rejection inspecting the video
    if length(onsetTimes) > 120  % If there are extra trials detected
        error('There are false onsets detected, inspect the video!')

    end
    %
    % saveas(gcf, [out_subfold, 'PLD_onsets_', participant, '.png']); % Save the figure as a PNG image
    % saveas(gcf, [outpath, '\\group_analysis\\','PLD_onsets_', participant, '.png']); % Save the figure as a PNG image


    %% Deleting false onsets

    if length(onsetTimes) > 120  % number of trials


        % Displaying frame markers for inspection

        % Determine the number of frames and landmarks
        numFrames = size(timeseries_mp, 2); % Number of frames
        numLandmarks = size(timeseries_mp, 1) / 3; % Number of landmarks (33 according to MediaPipe)

        % Labels (33) based on the documentation
        landmarkLabels = {
            'nose', 'left eye (inner)', 'left eye', 'left eye (outer)', 'right eye (inner)',...
            'right eye', 'right eye (outer)', 'left ear', 'right ear', ...
            'mouth (left)', 'mouth (right)', 'left shoulder', 'right shoulder', 'left elbow',...
            'right elbow', 'left wrist', 'right wrist', 'left pinky', 'right pinky', ...
            'left index', 'right index', 'left thumb', 'right thumb', 'left hip', 'right hip',...
            'left knee', 'right knee', 'left ankle', 'right ankle', ...
            'left heel', 'right heel', 'left foot index', 'right foot index'
            };

        % Define colors for each body point
        pointColors = lines(numLandmarks);

        % Define clusters of body parts
        clusters = {[1:10], [11, 13, 15, 17, 19, 21], [12, 14, 16, 18, 20, 22], [23, 24],...
            [12, 11], [24, 26, 28], [23, 25, 27], [28, 30, 32], [27, 29, 31]};

        % Initialize the figure and set axis limits based on the data
        figure;
        minX = min(min(timeseries_mp(1:3:end, :)));
        maxX = max(max(timeseries_mp(1:3:end, :)));

        minY = min(min(timeseries_mp(2:3:end, :)));
        maxY = max(max(timeseries_mp(2:3:end, :)));

        minZ = min(min(timeseries_mp(3:3:end, :)));
        maxZ = max(max(timeseries_mp(3:3:end, :)));

        axisLimits = [minX, maxX, minY, maxY, minZ, maxZ];

        % Loop through the markers to visualize events in body shape
        for i = 1:length(onsetFrames)

            landmarks = timeseries_mp(:, onsetFrames(i));

            % Split the data into X, Y, and Z coordinates for each landmark
            x_coords = landmarks(1:3:end);
            y_coords = landmarks(2:3:end);
            z_coords = landmarks(3:3:end);

            % Create a 3D plot of the landmarks for each body part with unique colors
            figure (onsetFrames(i))
            scatter3(x_coords, y_coords, z_coords, 30, pointColors, 'filled');

            % Add labels for each body point
            for foton = [6, 17]
                text(x_coords(foton), y_coords(foton), z_coords(foton), landmarkLabels{foton},...
                    'Color', pointColors(foton, :));
            end

            % Customize the plot appearance (e.g., title, labels, etc.)
            title('Pose Landmark Detection');
            subtitle(['Frame ', num2str(onsetFrames(i))])
            subtitle(['Frame: ', num2str(onsetFrames(i)), ' at ', num2str(timestamps_mp(onsetFrames(i))), ' [s]'])
            xlabel('X-coordinate');
            ylabel('Y-coordinate');
            zlabel('Z-coordinate');

            % Set axis limits for consistent scaling
            axis(axisLimits);

            % View with the desired orientation
            % view(3) %for 3D
            view(0, -90);

            % Connect body parts within clusters with lines
            for b_part = 1:length(clusters)
                cluster_points = clusters{b_part};
                for foton = 1:(length(cluster_points) - 1)
                    idx1 = cluster_points(foton);
                    idx2 = cluster_points(foton + 1);
                    line([x_coords(idx1), x_coords(idx2)],...
                        [y_coords(idx1), y_coords(idx2)],...
                        [z_coords(idx1), z_coords(idx2)], 'Color', pointColors(idx1, :));
                end
            end

            % Display the plot
            grid on;
            axis equal;

        end

    end


    %%

    % % After inspection, delete bad trial frames here:
    % framesToDelete = [];
    %
    % % Find indices of frames to delete
    % indicesToDelete = find(ismember(onsetFrames, framesToDelete));
    %
    % % Delete corresponding frames and onsetTimes
    % onsetFrames(indicesToDelete) = []; % MediaPipe crush
    % onsetTimes(indicesToDelete) = [];

    % This part of the code is only used in case that the participants raised
    % their hands for some reason during the experiment, which will be detected
    % as an onset. Therefore, for the ones that have more than 120 shots
    % detected, one can plot the frames of those timepoints and compare whith
    % the video recordings to know what happened. Usually one can see that most
    % of those false detected osets are in the trial break, when the
    % participants are sitting down, by drinking water for example, which is
    % visible in the PLD frame.


    %% Assigning labels to the markers by table

    % Load the labels from an Excel file
    labelsTable = readtable([out_subfold, participant, '.xlsx']);

    % Extract block names and shot labels
    blockNames = labelsTable.Properties.VariableNames;
    shotLabels = table2array(labelsTable);

    % Create events structure
    numEvents = numel(shotLabels);
    events_MP = struct('type', cell(1, numEvents), 'latency', zeros(1, numEvents), 'urevent', zeros(1, numEvents), 'duration', ones(1, numEvents));

    % Assign 'hit', 'miss', or 'none' for each event based on the labels table
    for shot = 1:numEvents
        events_MP(shot).latency = onsetFrames(shot);
        events_MP(shot).time = onsetTimes(shot);

        % Extract block name for the current shot
        currentBlock = blockNames{ceil(shot / size(shotLabels, 1))};

        % Check label and assign type accordingly
        if shotLabels(shot) == 1
            events_MP(shot).type = 'hit_MP';
        elseif shotLabels(shot) == 2
            events_MP(shot).type = 'miss_MP';
        else
            error('Invalid label. Label must be 1 (hit) or 2 (miss).');
        end

        events_MP(shot).urevent = shot;
        events_MP(shot).duration = 1; % samples
        events_MP(shot).block = currentBlock; % Add block information
    end


    %% Extra: Interpolation procedure illustration

    % Time series plot for average body data
    figure('units','normalized','outerposition', [0 0 1 1]);
    subplot 211
    plot(timestamps_eeg, mean(timeseries_eeg, 1), 'Color', 	"#A2142F");
    xlabel('Time [s]', 'FontSize', 11);
    ylabel('Amplitude [\muV]', 'FontSize', 11, 'Color', 'k');
    title('EEG raw data', 'FontSize', 12);
    subtitle(['Sub. [', num2str(sub), ']',' / Mean time series across electrodes'], 'FontSize', 11.5);
    legend('EEG raw signal', 'Onset times', 'Location', 'southwest');  % Show the legend
    axis tight;
    grid on;  % Display the plot

    subplot 212
    plot(timestamps_mp, ave_channels);
    xlabel('Time [s]', 'FontSize', 11); % Customize the plot
    ylabel('Coordinates [x, y, z]', 'FontSize', 11);
    title('Interpolated PLD raw data', 'FontSize', 12);
    subtitle(['Sub. [', num2str(sub), ']',' / Mean time series across coordinates'], 'FontSize', 11.5);
    legend('Body raw signal', 'Onset times', 'Location', 'southwest');  % Show the legend
    axis tight;
    grid on;  % Display the plot


    % saveas(gcf, [out_subfold, 'PLD_interpolation_', participant, '.png']); % Save the figure as a PNG image
    % saveas(gcf, [outpath, '\\group_analysis\\','PLD_interpolation_', participant, '.png']); % Save the figure as a PNG image


    %% Saving events

    mkdir(outpath, participant); % creating participant sub-folders

    % % Save it in .mat file
    % save([out_subfold, 'events_MP_', participant,'.mat'],'events_MP', 'onsetTimes', 'onsetFrames', ...
    %     'timeseries_mp', 'timestamps_mp', 'sampling_rate_eeg', 'sampling_rate_mp');


    %%

    disp([participant, ' finalized!']);

end