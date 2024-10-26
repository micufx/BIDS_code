clc, clear, close all;

%% EEG condition analysis

mainpath = 'C:\Users\micua\Desktop\eeglab2023.0\'; % eeglab folder
path = 'C:\Users\micua\OneDrive - Benemérita Universidad Autónoma de Puebla\NCP_Basketball\MediaPipe\';
outpath = 'C:\\Users\\micua\\OneDrive - Benemérita Universidad Autónoma de Puebla\\Oldenburg_University\\Thesis\\data_hoops\\';
files = dir( fullfile( path,'\*.xdf')); % listing data sets

nochans = {'AccX','AccY','AccZ','GyroX','GyroY','GyroZ', ... % channels to be ignored
    'QuatW','QuatX','QuatY','QuatZ'};

conditions = {'hit', 'miss'};
num_conditions = 2; % (Conditions: 1=hit 2=miss)

load([outpath, 'means_Cz_hits.mat']); % Loading channels file
load([outpath, 'means_Cz_miss.mat']); % Loading channels file


%% Distribution box plot (Mean comparison across participants)

% Example data - replace these with your actual means and SEs
mean_hits = mean(means_Cz_hits{1:end-1, :});  % Mean of hits for each bin
se_hits = std(means_Cz_hits{1:end-1, :}) / sqrt(size(means_Cz_hits, 1) - 1); % SE for hits

mean_misses = mean(means_Cz_miss{1:end-1, :}); % Mean of misses for each bin
se_misses = std(means_Cz_miss{1:end-1, :}) / sqrt(size(means_Cz_miss, 1) - 1); % SE for misses

% The bins you have - this should match the number of bins you actually have
bin_labels = {'[-1500:-1400 ms]', '[-1400:-1300 ms]', '[-1300:-1200 ms]',...
    '[-1200:-1100 ms]', '[-1100:-1000 ms]', '[-1000:-900 ms]',...
    '[-900:-800 ms]', '[-800:-700 ms]', '[-700:-600 ms]', ...
    '[-600:-500 ms]', '[-500:-400 ms]',...
    '[-400:-300ms]', '[-300:-200 ms]', '[-200:-100 ms]', '[-100:0 ms]'};


% Create figure
figure('units','normalized','outerposition', [0 0 1 1]);
hold on;

% Define the number of bins
nBins = numel(bin_labels);

% Calculate the width for each bar group
groupWidth = min(0.8, nBins/(nBins + 1.5));

% Plot bars and error bars for each bin
for i = 1:nBins
    % Set the positions for hits and misses bars in each bin
    positions_hits = i - groupWidth / 4;
    positions_misses = i + groupWidth / 4;
    
    % Plot bars for hits
    bar(positions_hits, mean_hits(i), groupWidth / 4, 'FaceColor', [0 .7 .7]);
    
    % Plot bars for misses
    bar(positions_misses, mean_misses(i), groupWidth / 4, 'FaceColor', 'r');
    
    % Error bars for hits
    errorbar(positions_hits, mean_hits(i), se_hits(i), 'k', 'linestyle', 'none');
    
    % Error bars for misses
    errorbar(positions_misses, mean_misses(i), se_misses(i), 'k', 'linestyle', 'none');
end

% Customizing the plot
set(gca, 'xtick', 1:nBins, 'xticklabel', bin_labels, 'XTickLabelRotation', 30, 'FontSize', 11);
ylabel('Amplitude [µV]', 'FontSize', 12);
title('Mean Amplitude Across Participants', 'FontSize', 12);
subtitle('Readiness Potential [Cz]', 'FontSize', 12);
legend({'Hits', 'Misses'}, 'Location', 'Best', 'FontSize', 10);
grid on;

hold off;


% Saving figures
saveas(gcf, [outpath, '\\group_analysis\\','box_mean_comparison', '.png']); % Save the figure as a PNG image
% saveas(gcf, [outpath,'box_mean_comparison', '.png']); % Save the figure as a PNG image


%% Save files

% Save it in .mat file
save([outpath, 'means_Cz_hits', '.xlsx'], 'means_Cz_hits');

% Save it in .mat file
save([outpath, 'means_Cz_miss', '.xlsx'], 'means_Cz_miss');

%%