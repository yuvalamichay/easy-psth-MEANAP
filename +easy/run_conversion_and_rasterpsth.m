% run_conversion_and_rasterpsth.m

% --- USER PARAMETERS: CONVERSION ---

% Path to the .mat file containing spike times or spike matrix.
spikeFile = 'OWT220207_1I_DIV63_HUB63_6UA_spikes.mat';

% Path to the .mat file containing the raw electrophysiological data.
rawFile   = 'OWT220207_1I_DIV63_HUB63_6UA.mat';

% Time window (in milliseconds) after each stimulus during which detected spikes are considered artifacts and will be removed.
artifact_window_ms = [0, 2];

% Sampling frequency (Hz) of the recording.
fs = 25000;

% Name of the spike detection method (should match the field name in the loaded spike data).
spikeMethod = 'bior1p5';

% Total number of channels in the recording. Used to allocate arrays and loop over channels.
numChannels = 60;

% Index of the channel to extract spike times from and analyze (1-based index).
channelToExtract = 5;

% --- USER PARAMETERS: RasterPSTH ---

% Time window (in seconds) around each event (stimulus) to display in the raster/PSTH.
rasterpsth_params.window = [-1 1];

% Bin width (in seconds) for the PSTH histogram. 
rasterpsth_params.psthBinWidth = 0.001;

% Standard deviation (in seconds) of the Gaussian smoothing kernel applied to the PSTH.
% Use 0 for no smoothing.
rasterpsth_params.psthSmoothWidth = 0.05;

% (Optional) Provide a Kilosort template to use for spike sorting visualization.
% Leave empty if not using Kilosort output.
rasterpsth_params.kilosortTemplate = [];

% Title text for the raster/PSTH plot.
rasterpsth_params.titleText = 'Raster/PSTH';

% (Optional) Vector or cell array to split trials into groups (e.g., by stimulus type).
% Leave empty if not splitting.
rasterpsth_params.splitBy = [];

% (Optional) Field to sort trials by (e.g., reaction time).
% Leave empty for no sorting.
rasterpsth_params.sortBy = [];

% (Optional) Custom colors for different split groups.
% Leave empty for default colors.
rasterpsth_params.splitByColours = [];

% (Optional) Set y-axis limits for the plot.
% Leave empty for automatic scaling.
rasterpsth_params.ylim = [];

% (Optional) Path to save the generated figure (e.g., 'output/myplot.pdf').
% Leave empty to not save automatically.
rasterpsth_params.saveFigure = '';
%% --- CONVERSION SECTION (from extracting/conversion_script.m) ---
S = load(spikeFile);
if isfield(S, 'spikeTimes')
    spiketimesconverted = S.spikeTimes;
elseif isfield(S, 'spikes')
    disp('Converting ''spikes'' matrix to ''spiketimesconverted'' struct format...');
    [row, col] = find(S.spikes);
    spiketimesconverted = cell(1, numChannels);
    for ch = 1:numChannels
        spike_samples = row(col == ch);
        spike_sec = spike_samples / fs;
        spiketimesconverted{ch} = struct(spikeMethod, spike_sec);
    end
else
    error('Neither ''spikeTimes'' nor ''spikes'' found in file.');
end

R = load(rawFile);
if isfield(R, 'dat')
    dat = double(R.dat);
else
    error('Raw data variable "dat" not found in %s', rawFile);
end
[num_samples, ~] = size(dat);

stimThreshold = -1000;
flat_window_ms = 1.5;
flat_thresh = 0.05;
min_interval_ms = 2500;
flat_search_window_ms = 100;
flat_window_samples = round(flat_window_ms * fs / 1000);
flat_search_samples = round(flat_search_window_ms * fs / 1000);

stim_times_sec = [];
for channel_idx = 1:numChannels
    trace = dat(:, channel_idx);
    idx = find(trace > stimThreshold);
    if isempty(idx), continue, end
    idx = idx(:);
    keep = [true; diff(idx) > round(0.010 * fs)];
    idx = idx(keep);
    for i = 1:length(idx)
        center_idx = idx(i);
        win_start = max(1, center_idx - flat_search_samples);
        win_end = min(num_samples, center_idx + flat_search_samples);
        win_trace = trace(win_start:win_end);
        abs_diff = [0; abs(diff(win_trace))];
        mov_abs_diff = movmean(abs_diff, flat_window_samples);
        flat_idx = find(mov_abs_diff < flat_thresh);
        if ~isempty(flat_idx)
            flat_onsets = flat_idx([true; diff(flat_idx) > 1]);
            flat_onsets_adj = flat_onsets - floor(flat_window_samples/2);
            flat_onsets_adj(flat_onsets_adj < 1) = 1;
            if ~isempty(flat_onsets_adj)
                keep_idx = [true; diff(flat_onsets_adj) > round(min_interval_ms * fs / 1000)];
                flat_onsets_adj = flat_onsets_adj(keep_idx);
            end
            for j = 1:length(flat_onsets_adj)
                stim_idx = win_start - 1 + flat_onsets_adj(j);
                stim_time_sec = stim_idx / fs;
                stim_times_sec = [stim_times_sec; stim_time_sec]; %#ok<AGROW>
            end
        end
    end
end
stim_times_sec = sort(stim_times_sec);
eventTimes = stim_times_sec(:); % [T x 1] vector, seconds

% --- EXTRACT SPIKE TIMES FOR ONE CHANNEL AND REMOVE ARTIFACT SPIKES ---
if isempty(spiketimesconverted{channelToExtract}) || ~isfield(spiketimesconverted{channelToExtract}, spikeMethod)
    spikeTimes = [];
else
    spikeTimes_sec = spiketimesconverted{channelToExtract}.(spikeMethod); % s
    % Remove spikes near any stim
    for stimIdx = 1:numel(eventTimes)
        stimTime = eventTimes(stimIdx); % s
        spikeTimes_sec = spikeTimes_sec(...
            spikeTimes_sec < (stimTime + artifact_window_ms(1)/1000) | ...
            spikeTimes_sec >= (stimTime + artifact_window_ms(2)/1000) ...
            );
    end
    spikeTimes = sort(spikeTimes_sec(:)); % [S x 1] vector, seconds
end

% Optionally save for re-use:
save(spikeFile, 'spikeTimes', 'eventTimes', '-append');
fprintf('Saved spikeTimes & eventTimes for channel %d.\n', channelToExtract);

%% --- CALL RasterPSTH ---
% Set up events cell array: label + eventTimes
events = {'stim', eventTimes};

% Build Name,Value argument pairs dynamically from rasterpsth_params struct
param_names = fieldnames(rasterpsth_params);
param_cell = {};
for i = 1:numel(param_names)
    param_cell = [param_cell, param_names{i}, rasterpsth_params.(param_names{i})];
end

% Call RasterPSTH
grammObject = easy.RasterPSTH(spikeTimes, events, param_cell{:});
