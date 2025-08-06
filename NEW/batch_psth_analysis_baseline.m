%% BATCH ANALYSIS OF KERNEL DENSITY ESTIMATION FOR PSTH (Multi-Baseline Corrected with Rmax/Half-Max)
% This script loops through all channels in a recording. For each channel,
% it generates a single 2x2 figure that includes:
% 1. A spike raster plot for the response window.
% 2. A diagnostic plot comparing the response PSTH against 50 pre-stimulus baseline PSTHs.
% 3. The smoothed response PSTH with key metrics (Rmax, Half-Rmax).
% 4. The adaptive kernel bandwidth over time for the response PSTH.
% All metrics, including a baseline-corrected AUC, are saved.

clear; clc;

%% --- USER PARAMETERS ---
spikeFile = 'OWT220207_1I_DIV63_PER37_6UA_Cspikes_L0_RP2.mat_Nortefact.mat';
rawFile   = 'OWT220207_1I_DIV63_PER37_6UA.mat';
fs = 25000;
spikeMethod = 'bior1p5';
numChannels = 60;
artifact_window_ms = [0, 2];
psth_window_s = [0, 0.02]; % Analysis window (0 to 20ms post-stimulus)
psth_bin_width_s = 0.001;
num_baseline_psths = 100; % Number of baseline PSTHs to calculate backwards
baseline_duration_s = 0.02; % Duration of each baseline window (20ms)

%% --- CHANNEL REMAPPING & SETUP ---
indices = [24 26 29 32 35 37, 21 22 25 30 31 36 39 40, 19 20 23 28 33 38 41 42, 16 17 18 27 34 43 44 45, 15 14 13 4 57 48 47 46, 12 11 8 3 58 53 50 49, 10 9 6 1 60 55 52 51, 7 5 2 59 56 54];
ids = [21 31 41 51 61 71, 12 22 32 42 52 62 72 82, 13 23 33 43 53 63 73 83, 14 24 34 44 54 64 74 84, 15 25 35 45 55 65 75 85, 16 26 36 46 56 66 76 86, 17 27 37 47 57 67 77 87, 28 38 48 58 68 78];
channelMap = containers.Map('KeyType','double','ValueType','double');
for i = 1:numel(indices), channelMap(indices(i)) = ids(i); end
timestamp = datestr(now, 'ddmmmyyyy_HHMMSS');
outputDir = ['PSTHanalysis_multi_baseline_corr(' timestamp ')'];
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
fprintf('Saving analysis plots to folder: %s\n', outputDir);

%% --- LOAD DATA & FIND STIMS ---
fprintf('Loading and processing data...\n');
S = load(spikeFile);
if isfield(S, 'spikeTimes'), spikeTimesConverted = S.spikeTimes;
elseif isfield(S, 'spikes'), fprintf('Converting ''spikes'' matrix...\n');
    [row, col] = find(S.spikes); spikeTimesConverted = cell(1, numChannels);
    for ch = 1:numChannels, spike_samples = row(col == ch); spike_sec = spike_samples / fs; spikeTimesConverted{ch} = struct(spikeMethod, spike_sec); end
else, error('Spike data not found.'); end
R = load(rawFile);
if isfield(R, 'dat'), dat = double(R.dat); else, error('Raw data not found.'); end
[num_samples, ~] = size(dat); stimThreshold = -1000; min_interval_ms = 2500; flat_search_window_ms = 100; flat_window_ms = 1.5; flat_thresh = 0.05;
flat_window_samples = round(flat_window_ms * fs / 1000); flat_search_samples = round(flat_search_window_ms * fs / 1000); stim_times_sec = [];
for channel_idx = 1:numChannels, trace = dat(:, channel_idx); idx = find(trace > stimThreshold); if isempty(idx), continue, end; idx = idx(:);
    keep = [true; diff(idx) > round(0.010 * fs)]; idx = idx(keep);
    for i = 1:length(idx), center_idx = idx(i); win_start = max(1, center_idx - flat_search_samples); win_end = min(num_samples, center_idx + flat_search_samples);
        win_trace = trace(win_start:win_end); abs_diff = [0; abs(diff(win_trace))]; mov_abs_diff = movmean(abs_diff, flat_window_samples);
        flat_idx = find(mov_abs_diff < flat_thresh);
        if ~isempty(flat_idx), flat_onsets = flat_idx([true; diff(flat_idx) > 1]); flat_onsets_adj = flat_onsets - floor(flat_window_samples/2); flat_onsets_adj(flat_onsets_adj < 1) = 1;
            if ~isempty(flat_onsets_adj), keep_idx = [true; diff(flat_onsets_adj) > round(min_interval_ms * fs / 1000)]; flat_onsets_adj = flat_onsets_adj(keep_idx); end
            for j = 1:length(flat_onsets_adj), stim_idx = win_start - 1 + flat_onsets_adj(j); stim_times_sec = [stim_times_sec; stim_idx / fs]; end
        end; end; end
stimTimes = sort(unique(stim_times_sec(:)));
fprintf('Found %d stimulation events.\n', length(stimTimes));

%% --- LOOP THROUGH CHANNELS FOR ANALYSIS ---
networkResponse = [];
valid_channel_count = 0;
for file_idx = 1:numChannels
    if ~isKey(channelMap, file_idx), continue; end
    channel_id = channelMap(file_idx);
    fprintf('\n--- Processing File Index %d (Channel ID %d) ---\n', file_idx, channel_id);

    if file_idx > numChannels || file_idx < 1 || isempty(spikeTimesConverted{file_idx}) || ~isfield(spikeTimesConverted{file_idx}, spikeMethod), all_spike_times_s = [];
    else, all_spike_times_s = spikeTimesConverted{file_idx}.(spikeMethod); end
    if isempty(all_spike_times_s), fprintf('No spikes found. Skipping.\n'); continue; end
    
    spikeTimes_cleaned_s = all_spike_times_s;
    for stimIdx = 1:numel(stimTimes), stimTime = stimTimes(stimIdx);
        spikeTimes_cleaned_s = spikeTimes_cleaned_s(spikeTimes_cleaned_s < (stimTime + artifact_window_ms(1)/1000) | spikeTimes_cleaned_s >= (stimTime + artifact_window_ms(2)/1000));
    end
    spikeTimes_cleaned_s = sort(spikeTimes_cleaned_s(:));

    % --- PSTH Calculation for Response Window ---
    [response, resp_metrics] = calculate_psth_metrics(spikeTimes_cleaned_s, stimTimes, psth_window_s, psth_bin_width_s);
    if isempty(response.psth_samples), fprintf('No spikes in PSTH window. Skipping.\n'); continue; end
    
    % --- PSTH Calculation for Multiple Baseline Windows ---
    fprintf('Calculating %d baseline PSTHs...\n', num_baseline_psths);
    baseline_aucs = zeros(num_baseline_psths, 1);
    all_baseline_psth_smooth = [];
    for i = 1:num_baseline_psths
        start_s = -(i * baseline_duration_s);
        end_s = -((i-1) * baseline_duration_s);
        current_baseline_window_s = [start_s, end_s];
        [~, base_metrics] = calculate_psth_metrics(spikeTimes_cleaned_s, stimTimes, current_baseline_window_s, psth_bin_width_s);
        baseline_aucs(i) = base_metrics.auc;
        if isempty(all_baseline_psth_smooth), all_baseline_psth_smooth = zeros(num_baseline_psths, length(base_metrics.psth_smooth)); end
        all_baseline_psth_smooth(i, :) = base_metrics.psth_smooth;
    end
    
    % --- Calculate Corrected AUC ---
    mean_baseline_auc = mean(baseline_aucs);
    auc_corrected = resp_metrics.auc - mean_baseline_auc;
    mean_baseline_psth = mean(all_baseline_psth_smooth, 1);
    fprintf('Response AUC: %.4f, Mean Baseline AUC: %.4f, Corrected AUC: %.4f\n', resp_metrics.auc, mean_baseline_auc, auc_corrected);

    % --- Find Half-Max Following Peak ---
    [Rmax, Rmax_idx] = max(resp_metrics.psth_smooth);
    halfRmax = Rmax / 2;
    halfRmax_idx = find(resp_metrics.psth_smooth(Rmax_idx:end) <= halfRmax, 1, 'first');
    if ~isempty(halfRmax_idx)
        halfRmax_idx = halfRmax_idx + Rmax_idx - 1;
        halfRmax_time_s = resp_metrics.time_vector_s(halfRmax_idx);
        halfRmax_val = resp_metrics.psth_smooth(halfRmax_idx);
    else
        halfRmax_time_s = NaN;
        halfRmax_val = NaN;
    end

    %% --- CONSOLIDATED 2x2 PLOT GENERATION ---
    fig = figure('Position',[100 100 1200 900], 'Visible', 'off');
    psth_window_ms = psth_window_s * 1000;
    sgtitle(sprintf('Channel %d | Peak Rate: %.1f Hz | Corrected AUC: %.3f', channel_id, resp_metrics.peak_firing_rate, auc_corrected), 'FontWeight', 'bold');

    % Plot 1: Raster Plot
    ax1 = subplot(2,2,1); hold on;
    for trial_idx = 1:length(response.spikeTimes_byEvent), trial_spikes_s = response.spikeTimes_byEvent{trial_idx}; if ~isempty(trial_spikes_s), plot(trial_spikes_s * 1000, trial_idx * ones(size(trial_spikes_s)), 'r.', 'MarkerSize', 5); end; end
    hold off; set(gca, 'YDir', 'reverse'); xlim(psth_window_ms); ylim([0 length(stimTimes)+1]); ylabel('Trial Number'); title('Spike Raster (Response)'); grid on;

    % Plot 2: Diagnostic Plot (Response vs. Baselines)
    ax2 = subplot(2,2,2); hold on;
    baseline_time_ms = (base_metrics.time_vector_s - current_baseline_window_s(1)) * 1000;
    for i = 1:num_baseline_psths, plot(baseline_time_ms, all_baseline_psth_smooth(i, :), 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5); end
    p1_diag = plot(resp_metrics.time_vector_s*1000, resp_metrics.psth_smooth, 'r-', 'LineWidth', 2);
    p2_diag = plot(baseline_time_ms, mean_baseline_psth, 'k-', 'LineWidth', 2);
    hold off; title('Diagnostic: Response vs. Baselines'); ylabel('Firing Rate (spikes/s)'); xlabel('Time from stimulus (ms)'); legend([p1_diag, p2_diag], 'Response', 'Mean Baseline', 'Location', 'Best'); grid on;

    % Plot 3: Smoothed PSTH and Metrics
    ax3 = subplot(2,2,3); hold on;
    edges_s = psth_window_s(1):psth_bin_width_s:psth_window_s(2);
    bar(edges_s * 1000, response.psth_histogram, 1, 'FaceColor',.7*[1 1 1],'EdgeColor',.8*[1 1 1], 'HandleVisibility','off');
    plot(resp_metrics.time_vector_s * 1000, resp_metrics.psth_smooth, 'r-', 'LineWidth', 2, 'DisplayName', 'Smoothed PSTH');
    plot(resp_metrics.peak_time_s*1000, resp_metrics.peak_firing_rate, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 7, 'DisplayName', 'R_{max}');
    text(resp_metrics.peak_time_s*1000, resp_metrics.peak_firing_rate, sprintf(' R_{max}: %.1f Hz', resp_metrics.peak_firing_rate), 'VerticalAlignment', 'bottom', 'HorizontalAlignment','left', 'Color', 'k', 'FontWeight','bold');
    if ~isnan(halfRmax_time_s)
        plot(halfRmax_time_s*1000, halfRmax_val, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 7, 'DisplayName', 'Half R_{max}');
        text(halfRmax_time_s*1000, halfRmax_val, sprintf(' Half R_{max} @ %.1f ms', halfRmax_time_s*1000), 'VerticalAlignment', 'top', 'HorizontalAlignment','left', 'Color', 'k', 'FontWeight','bold');
    end
    hold off; ylabel('Firing Rate (spikes/s)'); xlabel('Time from stimulus (ms)'); title('Smoothed Response PSTH & Metrics'); legend('Location', 'northeast'); grid on;

    % Plot 4: Kernel Bandwidth
    ax4 = subplot(2,2,4);
    plot(resp_metrics.time_vector_s * 1000, resp_metrics.kernel_bandwidth_s * 1000, 'b-', 'LineWidth', 2);
    ylabel('Bandwidth (ms)'); xlabel('Time from stimulus (ms)'); title('Adaptive Kernel Bandwidth (Response)'); grid on;
    
    linkaxes([ax1, ax2, ax3, ax4], 'x'); xlim(psth_window_ms);
    plot_filename = fullfile(outputDir, sprintf('Analysis_Plot_Chan_%d.png', channel_id));
    saveas(fig, plot_filename); fprintf('Saved consolidated plot to %s\n', plot_filename); close(fig);

  %% --- SAVE RESULTS ---
    valid_channel_count = valid_channel_count + 1;
    networkResponse(valid_channel_count).channel_id = channel_id;
    networkResponse(valid_channel_count).file_index = file_idx;
    networkResponse(valid_channel_count).auc_response = resp_metrics.auc;
    networkResponse(valid_channel_count).auc_baseline_mean = mean_baseline_auc;
    networkResponse(valid_channel_count).auc_corrected = auc_corrected;
    networkResponse(valid_channel_count).peak_firing_rate_hz = resp_metrics.peak_firing_rate;
    networkResponse(valid_channel_count).peak_time_ms = resp_metrics.peak_time_s * 1000;
    networkResponse(valid_channel_count).halfRmax_time_ms = halfRmax_time_s * 1000;
end

%% --- SAVE AGGREGATE RESULTS TO OUTPUT DIRECTORY ---
if ~isempty(networkResponse)
    output_filename = fullfile(outputDir, sprintf('networkResponse_all_channels_%s.mat', timestamp));
    save(output_filename, 'networkResponse');
    fprintf('\nSaved aggregate network response metrics for %d channels to %s\n', valid_channel_count, output_filename);
else, fprintf('\nNo valid channels with spikes were found to save.\n'); end
fprintf('Script finished.\n');
