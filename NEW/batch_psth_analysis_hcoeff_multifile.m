%% BATCH ANALYSIS FOR MULTIPLE FILES (H-COEFFICIENT & KERNEL DENSITY ESTIMATION)
% This script loops through a user-defined list of spike and raw file pairs.
% For each pair, it creates a separate output folder and runs a complete
% PSTH analysis using the h-coefficient method.
% It generates a 2x2 summary plot for each channel and saves all metrics
% into a single aggregate results file per file pair.

clear; clc;

%% --- USER PARAMETERS ---

% --- List of file pairs to analyze ---
% Add as many structs to this cell array as needed.
% Each struct must contain a 'spikeFile' and a 'rawFile'.
filePairs = { ...
    struct('spikeFile', 'OWT220207_1I_DIV63_HUB63_6UA_Cspikes_L0_RP2.mat_Nortefact.mat', 'rawFile', 'OWT220207_1I_DIV63_HUB63_6UA.mat') ...
    ,struct('spikeFile', 'OWT220207_2B_DIV63_HUB24_6UA_Cspikes_L0_RP2.mat_Nortefact.mat', 'rawFile', 'OWT220207_2B_DIV63_HUB24_6UA.mat') ...
    ,struct('spikeFile', 'OWT220207_2D_DIV63_HUB73_6UA_Cspikes_L0_RP2.mat_Nortefact.mat', 'rawFile', 'OWT220207_2D_DIV63_HUB73_6UA.mat') ...
    ,struct('spikeFile', 'OWT220207_2E_DIV63_HUB36_6UA_Cspikes_L0_RP2.mat_Nortefact.mat', 'rawFile', 'OWT220207_2E_DIV63_HUB36_6UA.mat') ...
    ,struct('spikeFile', 'OWT220207_2I_DIV63_HUB24_6UA_Cspikes_L0_RP2.mat_Nortefact.mat', 'rawFile', 'OWT220207_2I_DIV63_HUB24_6UA.mat') ...
    ,struct('spikeFile', 'OWT220207_1I_DIV63_PER37_6UA_Cspikes_L0_RP2.mat_Nortefact.mat', 'rawFile', 'OWT220207_1I_DIV63_PER37_6UA.mat') ...
    ,struct('spikeFile', 'OWT220207_2B_DIV63_PER52_6UA_Cspikes_L0_RP2.mat_Nortefact.mat', 'rawFile', 'OWT220207_2B_DIV63_PER52_6UA.mat') ...
    ,struct('spikeFile', 'OWT220207_2D_DIV63_PER58_6UA_Cspikes_L0_RP2.mat_Nortefact.mat', 'rawFile', 'OWT220207_2D_DIV63_PER58_6UA.mat') ...
    ,struct('spikeFile', 'OWT220207_2E_DIV63_PER71_6UA_Cspikes_L0_RP2.mat_Nortefact.mat', 'rawFile', 'OWT220207_2E_DIV63_PER71_6UA.mat') ...
    ,struct('spikeFile', 'OWT220207_2I_DIV63_PER83_6UA_Cspikes_L0_RP2.mat_Nortefact.mat', 'rawFile', 'OWT220207_2I_DIV63_PER83_6UA.mat') ...
    % --- Add more file pairs below ---
    % Example:
    % ,struct('spikeFile', 'path/to/your/second_spike_file.mat', 'rawFile', 'path/to/your/second_raw_file.mat') ...
    % ,struct('spikeFile', 'path/to/your/third_spike_file.mat', 'rawFile', 'path/to/your/third_raw_file.mat') ...
};

% --- Analysis Parameters ---
fs = 25000;
spikeMethod = 'bior1p5';
numChannels = 60;
artifact_window_ms = [0, 2];
psth_window_s = [0, 0.02];      % Analysis window (e.g., 0 to 20ms post-stimulus)
psth_bin_width_s = 0.001;

% --- Channel Remapping (assumes this is consistent across files) ---
indices = [24 26 29 32 35 37, 21 22 25 30 31 36 39 40, 19 20 23 28 33 38 41 42, 16 17 18 27 34 43 44 45, 15 14 13 4 57 48 47 46, 12 11 8 3 58 53 50 49, 10 9 6 1 60 55 52 51, 7 5 2 59 56 54];
ids = [21 31 41 51 61 71, 12 22 32 42 52 62 72 82, 13 23 33 43 53 63 73 83, 14 24 34 44 54 64 74 84, 15 25 35 45 55 65 75 85, 16 26 36 46 56 66 76 86, 17 27 37 47 57 67 77 87, 28 38 48 58 68 78];
channelMap = containers.Map('KeyType','double','ValueType','double');
for i = 1:numel(indices), channelMap(indices(i)) = ids(i); end


%% --- MAIN LOOP TO PROCESS EACH FILE PAIR ---
for k = 1:numel(filePairs)
    % --- Get current file pair ---
    spikeFile = filePairs{k}.spikeFile;
    rawFile = filePairs{k}.rawFile;

    fprintf('\n\n============================================================\n');
    fprintf('STARTING ANALYSIS FOR FILE PAIR %d of %d:\n', k, numel(filePairs));
    fprintf('  Spike File: %s\n', spikeFile);
    fprintf('  Raw File: %s\n', rawFile);
    fprintf('============================================================\n');

    % --- Generate a unique output directory for this file pair ---
    [~, baseName, ~] = fileparts(spikeFile);
    timestamp = datestr(now, 'ddmmmyyyy_HH:MM');
    outputDir = sprintf('PSTH_HCoeff_Analysis_%s_%s', baseName, timestamp);
    if ~exist(outputDir, 'dir'), mkdir(outputDir); end
    fprintf('Saving analysis plots to folder: %s\n', outputDir);

    % --- LOAD DATA & FIND STIMS ---
    fprintf('Loading and processing data...\n');
    S = load(spikeFile);
    if isfield(S, 'spikeTimes'), spikeTimesConverted = S.spikeTimes;
    elseif isfield(S, 'spikes'), fprintf('Converting ''spikes'' matrix...\n');
        [row, col] = find(S.spikes); spikeTimesConverted = cell(1, numChannels);
        for ch = 1:numChannels, spike_samples = row(col == ch); spike_sec = spike_samples / fs; spikeTimesConverted{ch} = struct(spikeMethod, spike_sec); end
    else, error('Spike data not found in file: %s', spikeFile); end

    R = load(rawFile);
    if isfield(R, 'dat'), dat = double(R.dat); else, error('Raw data not found in file: %s', rawFile); end

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

    % --- LOOP THROUGH CHANNELS FOR ANALYSIS ---
    networkResponse = []; % Reset results for each file pair
    valid_channel_count = 0; % Reset counter for each file pair
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

        out = WithinRanges(spikeTimes_cleaned_s, stimTimes + psth_window_s, (1:length(stimTimes))', 'matrix');
        spikeTimes_byEvent = arrayfun(@(n) spikeTimes_cleaned_s(logical(out(:,n))) - stimTimes(n), 1:length(stimTimes), 'uni', 0)';
        psth_samples = cell2mat(spikeTimes_byEvent);

        if isempty(psth_samples), fprintf('No spikes in PSTH window. Skipping.\n'); continue; end

        % --- H-COEFFICIENT & SMOOTHED PSTH CALCULATION ---
        fprintf('Calculating h-coefficient and smoothed PSTH...\n');
        [hcoeffs, hcoeffs2D, ~, bslshuff, mfr] = hcoeff_wrapper(spikeTimes_byEvent, all_spike_times_s, psth_window_s, psth_bin_width_s);
        fprintf('h-coeff (max/min): %.2f / %.2f | Mean Firing Rate: %.2f Hz\n', hcoeffs(1), hcoeffs(2), mfr);

        L = 1000; t_s = linspace(psth_window_s(1), psth_window_s(2), L);
        [yv_pdf, tv_s, optw_variable_s] = ssvkernel(psth_samples, t_s);
        num_trials = length(stimTimes);
        avg_spikes_per_trial = length(psth_samples) / num_trials;
        response_psth_smooth = yv_pdf * avg_spikes_per_trial;

        % --- CALCULATE AUC (RESPONSE & CORRECTED) ---
        auc_response = trapz(tv_s, response_psth_smooth);
        
        baseline_aucs = [];
        if ~isempty(bslshuff) && isfield(bslshuff, 'kern')
            baseline_aucs = zeros(length(bslshuff.kern), 1);
            for i = 1:length(bslshuff.kern)
                if ~isempty(bslshuff.kern{i})
                    baseline_aucs(i) = trapz(bslshuff.kern{i}(1,:), bslshuff.kern{i}(2,:));
                end
            end
        end
        mean_baseline_auc = mean(baseline_aucs);
        auc_corrected = auc_response - mean_baseline_auc;
        fprintf('Response AUC: %.4f, Mean Shuffled Baseline AUC: %.4f, Corrected AUC: %.4f\n', auc_response, mean_baseline_auc, auc_corrected);

        % --- FIND RMAX AND HALF-RMAX ---
        [Rmax, Rmax_idx] = max(response_psth_smooth);
        Rmax_time_s = tv_s(Rmax_idx);
        halfRmax = Rmax / 2;
        halfRmax_idx = find(response_psth_smooth(Rmax_idx:end) <= halfRmax, 1, 'first');
        if ~isempty(halfRmax_idx)
            halfRmax_idx = halfRmax_idx + Rmax_idx - 1;
            halfRmax_time_s = tv_s(halfRmax_idx);
            halfRmax_val = response_psth_smooth(halfRmax_idx);
        else
            halfRmax_time_s = NaN;
            halfRmax_val = NaN;
        end

        % --- CONSOLIDATED 2x2 PLOT GENERATION ---
        fig = figure('Position',[100 100 1200 900], 'Visible', 'off');
        psth_window_ms = psth_window_s * 1000;
        sgtitle(sprintf('Channel %d | h_max: %.2f (a=%d, b=%d, c=%d) | Corrected AUC: %.3f', channel_id, hcoeffs(1), hcoeffs2D(1,1), hcoeffs2D(1,2), hcoeffs2D(1,3), auc_corrected), 'FontWeight', 'bold');

        ax1 = subplot(2,2,1); hold on;
        for trial_idx = 1:length(spikeTimes_byEvent), trial_spikes_s = spikeTimes_byEvent{trial_idx}; if ~isempty(trial_spikes_s), plot(trial_spikes_s * 1000, trial_idx * ones(size(trial_spikes_s)), 'r.', 'MarkerSize', 5); end; end
        hold off; set(gca, 'YDir', 'reverse'); xlim(psth_window_ms); ylim([0 length(stimTimes)+1]); ylabel('Trial Number'); title('Spike Raster (Response)'); grid on;

        ax2 = subplot(2,2,2); hold on;
        if ~isempty(bslshuff) && isfield(bslshuff, 'kern')
            for i = 1:length(bslshuff.kern)
                if ~isempty(bslshuff.kern{i}), plot(bslshuff.kern{i}(1,:)*1000, bslshuff.kern{i}(2,:), 'Color', [0.8 0.8 0.8], 'HandleVisibility', 'off'); end
            end
        end
        p1_diag = plot(tv_s*1000, response_psth_smooth, 'r-', 'LineWidth', 2);
        p2_diag = plot(NaN,NaN,'Color',[0.8 0.8 0.8],'LineWidth',2); % Dummy for legend
        p_mfr = yline(mfr, 'k:', 'LineWidth', 1.5);
        text(psth_window_ms(1), mfr, sprintf(' MFR: %.2f Hz', mfr), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Color', 'k', 'FontWeight','bold');
        hold off; title('H-Coefficient Diagnostic Plot'); ylabel('Firing Rate (spikes/s)'); xlabel('Time from stimulus (ms)'); legend([p1_diag, p2_diag, p_mfr], 'Response', 'Shuffled Baselines', 'Mean Firing Rate', 'Location', 'Best'); grid on;

        ax3 = subplot(2,2,3); hold on;
        edges_s = psth_window_s(1):psth_bin_width_s:psth_window_s(2);
        b = histc(psth_samples, edges_s);
        bar(edges_s * 1000, b / (num_trials * psth_bin_width_s), 1, 'FaceColor',.7*[1 1 1],'EdgeColor',.8*[1 1 1], 'HandleVisibility','off');
        plot(tv_s * 1000, response_psth_smooth, 'r-', 'LineWidth', 2, 'DisplayName', 'Smoothed PSTH');
        plot(Rmax_time_s*1000, Rmax, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 7, 'DisplayName', 'R_{max}');
        text(Rmax_time_s*1000, Rmax, sprintf(' R_{max}: %.1f Hz', Rmax), 'VerticalAlignment', 'bottom', 'HorizontalAlignment','left', 'Color', 'k', 'FontWeight','bold');
        if ~isnan(halfRmax_time_s)
            plot(halfRmax_time_s*1000, halfRmax_val, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 7, 'DisplayName', 'Half R_{max}');
            text(halfRmax_time_s*1000, halfRmax_val, sprintf(' Half R_{max} @ %.1f ms', halfRmax_time_s*1000), 'VerticalAlignment', 'top', 'HorizontalAlignment','left', 'Color', 'k', 'FontWeight','bold');
        end
        hold off; ylabel('Firing Rate (spikes/s)'); xlabel('Time from stimulus (ms)'); title('Smoothed Response PSTH & Metrics'); legend('Location', 'northeast'); grid on;

        ax4 = subplot(2,2,4);
        plot(tv_s * 1000, optw_variable_s * 1000, 'b-', 'LineWidth', 2);
        ylabel('Bandwidth (ms)'); xlabel('Time from stimulus (ms)'); title('Adaptive Kernel Bandwidth (Response)'); grid on;

        linkaxes([ax1, ax2, ax3, ax4], 'x'); xlim(psth_window_ms);
        plot_filename = fullfile(outputDir, sprintf('Analysis_Plot_Chan_%d.png', channel_id));
        saveas(fig, plot_filename); fprintf('Saved consolidated plot to %s\n', plot_filename); close(fig);

        % --- SAVE RESULTS ---
        valid_channel_count = valid_channel_count + 1;
        networkResponse(valid_channel_count).channel_id = channel_id;
        networkResponse(valid_channel_count).file_index = file_idx;
        networkResponse(valid_channel_count).mfr = mfr;
        networkResponse(valid_channel_count).auc_response = auc_response;
        networkResponse(valid_channel_count).auc_baseline_mean = mean_baseline_auc;
        networkResponse(valid_channel_count).auc_corrected = auc_corrected;
        networkResponse(valid_channel_count).peak_firing_rate_hz = Rmax;
        networkResponse(valid_channel_count).peak_time_ms = Rmax_time_s * 1000;
        networkResponse(valid_channel_count).halfRmax_time_ms = halfRmax_time_s * 1000;
        networkResponse(valid_channel_count).h_coeff_max = hcoeffs(1);
        networkResponse(valid_channel_count).h_coeff_min = hcoeffs(2);
        networkResponse(valid_channel_count).h_coeff_2D = hcoeffs2D;
    end

    % --- SAVE AGGREGATE RESULTS TO OUTPUT DIRECTORY ---
    if ~isempty(networkResponse)
        output_filename = fullfile(outputDir, sprintf('networkResponse_all_channels_%s.mat', timestamp));
        save(output_filename, 'networkResponse');
        fprintf('\nSaved aggregate network response metrics for %d channels to %s\n', valid_channel_count, output_filename);
    else, fprintf('\nNo valid channels with spikes were found to save for this file pair.\n'); end

    fprintf('\n--- Finished analysis for: %s ---\n', spikeFile);
end

fprintf('\n\nAll file pairs have been processed.\n');