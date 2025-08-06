%% mfr version

%% BATCH ANALYSIS OF KERNEL DENSITY ESTIMATION FOR PSTH
% This script loops through all channels in a recording. For each channel,
% it generates a single 2x2 figure that includes:
% 1. A spike raster plot.
% 2. The h-coefficient diagnostic plot (response vs. shuffled baseline), with Rmax, half Rmax, and AUC markers.
% 3. The smoothed PSTH.
% 4. The adaptive kernel bandwidth over time.
% All metrics are saved into a single aggregate struct.

clear; clc;

%% --- USER PARAMETERS ---
spikeFile = 'OWT220207_1I_DIV63_HUB63_6UA_Cspikes_L0_RP2.mat_Nortefact.mat';
rawFile   = 'OWT220207_1I_DIV63_HUB63_6UA.mat';
fs = 25000;
spikeMethod = 'bior1p5';
numChannels = 60;
artifact_window_ms = [0, 2];
psth_window_s = [0, 0.02];
psth_bin_width_s = 0.001;

% --- Stimulation Detection Parameters (from batch_latency, for "longblank" method) ---
min_blanking_duration_ms = 1.5; % This corresponds to the old 'flat_window_ms'
min_interval_ms = 2500;         % This is the same.


%% --- CHANNEL REMAPPING & SETUP ---
indices = [24 26 29 32 35 37, 21 22 25 30 31 36 39 40, 19 20 23 28 33 38 41 42, 16 17 18 27 34 43 44 45, 15 14 13 4 57 48 47 46, 12 11 8 3 58 53 50 49, 10 9 6 1 60 55 52 51, 7 5 2 59 56 54];
ids = [21 31 41 51 61 71, 12 22 32 42 52 62 72 82, 13 23 33 43 53 63 73 83, 14 24 34 44 54 64 74 84, 15 25 35 45 55 65 75 85, 16 26 36 46 56 66 76 86, 17 27 37 47 57 67 77 87, 28 38 48 58 68 78];
channelMap = containers.Map('KeyType','double','ValueType','double');
for i = 1:numel(indices), channelMap(indices(i)) = ids(i); end
timestamp = datestr(now, 'ddmmmyyyy_HH:MM');
outputDir = ['PSTHanalysis(' timestamp ')'];
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

% --- Detect stimulation times using the "longblank" external function ---
fprintf('Detecting stimulation events...\n');
% Note the tildes (~) for unused parameters in the new function
[stimTimes_ms, ~, ~] = detect_stim_times(rawFile, numChannels, fs, ...
    [], min_blanking_duration_ms, [], min_interval_ms, []);

% Convert stim times to seconds for consistency with the rest of the script
stimTimes = stimTimes_ms / 1000;
stimTimes = sort(unique(stimTimes(:))); % Ensure sorted unique values

fprintf('Found %d stimulation events.\n', length(stimTimes));

%% --- STEP 2: LOOP THROUGH CHANNELS FOR ANALYSIS ---
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

    out = WithinRanges(spikeTimes_cleaned_s, stimTimes + psth_window_s, (1:length(stimTimes))', 'matrix');
    spikeTimes_byEvent = arrayfun(@(n) spikeTimes_cleaned_s(logical(out(:,n))) - stimTimes(n), 1:length(stimTimes), 'uni', 0)';
    psth_samples = cell2mat(spikeTimes_byEvent);

    if isempty(psth_samples), fprintf('No spikes in PSTH window. Skipping.\n'); continue; end
    fprintf('Processing %d spikes from Channel ID %d.\n', length(psth_samples), channel_id);

    %% --- H-COEFFICIENT CALCULATION ---
    fprintf('Calculating h-coefficient...\n');
    % MODIFICATION: Capture mfr output from wrapper
    [hcoeffs, hcoeffs2D, ~, bslshuff, mfr] = hcoeff_wrapper(spikeTimes_byEvent, all_spike_times_s, psth_window_s, psth_bin_width_s);
    fprintf('h-coeff (max/min): %.2f / %.2f\n', hcoeffs(1), hcoeffs(2));

    %% --- CONSOLIDATED 2x2 PLOT GENERATION ---
    fig = figure('Position',[100 100 1200 900], 'Visible', 'off');
    psth_window_ms = psth_window_s * 1000;
    
    % --- MODIFICATION 1: Main title for the entire figure ---
    sgtitle(sprintf('Channel %d | h_max: %.2f (a = %d, b = %d, c = %d)', channel_id, hcoeffs(1), hcoeffs2D(1,1), hcoeffs2D(1,2), hcoeffs2D(1,3)), 'FontWeight', 'bold');

    % --- Plot 1: Raster Plot ---
    ax1 = subplot(2,2,1);
    hold on;
    for trial_idx = 1:length(spikeTimes_byEvent)
        trial_spikes_s = spikeTimes_byEvent{trial_idx};
        if ~isempty(trial_spikes_s)
            plot(trial_spikes_s * 1000, trial_idx * ones(size(trial_spikes_s)), 'r.', 'MarkerSize', 5);
        end
    end
    hold off;
    set(gca, 'YDir', 'reverse'); 
    xlim(psth_window_ms);
    ylim([0 length(spikeTimes_byEvent)+1]);
    ylabel('Trial Number');
    title('Spike Raster');
    grid on;

    % --- Plot 2: H-Coefficient Diagnostic Plot ---
    ax2 = subplot(2,2,2);
    hold on;
    if ~isempty(bslshuff) && isfield(bslshuff, 'kern')
        % Plot all the shuffled baseline kernels in gray
        for i = 1:length(bslshuff.kern)
            if ~isempty(bslshuff.kern{i})
                plot(bslshuff.kern{i}(1,:)*1000, bslshuff.kern{i}(2,:), 'Color', [0.8 0.8 0.8], 'HandleVisibility', 'off');
            end
        end
    end
    % Calculate smoothed PSTH for plotting
    L = 1000; t_s = linspace(psth_window_s(1), psth_window_s(2), L);
    [yv_pdf, tv_s, optw_variable_s] = ssvkernel(psth_samples, t_s);
    num_trials = length(stimTimes);
    avg_spikes_per_trial = length(psth_samples) / num_trials;
    yv = yv_pdf * avg_spikes_per_trial;
    % Plot the actual response kernel in red
    p1_diag = plot(tv_s*1000, yv, 'r-', 'LineWidth', 2);
    p2_diag = plot(NaN,NaN,'Color',[0.8 0.8 0.8],'LineWidth',2); % Dummy plot for legend
    
    % MODIFICATION: Plot mean firing rate (mfr) as a dotted line
    p4_diag = yline(mfr, 'k:', 'LineWidth', 1.5);

    hold off;
    title('H-Coeff Diagnostic');
    ylabel('Firing Rate (spikes/s)');
    legend([p1_diag, p2_diag, p4_diag], 'Response', 'Baseline', 'Mean Firing Rate', 'Location', 'Best');
    grid on;

    % --- Plot 3: Smoothed PSTH and Metrics with Rmax, Half Rmax, AUC ---
    ax3 = subplot(2,2,3);
    hold on;
    % Raw histogram bars
    edges_s = psth_window_s(1):psth_bin_width_s:psth_window_s(2);
    b = histc(psth_samples, edges_s);
    bar(edges_s * 1000, b / (num_trials * psth_bin_width_s), 1, 'FaceColor',.7*[1 1 1],'EdgeColor',.8*[1 1 1], 'HandleVisibility','off');
    % Smoothed PSTH
    plot(tv_s * 1000, yv, 'r-', 'LineWidth', 2);

    % --- Find Global Maximum and Half-Max Following Peak ---
    [Rmax, Rmax_idx] = max(yv);
    Rmax_time_s = tv_s(Rmax_idx);
    halfRmax = Rmax / 2;
    halfRmax_idx = find(yv(Rmax_idx:end) <= halfRmax, 1, 'first');
    if ~isempty(halfRmax_idx)
        halfRmax_idx = halfRmax_idx + Rmax_idx - 1;
        halfRmax_time_s = tv_s(halfRmax_idx);
        halfRmax_val = yv(halfRmax_idx);
    else
        halfRmax_time_s = NaN;
        halfRmax_val = NaN;
    end

    % AUC calculation
    auc = trapz(tv_s, yv);

    % Plot Rmax marker
    plot(Rmax_time_s*1000, Rmax, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 7, 'DisplayName', 'R_{max}');
    text(Rmax_time_s*1000, Rmax, sprintf(' R_{max}: %.1f Hz', Rmax), 'VerticalAlignment', 'bottom', 'HorizontalAlignment','left', 'Color', 'k', 'FontWeight','bold');

    % Plot Half Rmax marker
    if ~isnan(halfRmax_time_s)
        plot(halfRmax_time_s*1000, halfRmax_val, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 7, 'DisplayName', 'Half R_{max}');
        text(halfRmax_time_s*1000, halfRmax_val, sprintf(' Half R_{max} @ %.1f ms', halfRmax_time_s*1000), 'VerticalAlignment', 'top', 'HorizontalAlignment','left', 'Color', 'k', 'FontWeight','bold');
    end

    % --- MODIFICATION 2: Show AUC on plot (right side) ---
    plot_ylim = get(gca, 'YLim');
    text(psth_window_ms(2), plot_ylim(2)*0.9, sprintf('AUC: %.3f', auc), 'FontWeight', 'bold', 'Color','k', 'HorizontalAlignment', 'right');
    
    hold off;
    ylabel('Firing Rate (spikes/s)');
    xlabel('Time from stimulus (ms)');
    title('Smoothed PSTH & Metrics (+ R_{max}, Half R_{max}, AUC)');
    legend('Smoothed PSTH', 'Location', 'Best');
    grid on;

    % --- Plot 4: Kernel Bandwidth ---
    ax4 = subplot(2,2,4);
    plot(tv_s * 1000, optw_variable_s * 1000, 'b-', 'LineWidth', 2);
    ylabel('Bandwidth (ms)');
    xlabel('Time from stimulus (ms)');
    title('Adaptive Kernel Bandwidth');
    grid on;
    
    % --- Final Figure Adjustments ---
    linkaxes([ax1, ax2, ax3, ax4], 'x');
    xlim(psth_window_ms);

    plot_filename = fullfile(outputDir, sprintf('Analysis_Plot_Chan_%d.png', channel_id));
    saveas(fig, plot_filename);
    fprintf('Saved consolidated plot to %s\n', plot_filename);
    close(fig);

  %% --- SAVE RESULTS ---
    valid_channel_count = valid_channel_count + 1;
    networkResponse(valid_channel_count).channel_id = channel_id;
    networkResponse(valid_channel_count).file_index = file_idx;
    networkResponse(valid_channel_count).auc = auc;
    networkResponse(valid_channel_count).peak_firing_rate_hz = Rmax;
    networkResponse(valid_channel_count).peak_time_ms = Rmax_time_s * 1000;
    networkResponse(valid_channel_count).halfRmax_time_ms = halfRmax_time_s * 1000;
    networkResponse(valid_channel_count).h_coeff_max = hcoeffs(1);
    networkResponse(valid_channel_count).h_coeff_min = hcoeffs(2);
    networkResponse(valid_channel_count).h_coeff_2D = hcoeffs2D;
    networkResponse(valid_channel_count).response_duration_ms = halfRmax_time_s * 1000; % duration after Rmax at half-max
    networkResponse(valid_channel_count).mfr = mfr; % Also save mfr to results
end

%% --- SAVE AGGREGATE RESULTS TO OUTPUT DIRECTORY ---
if ~isempty(networkResponse)
    output_filename = fullfile(outputDir, sprintf('networkResponse_all_channels_%s.mat', timestamp));
    save(output_filename, 'networkResponse');
    fprintf('\nSaved aggregate network response metrics for %d channels to %s\n', valid_channel_count, output_filename);
else, fprintf('\nNo valid channels with spikes were found to save.\n'); end
fprintf('Script finished.\n');
