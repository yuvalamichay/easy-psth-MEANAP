% ADD 'bslshuff' to the function definition
function [hcoeffs, hcoeffs2D, kern, bslshuff] = hcoeff_wrapper(peristim_cel, all_spikes_vec, psth_window_s)
%Hcoeff_wrapper Adapts per-channel data to be used by the original hcoeff function.

% --- Parameters for hcoeff ---
nstims = length(peristim_cel);
rsp = psth_window_s * 1000;
fast = 1;

% --- Calculate Mean Firing Rate (mfr) ---
if isempty(all_spikes_vec) || length(all_spikes_vec) < 2
    mfr = 0;
else
    total_duration_s = all_spikes_vec(end) - all_spikes_vec(1);
    mfr = length(all_spikes_vec) / total_duration_s;
end

if mfr == 0 || nstims == 0
    hcoeffs = [NaN, NaN];
    hcoeffs2D = NaN(2, 3);
    kern = [];
    bslshuff = []; % ADDED: Return empty if no data
    return;
end

% --- Prepare data for hcoeff function ---
peristim_ms = cellfun(@(x) x(:)' * 1000, peristim_cel, 'UniformOutput', false);
rawdata_ms = all_spikes_vec(:)' * 1000;

% --- Call the original hcoeff function ---
% MODIFIED: Update the call to receive bslshuff
[kern, hcoeffs, hcoeffs2D, bslshuff] = hcoeff(peristim_ms, rawdata_ms, mfr, nstims, rsp, fast);

end
