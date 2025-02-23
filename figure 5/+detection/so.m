function [event, data] = so(ch_data, fs, filter_lower_band, filter_upper_band, positive_peak, negetive_peak, min_duration, max_duration)
% SO Detects delta wave events in neural signal data.
%
% This function filters the input channel data to extract delta frequency components 
% and detects delta waves using a modified version of the FMA.FindDeltaWaves function.
% The FindDeltaWaves function is copied from the FMA toolbox and has been corrected 
% for use in my analysis.
%
% Author: AmirMohammad Azarmehri
%
% Inputs:
%   ch_data           - Raw neural signal data (vector).
%   fs                - Sampling frequency (Hz).
%   filter_lower_band - Lower cutoff frequency for primary filtering (Hz).
%   filter_upper_band - Upper cutoff frequency for primary filtering (Hz).
%   positive_peak     - Threshold for positive peaks.
%   negetive_peak     - Threshold for negative peaks.
%   min_duration      - Minimum duration (in ms) for a valid delta wave event.
%   max_duration      - Maximum duration (in ms) for a valid delta wave event.
%
% Outputs:
%   event - An object containing the start, stop, and peak indices of detected delta waves.
%   data  - The filtered neural signal data (converted to single precision).

    % Apply FIR filtering to extract the primary frequency band of interest.
    % The last parameter (0) may indicate a specific filtering mode or flag.
    data = filter.fir(ch_data, fs, filter_lower_band, filter_upper_band);
    
    % Apply FIR filtering to extract the "fast" component (delta band: 0.5-4 Hz).
    data_fast = filter.fir(ch_data, fs, 0.5, 4);
    
    % Create a time vector corresponding to the filtered data.
    time = (1:length(data)) / fs;
    
    % Combine time and data into a two-column matrix for delta wave detection.
    % Each row consists of [time, data_value].
    filtered = [time.' data.'];
    filtered_fast = [time.' data_fast.'];
    
    % Detect delta waves using the modified FindDeltaWaves function.
    % The function uses both the primary filtered data and the fast (delta band) component.
    % Thresholds for positive and negative peaks and duration limits (in ms) are specified.
    deltas = FindDeltaWaves(filtered, filtered_fast, fs, ...
        'thresholds', [positive_peak positive_peak negetive_peak negetive_peak], ...
        'durations', [min_duration max_duration]);
    
    % Convert detected delta wave times to sample indices.
    start = round(deltas(:,1) * fs);
    stop = round(deltas(:,3) * fs);
    
    % Find the peak index within each detected delta wave event.
    peak = zeros(length(start), 1);
    for i = 1:length(start)
        % Identify the maximum value in the data segment between start and stop.
        [~, p_ind] = max(data(start(i):stop(i)));
        % Adjust the index relative to the entire data vector.
        peak(i) = p_ind + start(i) - 1;
    end
    
    % Convert the filtered data to single precision.
    data = single(data);
    
    % Create an event object containing the start, stop, and peak indices.
    event = detection.event(start, stop, peak);
end
