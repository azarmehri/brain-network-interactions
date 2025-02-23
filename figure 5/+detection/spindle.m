function [event, data] = spindle(ch_data, fs, filter_lower_band, filter_upper_band, positive_peak, min_duration, max_duration)
% SPINDLE Detects spindle events in neural signal data.
%
% This function filters the input channel data, detects spindle events based on
% duration and threshold criteria using the FMA.FindSpindles method, and then 
% extracts the start, stop, and peak indices of the detected spindles.
%
% Author: AmirMohammad Azarmehri
%
% Inputs:
%   ch_data           - Raw neural signal data (vector).
%   fs                - Sampling frequency (Hz).
%   filter_lower_band - Lower cutoff frequency for filtering (Hz).
%   filter_upper_band - Upper cutoff frequency for filtering (Hz).
%   positive_peak     - Threshold multiplier for spindle detection.
%   min_duration      - Minimum duration (in ms) for a valid spindle event.
%   max_duration      - Maximum duration (in ms) for a valid spindle event.
%
% Outputs:
%   event - An object containing the start, stop, and peak indices of detected spindles.
%   data  - The filtered neural signal data, converted to single precision.

    % Apply FIR filtering to the input channel data using specified cutoff frequencies.
    % The last parameter (0) may indicate a specific filtering mode or flag.
    data = filter.fir(ch_data, fs, filter_lower_band, filter_upper_band);
    
    % Create a time vector corresponding to the filtered data.
    time = (1:length(data)) / fs;
    
    % Combine time and data into a two-column matrix for spindle detection.
    % Each row: [time, data_value]
    filtered = [time.' data.'];
    
    % Attempt to detect spindles using the FMA.FindSpindles function.
    % The function is provided with duration limits, a threshold, and peak detection flag.
    try
        spindles = FindSpindles(filtered, ...
            'durations', [min_duration max_duration], ...
            'threshold', positive_peak, 'peak', 1);
    catch
        % If an error occurs during detection, set spindles to empty.
        spindles = [];
    end
    
    % If no spindles are detected, return an empty event and the filtered data.
    if isempty(spindles)
        event = [];
        data = single(data);
        return
    end
    
    % Convert detected spindle times to sample indices.
    % Assumes spindles(:,1) contains start times and spindles(:,3) contains stop times.
    start = round(spindles(:,1) * fs);
    stop = round(spindles(:,3) * fs);
    
    % Ensure the stop index of the last event does not exceed the length of the data.
    if stop(end) > length(data)
        stop(end) = length(data);
    end
    
    % Find the peak index within each detected spindle event.
    peak = zeros(length(start), 1);
    for i = 1:length(start)
        % Find the index of the maximum value in the data segment between start and stop.
        [~, p_ind] = max(data(start(i):stop(i)));
        % Adjust the peak index relative to the entire data vector.
        peak(i) = p_ind + start(i) - 1;
    end
    
    % Convert the filtered data to single precision.
    data = single(data);
    
    % Create an event object using the detected start, stop, and peak indices.
    event = detection.event(start, stop, peak);
end
