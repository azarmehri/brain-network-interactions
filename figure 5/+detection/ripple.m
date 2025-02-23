function [event, data] = ripple(ch_data, fs, filter_lower_band, filter_upper_band, positive_peak, min_duration)
% RIPPLE Detects ripple events in neural signal data.
%
% This function applies an FIR filter to the input channel data, computes the
% analytic signal to obtain its amplitude envelope, and then uses threshold
% crossings on a smoothed power signal to identify potential ripple events.
% It further discards events that are too short and determines the peak index 
% within each event.
%
% Author: AmirMohammad Azarmehri
%
% Inputs:
%   ch_data           - Raw neural signal data (vector).
%   fs                - Sampling frequency in Hz.
%   filter_lower_band - Lower cutoff frequency for filtering (Hz).
%                       (Set >0 for bandpass/highpass filtering.)
%   filter_upper_band - Upper cutoff frequency for filtering (Hz).
%                       (Set >0 for bandpass/lowpass filtering.)
%   positive_peak     - Multiplier for the standard deviation added to the
%                       median power to set the detection threshold.
%   min_duration      - Minimum duration (in milliseconds) for a valid event.
%
% Outputs:
%   event - An Event object containing the start, stop, and peak indices of 
%           detected ripple events.
%   data  - The filtered neural signal data (converted to single precision).

    % Apply FIR filtering to the channel data with the specified cutoff frequencies.
    data = filter.fir(ch_data, fs, filter_lower_band, filter_upper_band);

    % Compute the amplitude envelope of the filtered signal using the Hilbert transform.
    amp = abs(hilbert(data));

    % Smooth the amplitude envelope twice using a moving average filter.
    % The smoothing window length is set to 0.05 * fs samples.
    power = smooth(smooth(amp), (0.05 * fs), 'moving');

    % Define a threshold for ripple detection:
    % threshold = median(power) + (positive_peak * standard deviation of power)
    threshold = median(power) + positive_peak * std(power);

    % Identify threshold crossings:
    % A positive crossing indicates the start of an event, and a negative crossing indicates the stop.
    crossings = diff(power > threshold);
    start = find(crossings > 0);  % Indices where the signal rises above the threshold.
    stop = find(crossings < 0);   % Indices where the signal falls below the threshold.

    % Ensure that the first stop occurs after the first start.
    if stop(1) < start(1)
        stop(1) = [];
    end

    % Ensure that the last start occurs before the last stop.
    if start(end) > stop(end)
        start(end) = [];
    end

    % Discard events that are too short.
    % Convert minimum duration from milliseconds to samples and compare with event duration.
    tooShort = stop - start < round(min_duration / 1000 * fs);
    start(tooShort) = [];
    stop(tooShort) = [];

    % Determine the peak index within each detected event.
    peak = zeros(length(start), 1);
    for i = 1:length(start)
        % Find the maximum value (peak) in the segment between start and stop.
        [~, p_ind] = max(data(start(i):stop(i)));
        % Adjust the index relative to the entire data vector.
        peak(i) = p_ind + start(i) - 1;
    end

    % Convert the filtered data to single precision.
    data = single(data);

    % Create an Event object containing the start, stop, and peak indices.
    event = detection.event(start, stop, peak);
end
