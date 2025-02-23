function filtered_data = fir(data, fs, filter_lower_band, filter_upper_band)
% FIR Filter Function from EEGLAB Library
% This function applies a finite impulse response (FIR) filter to the input data.
% It can perform bandpass, highpass, or lowpass filtering based on the provided cutoff frequencies.
%
% Inputs:
%   data              - The input signal (vector or matrix) to be filtered.
%   fs                - Sampling frequency (in Hz).
%   filter_lower_band - Lower cutoff frequency (Hz); set >0 for bandpass or highpass filtering.
%   filter_upper_band - Upper cutoff frequency (Hz); set >0 for bandpass or lowpass filtering.
%
% Output:
%   filtered_data     - The filtered signal.
%
% The function determines an appropriate filter order based on the sampling frequency 
% and the cutoff frequencies, designs the FIR filter using fir1, and applies the filter 
% with zero-phase distortion using filtfilt.

    % Calculate the Nyquist frequency (half of the sampling frequency)
    nyq = fs * 0.5;
    
    % Minimum factor used in filter order calculation (ensures sufficient order)
    minfac = 3;
    
    % Minimum filter order (not directly used here but can serve as a guideline)
    min_filtorder = 15;
    
    % Transition band width factor (used for validating cutoff frequencies)
    trans = 0.15;

    % Determine the filter order based on the provided cutoff frequency.
    % The order is scaled with the sampling frequency and the relevant cutoff.
    if filter_lower_band > 0
        filtorder = minfac * fix(fs / filter_lower_band);
    elseif filter_upper_band > 0
        filtorder = minfac * fix(fs / filter_upper_band);
    end

    % Validate that the high cutoff frequency is not too close to the Nyquist frequency.
    % The check accounts for an additional transition band.
    if (1 + trans) * filter_upper_band / nyq > 1
        error('High cutoff frequency too close to Nyquist frequency');
    end

    % Design the FIR filter based on the type of filtering required:
    % - Bandpass: both lower and upper cutoff frequencies are provided.
    % - Highpass: only the lower cutoff frequency is provided.
    % - Lowpass: only the upper cutoff frequency is provided.
    if filter_lower_band > 0 && filter_upper_band > 0    % Bandpass filter
        filtwts = fir1(filtorder, [filter_lower_band, filter_upper_band] ./ (fs / 2));
    elseif filter_lower_band > 0                          % Highpass filter
        filtwts = fir1(filtorder, filter_lower_band ./ (fs / 2), 'high');
    elseif filter_upper_band > 0                          % Lowpass filter
        filtwts = fir1(filtorder, filter_upper_band ./ (fs / 2));
    else
        error('You must provide a non-0 low or high cut-off frequency');
    end

    % Check that the data length is sufficient for filtering.
    % The data should be at least three times longer than the filter order.
    if filtorder * 3 > length(data)
        error('Data length must be at least 3 times the filter order.');
    end

    % Apply zero-phase forward and reverse filtering using filtfilt to prevent phase distortion.
    filtered_data = filtfilt(filtwts, 1, data);
end