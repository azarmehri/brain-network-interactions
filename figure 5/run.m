% -------------------------------------------------------------------------
% Author:  AmirMohammad Azarmehri
% Date:    24-Feb-2025
%
% Description:
%   This script sets up and runs a neural mass model simulation by calling 
%   the function "brain_network.m". The model describes interacting brain 
%   regions (e.g., CA3, CA1, cortex, and thalamic nuclei). The simulation 
%   results are then processed to detect and characterize events such as 
%   ripples, spindles, and slow oscillations. A series of plots are 
%   generated to illustrate the dynamics and alignments of these events.
%
% Dependencies:
%   - brain_network.m (the ODE system defining the neural mass model)
%   - "library/FMA", "library/CircStat", "library/FieldTrip", 
%     "library/FieldTrip/external/eeglab", "library/EEGLAB/functions" 
%     for various detection and analysis functions
%   - detection.ripple(), detection.spindle(), detection.so() for event 
%     detection
%   - filter.fir(), analysis.pac(), and plotting functions
%
% Output:
%   - Figures illustrating detected events and average signals time-locked 
%     to those events.
%   - Printed statistics on the console (e.g., median delays, p-values).
%
% -------------------------------------------------------------------------

%% Initialize the environment and add required paths
clc
clear
close all

addpath(genpath("library/FMA"));
addpath(genpath("library/CircStat"));
addpath('library/FieldTrip', 'library/FieldTrip/external/eeglab')
addpath(genpath('library/EEGLAB/functions'))

%% Define simulation parameters
fs = 1000;             % Sampling frequency (Hz)
duration = 1000;       % Simulation duration (s)
tspan = linspace(0, duration, duration*fs+1);  % Time vector

% Initial conditions for the 17 state variables
y0 = zeros(17, 1);

% ODE solver options
options = odeset("RelTol", 1e-5, "AbsTol", ones(1,17)*1e-5, "MaxStep", duration/10);

%% Solve the ODE system
[t, y] = ode45(@brain_network, tspan, y0, options);

%% Remove transient (first second) of data
y(1:fs, :) = [];
t(1:fs, :) = [];

%% Extract indices of interest
CA1_e = 3;
CX_e = 7;
CXp_e = 10;
RE_e = 17;

%% Assign signals for convenience
mpfc_data = y(:, CX_e)';
mpfcp_data = y(:, CXp_e)';
re_data = y(:, RE_e)';
hipp_data = y(:, CA1_e)';

%% Detect ripple events in hippocampal signal
[ripple_event, ripple_data] = detection.ripple(hipp_data, fs, 150, 200, 4, 30);

%% Detect spindle events in different signals
[re_spindle_event, re_spindle_data] = detection.spindle(re_data, fs, 7, 15, 2.2, 500, 5000);
[mpfc_spindle_event, mpfc_spindle_data] = detection.spindle(mpfc_data, fs, 7, 15, 2.2, 500, 5000);
[mpfcp_spindle_event, mpfcp_spindle_data] = detection.spindle(mpfcp_data, fs, 7, 15, 2.2, 500, 5000);

%% Detect slow oscillation (SO) events
[~, mpfc_so_data] = detection.so(-mpfc_data, fs, 0.5, 2, 2, 1, 500, 2000);
[mpfcp_so_event, mpfcp_so_data] = detection.so(-mpfcp_data, fs, 0.5, 2, 2, 1, 500, 2000);
[re_so_event, re_so_data] = detection.so(-re_data, fs, 0.5, 2, 2, 1, 500, 2000);

% Restore sign (because detection was run on inverted data)
mpfc_so_data  = -mpfc_so_data;
mpfcp_so_data = -mpfcp_so_data;
re_so_data    = -re_so_data;

%% (Figure 5C) Plot histogram of RE_so_event vs. mpfcp_so_event
plot_hist(re_so_event.peaks, mpfcp_so_event.peaks, fs, [-.5 .5], 5, 0.02, [0 100], ...
    "Time to mPFC slow wave (s)", "C")

%% (Figure 5D) Plot histogram of RE spindle onset vs. mPFC spindle onset
plot_hist(re_spindle_event.starts, mpfcp_spindle_event.starts, fs, [-2 2], 5, 0.1, [0 40], ...
    "Time to mPFC spindle onset (s)", "D")

%% (Figure 5E & 5F) Time-frequency representation (TFR) around SO events
plot_tfr(mpfcp_data, re_data, mpfcp_so_event, fs)

%% (Figure 5G) Compare average signals time-locked to slow oscillation peaks
hipp_data_filter = filter.fir(hipp_data, fs, 0.5, 2);

figure("Units","centimeters","Position",[1.5 1.5 15 6]);
plot_mean_all(mpfcp_so_data, "CXp_e", ...
    re_so_data, "REU_e", hipp_data_filter, mpfcp_so_event.peaks, fs, [-1 1], "SO ref.");

%% (Figure 5H) Circular histogram of ripple peak phases relative to mPFC slow wave
ph = angle(hilbert(mpfcp_so_data));
phases = ph(ripple_event.peaks);
figure("Units","inches","Position",[1 1 4 4]);
plots.circ_hist(phases, "k", 14);

%% (Figure 5I) Plot mean signals of mPFCp, RE, hippocampus time-locked to ripple peaks
fig3 = figure("Units","centimeters","Position",[1.5 1.5 15 19]);
plot_mean(mpfcp_data, "CXp_e", ...
    re_data, "REU_e", hipp_data, ripple_event.peaks, fs, [-1 1], "Ripple peak");

%% (Figure 5J) Plot mean signals time-locked to mPFC spindle onset and ripple events
ripples = false(1, length(hipp_data));
ripples(ripple_event.peaks) = 1;

figure("Units","centimeters","Position",[1.5 1.5 15 19]);
plot_mean1(mpfcp_data, "CXp_e", ...
    re_data, "REU_e", ripples, mpfc_spindle_event.starts, fs, [-1 1], "CX spindle onset");


%% ------------------------------------------------------------------------
%  Below are the local helper functions for data alignment, averaging, and
%  plotting. They are organized to facilitate analysis and visualization of
%  detected events, histograms, time-frequency analyses, and mean responses.
% ------------------------------------------------------------------------

function plot_hist(re_spindle_event, mpfcp_spindle_event, fs, range, count, bin, yrange, x_label, name)
    % PLOT_HIST
    %   Plots a histogram of time delays between two sets of events.
    %
    %   Inputs:
    %       re_spindle_event    - Timestamps of events in the RE
    %       mpfcp_spindle_event - Timestamps of events in mPFC
    %       fs                  - Sampling frequency
    %       range               - x-axis range to plot
    %       count               - Number of tick marks for x-axis
    %       bin                 - Histogram bin width
    %       yrange              - y-axis limits
    %       x_label            - X-axis label
    %       name               - Name or label for printed output
    %
    mpfc_delay = calc_nearest_events_delay(re_spindle_event/fs, mpfcp_spindle_event/fs);
    mpfc_delay(abs(mpfc_delay) > range(2)) = [];

    disp_range_ticks = linspace(range(1),range(2),count);
    disp_range_lables = string(disp_range_ticks);
    bins = round(diff(range)/bin);

    f = figure("Units","centimeters", "Position", [3 3 10 5]);
    nt = nexttile();

    histogram(mpfc_delay, bins, "BinLimits", range, "FaceColor", "k");
    xline(0, 'r', 'LineWidth', 1.5)

    med = median(mpfc_delay);
    [p,h,stats] = signrank(mpfc_delay);

    annotation("arrow", "Position", [(med-range(1))/diff(range) nt.Position(2)+nt.Position(4) 0 -.07], ...
        "Color", "k", "LineWidth", 1, "HeadStyle", "cback1", "HeadLength", 5, "HeadWidth", 5);
    annotation("textbox", "Position", [nt.Position(1) nt.Position(2)+nt.Position(4) ...
        (med-range(1))/diff(range)-nt.Position(1)+0.025 0.1], ...
        "String", "REU", "Color", "k", "FontSize", 12, "EdgeColor", "none", ...
        "Margin", 0, "HorizontalAlignment", "right", "VerticalAlignment", "bottom");

    xticks(disp_range_ticks);
    xticklabels(disp_range_lables);
    xlim(range);
    ylim(yrange);
    xlabel(x_label);
    ylabel("N")

    nt.Box="off";
    nt.FontSize=12;

    disp(name)
    fprintf("n = %d\n", length(mpfc_delay))
    fprintf("median = %d, mad = %d\n", med, mad(mpfc_delay))
    fprintf("re: p = %d, z = %d\n", p, stats.zval)
end


function delay = calc_nearest_events_delay(event1, event2)
    % CALC_NEAREST_EVENTS_DELAY
    %   Computes the time difference between events in event1 and the
    %   nearest event in event2.
    %
    %   Inputs:
    %       event1 - Vector of event timestamps
    %       event2 - Vector of event timestamps
    %
    %   Output:
    %       delay  - Vector of delays (event1 minus closest event2)
    %
    [~, ind] = min(abs(event1 - event2'));
    delay = event1(ind) - event2;
end


function [res, t] = mean_data(data, event, fs, range)
    % MEAN_DATA
    %   Extracts time-locked segments around each event from 'data' and
    %   computes the average.
    %
    %   Inputs:
    %       data  - Signal vector
    %       event - Vector of event indices
    %       fs    - Sampling frequency
    %       range - [start end] time window in seconds
    %
    %   Outputs:
    %       res - The mean time-locked signal
    %       t   - Time vector (seconds)
    %
    ind = (range(1)*fs : range(2)*fs);
    res = zeros(size(ind));

    event(event < range(2)*fs) = [];
    event(event > length(data) + range(1)*fs) = [];

    for i = 1:length(event)
        res = res + (data(ind + event(i)) / length(event));
    end

    t = ind / fs;
end


function plot_mean(data1, title1, data2, title2, hipp_data, event, fs, range, fig_title)
    % PLOT_MEAN
    %   Plots the mean signals of three data sets (data1, data2, and 
    %   hipp_data) time-locked to the same set of events.
    %
    [data1_mean, time1] = mean_data(data1, event, fs, range);
    [data2_mean, time2] = mean_data(data2, event, fs, range);
    [data3_mean, time3] = mean_data(hipp_data, event, fs, range);

    t = tiledlayout(3, 1, "Padding", "compact");
    title(t, fig_title)

    nexttile;
    plot(time1, data1_mean, "k", "LineWidth", 2)
    ylabel(title1);
    xticks((range(1):0.5:range(2)));
    pax = gca;
    pax.FontSize = 12;

    nexttile;
    plot(time2, data2_mean, "k", "LineWidth", 2)
    ylabel(title2);
    xticks((range(1):0.5:range(2)));
    pax = gca;
    pax.FontSize = 12;

    nexttile;
    plot(time3, data3_mean, "k", "LineWidth", 2)
    ylabel("CA1_e");
    xticks((range(1):0.5:range(2)));
    pax = gca;
    pax.FontSize = 12;
end


function plot_mean_all(data1, title1, data2, title2, hipp_data, event, fs, range, fig_title)
    % PLOT_MEAN_ALL
    %   Plots mean signals (data1, data2, hipp_data) in a single axes, 
    %   time-locked to 'event'.
    %
    [data1_mean, time1] = mean_data(data1, event, fs, range);
    [data2_mean, time2] = mean_data(data2, event, fs, range);
    [data3_mean, time3] = mean_data(hipp_data, event, fs, range);

    title(fig_title)
    plot(time1, data1_mean, "LineWidth", 2)
    hold on
    plot(time2, data2_mean, "LineWidth", 2)
    plot(time3, data3_mean, "LineWidth", 2)

    legend(title1, title2, "CA1_e", "Location", "best")
    xticks((range(1):0.5:range(2)));
    pax = gca;
    pax.FontSize = 12;
end


function plot_mean1(data1, title1, data2, title2, ripples, event, fs, range, fig_title)
    % PLOT_MEAN1
    %   Similar to 'plot_mean' but also plots ripple probability in a third
    %   tile as a histogram of ripple events, aligned to 'event'.
    %
    [data1_mean, time1] = mean_data(data1, event, fs, range);
    [data2_mean, time2] = mean_data(data2, event, fs, range);

    t = tiledlayout(3, 1, "Padding", "compact");
    title(t, fig_title)

    nexttile;
    plot(time1, data1_mean, "k", "LineWidth", 2)
    ylabel(title1);
    xticks((range(1):0.5:range(2)));
    pax = gca;
    pax.FontSize = 12;

    nexttile;
    plot(time2, data2_mean, "k", "LineWidth", 2)
    ylabel(title2);
    xticks((range(1):0.5:range(2)));
    pax = gca;
    pax.FontSize = 12;

    nexttile;
    plots.hist_avg(event, ripples, 1, 0.05, fs)
    xlim(range)
    ylabel("Ripples probability")
    xlabel("Time (s)")
    xticks((range(1):0.5:range(2)));
    pax = gca;
    pax.FontSize = 12;
end


function plot_tfr(mpfc_data, re_data, so_event, fs)
    % PLOT_TFR
    %   Computes and plots time-frequency representations (TFR) and 
    %   phase-amplitude coupling (PAC) for mPFC and reunions data, time-locked 
    %   to slow oscillation (SO) events. Also plots histograms of coupling 
    %   strengths.
    %
    raw_freq = [7 15];    % Spindle band
    event_freq = [0.5 2]; % SO band

    [si_mpfc, mean_data_mpfc, ~, mean_tfr_mpfc, freq_mpfc, time_mpfc] = ...
        analysis.pac(mpfc_data, raw_freq, mpfc_data, event_freq, so_event.peaks, fs, 2, [-3.5 -2.5], false);
    [si_re, mean_data_re, ~, mean_tfr_re, freq_re, time_re] = ...
        analysis.pac(re_data, raw_freq, mpfc_data, event_freq, so_event.peaks, fs, 2, [-3.5 -2.5], false);

    label = 'Power';

    figure("Units","centimeters","Position",[1.5 1.5 30 19]);

    % -------------------- mPFC TFR --------------------
    fp = [0 2/3 .5 1/3];
    padding = [2.5 5 5 5]/100;  % [top right bottom left]
    aa = fp + [padding(4) padding(3) -padding(2)-padding(4) -padding(1)-padding(3)];

    annotation("textbox", "Position", [fp(1) fp(2)+fp(4)-padding(1) padding(4) padding(1)], ...
        "String", "E", "FontSize", 14, "FontWeight", "bold", "EdgeColor", "none", ...
        "Margin", 0, "HorizontalAlignment", "left", "VerticalAlignment", "top");
    annotation("textbox", "Position", [aa(1) fp(2)+fp(4)-padding(1) aa(3) padding(1)], ...
        "String", "mPFC PAC", "FontSize", 12, "EdgeColor", "none", "Margin", 0, ...
        "HorizontalAlignment", "center", "VerticalAlignment", "bottom");

    ax = axes("Position", aa);
    imagesc(time_mpfc, freq_mpfc, mean_tfr_mpfc)
    set(gca,'YDir', 'normal');
    colormap("jet");
    cb = colorbar;
    ylabel(cb,label);
    hold on
    plot(time_mpfc, mean_data_mpfc/max(abs(mean_data_mpfc))*(freq_mpfc(end)-freq_mpfc(1))/2+mean(freq_mpfc), ...
        "LineWidth", 2, "Color", "k")
    ylabel("Frequency (Hz)")
    xticks([-2 -1 0 1 2]);
    xticklabels(["-2", "-1", "0", "1", "2 s"]);
    xlim([-2 2]);
    ax.Box="off";
    ax.FontSize=12;

    % mPFC phase histogram
    fp = [.5 2/3+.07 .21 .21];
    padding = [0 0 0 0]/100;
    aa = fp + [padding(4) padding(3) -padding(2)-padding(4) -padding(1)-padding(3)];
    axes("Position", aa);
    plots.circ_hist(angle(si_mpfc), "r", 12);

    % -------------------- Reuniens TFR --------------------
    fp = [0 1/3 .5 1/3];
    padding = [2.5 5 5 5]/100;
    aa = fp + [padding(4) padding(3) -padding(2)-padding(4) -padding(1)-padding(3)];

    annotation("textbox", "Position", [aa(1) fp(2)+fp(4)-padding(1) aa(3) padding(1)], ...
        "String", "Reuniens PAC", "FontSize", 12, "EdgeColor", "none", "Margin", 0, ...
        "HorizontalAlignment", "center", "VerticalAlignment", "bottom");

    ax = axes("Position", aa);
    imagesc(time_re, freq_re, mean_tfr_re)
    set(gca,'YDir', 'normal');
    colormap("jet");
    cb = colorbar;
    ylabel(cb,label);
    hold on
    plot(time_re, mean_data_re/max(abs(mean_data_re))*(freq_re(end)-freq_re(1))/2+mean(freq_re), ...
        "LineWidth", 2, "Color", "k")
    ylabel("Frequency (Hz)")
    xticks([-2 -1 0 1 2]);
    xticklabels(["-2", "-1", "0", "1", "2 s"]);
    xlim([-2 2]);
    ax.Box="off";
    ax.FontSize=12;

    % Reuniens phase histogram
    fp = [.5 1/3+.07 .21 .21];
    padding = [0 0 0 0]/100;
    aa = fp + [padding(4) padding(3) -padding(2)-padding(4) -padding(1)-padding(3)];
    axes("Position", aa);
    plots.circ_hist(angle(si_re), "k", 12);

    % -------------------- Coupling strength histogram --------------------
    fp = [0 0 .5 1/3];
    padding = [2.5 5 5 5]/100;
    aa = fp + [padding(4) padding(3) -padding(2)-padding(4) -padding(1)-padding(3)];

    annotation("textbox", "Position", [fp(1) fp(2)+fp(4)-padding(1) padding(4) padding(1)], ...
        "String", "F", "FontSize", 14, "FontWeight", "bold", "EdgeColor", "none", ...
        "Margin", 0, "HorizontalAlignment", "left", "VerticalAlignment", "top");
    annotation("textbox", "Position", [aa(1) fp(2)+fp(4)-padding(1) aa(3) padding(1)], ...
        "String", "SO-Spindle coupling strength", "FontSize", 12, "EdgeColor", "none", ...
        "Margin", 0, "HorizontalAlignment", "center", "VerticalAlignment", "bottom");

    ax = axes("Position", aa);
    hold on;
    h1 = histogram(abs(si_re), 20, "FaceColor", "k", "EdgeColor", "k", "Normalization", "probability");
    h2 = histogram(abs(si_mpfc), 20, "FaceColor", "r", "EdgeColor", "r", "Normalization", "probability");
    xline(median(abs(si_mpfc)), "Color", "r", "LineWidth", 2);
    xline(median(abs(si_re)), "Color", "k", "LineWidth", 2);
    ylabel("P");
    legend([h1 h2], ["Reuniens","mPFC"], "Location", "best")

    ax.Box="off";
    ax.FontSize=12;

    fprintf("mpfc: median = %.2f, mad = %.2f\n", median(abs(si_mpfc)), mad(abs(si_mpfc), 1))
    fprintf("re: median = %.2f, mad = %.2f\n", median(abs(si_re)), mad(abs(si_re), 1))
end
