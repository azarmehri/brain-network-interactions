
function [si, meanData, meanRawData, meanTFR, freq, time] = pac(raw_data, raw_freq, event_data, event_freq, event, fs, window, pre_event, baseline)
win_total = 4;
padding_window = 9;
pre_event_start = pre_event(1);
pre_event_stop = pre_event(2);


% remove events
wt = floor(win_total*fs);
event(event <= wt | event >= length(raw_data)-wt) = [];
event_len = length(event);

% zero padding len
p = round(padding_window*fs);

% get all data
events_data = zeros(event_len, 2*wt+1);
raws_data = zeros(event_len, 2*wt+1);
for i = 1:event_len
    events_data(i,:) = event_data(event(i)-wt : event(i)+wt);
    raws_data(i,:) = raw_data(event(i)-wt : event(i)+wt);
end

% dummy run
[~, time, freq] = analysis.tfr(raws_data(1,:), raw_freq, fs);
time = time - wt/fs;

% find times
t1 = find(time<=-window, 1, 'last');
t2 = find(time>=window, 1, 'first');
t3 = find(time<=pre_event_start, 1, 'last');
t4 = find(time>=pre_event_stop, 1, 'first');

time = time(t1:t2);

% define vars for parfor
si = zeros(size(event));
meanRawData = zeros(1, t2-t1+1);
meanData = zeros(1, t2-t1+1);
meanTFR = zeros(length(freq), t2-t1+1);

parfor i = 1:event_len
    ed = events_data(i,:);
    rd = raws_data(i,:);

    [power, ~, ~] = analysis.tfr(rd, raw_freq, fs);

    if baseline
        pm = mean(power(:, t3:t4), 2);
        bi = (power - pm) ./ pm;
    else
        bi = power;
    end

    rp = mean(bi, 1);
    rp(isnan(rp)) = 0;

    rf = filter.fir(padarray(rp, [0 p], 0), fs, event_freq(1), event_freq(2));
    ef = filter.fir(padarray(ed, [0 p], 0), fs, event_freq(1), event_freq(2));

    rf = rf(p+t1 : p+t2);
    ef = ef(p+t1 : p+t2);

    ph = angle(hilbert(ef)) - angle(hilbert(rf));

    si(i) = mean(exp(1i * ph));
    meanTFR = meanTFR + bi(:, t1:t2) / event_len;
    meanRawData = meanRawData + ed(t1:t2) / event_len;
    meanData = meanData + ef / event_len;
end

end

