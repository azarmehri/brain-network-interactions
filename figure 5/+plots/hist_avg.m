function hist_avg(s, r, range, bin, fs)

ind = (-range*fs : range*fs-1);
ind = reshape(ind, round(bin*fs), []);

start = s(s > range*fs & s < length(r)-range*fs);
event_rate = zeros(length(start), size(ind, 2));

for i = 1:length(start)
    rc = sum(r(ind + start(i)));
    event_rate(i, :) = rc;
end

p = mean(event_rate)/bin;
p = smoothdata(p, "gaussian", 5);
c = mean(ind)/fs;

bar(c, p)

end

