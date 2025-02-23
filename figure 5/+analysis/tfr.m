function [ps, time, freq] = tfr(data, base_frequency, fs)
cfg = [];
cfg.output       = 'pow';
cfg.channel      = 'ROI';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = (base_frequency(1):0.25:base_frequency(2));
cfg.t_ftimwin    = 5./cfg.foi;  % 5 cycles per time window
cfg.trials       = 1;
cfg.channel      = 1;
cfg.fsample      = fs;
cfg.toi          = (0:length(data)-1)/fs;
cfg.pad          = 'nextpow2';

eeg = pop_importdata('data', data, 'srate', fs);
eeg.nbchan = 1;
eeg.trials = 1;
fteeg = eeglab2fieldtrip(eeg, 'preprocessing');

fa = ft_freqanalysis(cfg, fteeg);

ps = squeeze(fa.powspctrm);
time = fa.time;
freq = fa.freq;
end


