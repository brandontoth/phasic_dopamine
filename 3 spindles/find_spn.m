function [spn_loc, spn_rate] = find_spn(eeg, labels)
%% constants
fs      = 1000; % define sampling rate in Hz
d_fs    = 200;  % downsampled rate in Hz
win     = 150;  % define the size of the moving window in ms
spn_th  = 2.8;  % value used to determine height of threshold
min_len = 0.5;  % spindles should be at least .5 s in duration
max_len = 3;    % spindles should not exceed 3 s in duration

%% design filter
spn_filt = designfilt('bandpassiir',  ... 
    'StopbandFrequency1', 3, 'PassbandFrequency1', 10, ...     
    'PassbandFrequency2', 15, 'StopbandFrequency2', 22, ...     
    'StopbandAttenuation1', 24, 'StopbandAttenuation2', 24, ...
    'SampleRate', d_fs);

%% process eeg
% downsample signal
d_eeg = resample(eeg, d_fs, fs);

% get filtered eeg
eeg_filt = filtfilt(spn_filt, d_eeg); 

% calculate RMS with a moving window to create envelope
eeg_env = envelope(eeg_filt, win, 'rms'); 

% raise to the third power to amplify differences at higher thresholds
eeg_env = eeg_env .^ 3; 

%% identify spindles
% define threshold for spindle detection
thresh = spn_th * mean(eeg_env);

% find candidate spindles by identifying peaks above the threshold
st_spn = (find(diff(eeg_env >= thresh) ==  1));
en_spn = (find(diff(eeg_env >= thresh) == -1));

% catch if recording started during a spindle
if st_spn(1) > en_spn(1)    
    en_spn = en_spn(2:end);
end

% catch if recording ends during a spindle
if numel(st_spn) ~= numel(en_spn) 
    st_spn = st_spn(1:end - 1); 
end

% candidate spindle locations
can_spn = [st_spn en_spn];

%% remove spindles that are out of range
len     = (can_spn(:, 2) - can_spn(:, 1)) / d_fs;
can_spn = can_spn(len >= min_len & len <= max_len, :);

% spindle locations
spn_loc = can_spn(:, 1);

%% get spindle rates for all nrem bouts
slp  = parse_states(labels, 5, d_fs);
nrem = slp.nrem_loc * d_fs;

spn_rate = nan(size(nrem, 1), 1);
for i = 1:size(nrem, 1)
    idx = find(spn_loc >= nrem(i, 1) & spn_loc <= nrem(i, 2));
    dur = (nrem(i, 2) - nrem(i, 1)) / d_fs;

    spn_rate(i) = numel(idx) / dur;
end

end