function ph_rem = phasic_rem
%% load data in dir
files = dir('*.mat');
for i = 1:length(files)
    if contains(files(i).name, 'Signal') || contains(lower(files(i).name), 'labels', 'IgnoreCase', true) ...
        || contains(files(i).name, 'EEG', 'IgnoreCase', true)
        load(files(i).name)
    end
end

%% downsample eeg
% fs     = 1000;
d_fs   = 1000;
% d_eeg  = resample(EEG, d_fs, fs);
d_eeg = EEG;

%% define 11-element box filter
n_filt = 11;
filt   = ones(n_filt);
filt   = filt / sum(filt);

%% get rem locs and preallocate arrays
slp          = parse_states(labels, 5, d_fs);
rem_eeg      = [];
sm_diff_list = [];
tr_idx_list  = [];

rem_len     = (slp.rem_loc(:, 2) - slp.rem_loc(:, 1)) * d_fs;
max_len     = max(rem_len);
eeg_seq     = nan(size(slp.rem_loc, 1), max_len);
sm_diff_seq = nan(size(slp.rem_loc, 1), max_len);
tr_idx_seq  = nan(size(slp.rem_loc, 1), max_len);

%% iterate through rem bouts
for i = 1:size(slp.rem_loc, 1)
    % get eeg for the current bout
    cur_rng  = [slp.rem_loc(i, 1) slp.rem_loc(i, 2)] * d_fs;
    cur_rem  = d_eeg(cur_rng(1):cur_rng(2));
    
    % filter and get eeg phase and amplitude
    eeg_filt = theta_filt(cur_rem, d_fs);
    eeg_env  = hilbert(eeg_filt);
    eeg_ph   = angle(eeg_env);
    eeg_amp  = abs(eeg_env);
    
    % get the indices for troughs
    x        = (0:length(eeg_ph) - 1) + cur_rng(1);
    tr_idx   = islocalmin(eeg_ph, 'MinProminence', 3);
    tr_idx   = x(tr_idx);
    
    % get the diff between troughs and smooth
    tr_diff  = diff(tr_idx);
    sm_diff  = conv(tr_diff, filt, 'same');

    % append everything into lists
    rem_eeg      = [rem_eeg      eeg_amp];
    sm_diff_list = [sm_diff_list sm_diff];
    tr_idx_list  = [tr_idx_list  tr_idx];

    % store each bout in array
    eeg_seq    (i, 1:length(eeg_amp)) = eeg_amp;
    sm_diff_seq(i, 1:length(sm_diff)) = sm_diff;
    tr_idx_seq (i, 1:length(tr_idx))  = tr_idx;
end
    
%% define thresholds
thr1 = prctile(sm_diff_list, 15); % 15th percentile of smoothed diff
thr2 = prctile(sm_diff_list, 10); % 10th percentile of smoothed diff
thr3 = mean   (rem_eeg);          % mean theta of all rem bouts

%% iterate through candidate phasic rem events
ph_rem = [];
for i = 1:size(eeg_seq, 1)
    % get sequences for the current bout
    tr_idx  = rmmissing(tr_idx_seq (i, :));
    sm_diff = rmmissing(sm_diff_seq(i, :));
    eeg_amp = rmmissing(eeg_seq    (i, :));
    
    % find candidate phasic events based on smoothed diff
    idx  = find(sm_diff <= thr1);
    st   = find(diff(idx) > 1) + 1;
    en   = [st(2:end) - 1 length(idx)];
    cand = [idx(st)' idx(en)'];

    % iterate through candidate events
    for j = 1:size(cand, 1)
        rng = tr_idx - slp.rem_loc(i, 1) * d_fs;
        dur = ((tr_idx(cand(j, 2)) - tr_idx(cand(j, 1))) / d_fs) * 1000;
        
        % phasic event if at least 900 ms, if the shortest smoothed diff is
        % less than the 10th percentile of all smoothed diff, and if the
        % bout theta power is greater than theta power for all bouts 
        if dur > 900 && min(sm_diff(cand(j, :))) < thr2 && mean(eeg_amp(rng(cand(j, 1)):rng(cand(j, 2)))) > thr3
            ph_rem = [ph_rem; tr_idx(cand(j, 1)) tr_idx(cand(j, 2))];
        end
    end
end

end

%% nested functions
    function filt_eeg = theta_filt(d_eeg, d_fs)
        order      = 4;  % IIR filter order
        fcutlow    = 5;  % lower bandpass freq
        fcuthigh   = 12; % upper bandpass freq
        [b, a]     = butter(order, [fcutlow fcuthigh] / (d_fs / 2), 'bandpass');
        filt_eeg   = filtfilt(b, a, d_eeg');
    end