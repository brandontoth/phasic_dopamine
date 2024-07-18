function prior_state_to_nrem_tr
%% constants
epoch_len = 5;
Fs        = 1000;
min_len   = 25 * Fs;
rec_len   = 43200000;

%% load data
files = dir('*.mat');
for i = 1:length(files)
    if contains(files(i).name, 'EEG') || contains(lower(files(i).name), 'labels', 'IgnoreCase', true)
        load(files(i).name)
    end
end

%% get sleep states
slp  = parse_states (labels, epoch_len, Fs);
rem  = slp.rem_loc  * Fs; rem_dur  = rem (:, 2) - rem (:, 1);
nrem = slp.nrem_loc * Fs; nrem_dur = nrem(:, 2) - nrem(:, 1);
wake = slp.wake_loc * Fs; wake_dur = wake(:, 2) - wake(:, 1);

% remove short bouts
wake(wake_dur < min_len, :) = []; wake_dur(wake_dur < min_len) = [];

%% load transient data
[name, folder] = uigetfile('E:\2 transitions\0 NAcc\Analysis\1 transients - all bouts');
loc = fullfile(folder, name);
load(loc); %#ok<LOAD> 

%% get the NREM that follows REM sleep
prem_nrem_tr  = nan(size(rem, 1), 1);
prem_nrem_dur = nan(size(rem, 1), 1);
for i = 1:size(rem, 1)
    rem_end  = rem(i, 2);
    nrem_idx = find(nrem(:, 1) > rem_end);
    
    if isempty(nrem_idx); break; end
    nrem_idx = nrem_idx(1);

    prem_nrem_tr(i)  = rate.nrem_tot_rate(nrem_idx);
    prem_nrem_dur(i) = nrem_dur(nrem_idx);
end

post_rem_dat = [rem_dur/Fs prem_nrem_tr prem_nrem_dur/Fs];
post_rem_dat = rmmissing(post_rem_dat);

%% get the NREM that follows wake
pwake_nrem_tr  = nan(size(wake, 1), 1);
pwake_nrem_dur = nan(size(wake, 1), 1);
for i = 1:size(wake, 1)
    wake_end    = wake(i, 2);
    if wake_end == rec_len; break; end

    nrem_idx = find(nrem(:, 1) == wake_end);
    nrem_idx = nrem_idx(1);

    pwake_nrem_tr(i)  = rate.nrem_tot_rate(nrem_idx);
    pwake_nrem_dur(i) = nrem_dur(nrem_idx);
end

post_wake_dat = [wake_dur/Fs pwake_nrem_tr pwake_nrem_dur/Fs];
post_wake_dat = rmmissing(post_wake_dat);

%% get the NREM that precedes REM sleep
[~, idx, ~] = intersect(slp.nrem_loc(:, 2), slp.nrem_to_rem / Fs);
ntr_loc = nrem(idx, :);
ntr_dur = ntr_loc(:, 2) - ntr_loc(:, 1);
ntr_tr  = rate.nrem_tot_rate(idx);

pre_rem_dat = [rem_dur/Fs ntr_tr' ntr_dur/Fs];

%% get the NREM that precedes wake
[~, idx, ~] = intersect(slp.nrem_loc(:, 2), slp.nrem_to_wake / Fs);
ntw_loc = nrem(idx, :);
ntw_dur = ntw_loc(:, 2) - ntw_loc(:, 1);
ntw_tr  = rate.nrem_tot_rate(idx);

[~, idx, ~]  = intersect(wake(:, 1), slp.nrem_to_wake);
wake_dur     = wake_dur (idx);

pre_wake_dat = [wake_dur/Fs ntw_tr' ntw_dur/Fs];

%% save data
save_loc     = 'E:\2 transitions\0 NAcc\Analysis\10 prior REM len v NREM transients\0 wt light all bouts';
[~, name, ~] = fileparts(pwd);
save([save_loc '\' name '_post_REM_tr'],  'post_rem_dat')
save([save_loc '\' name '_post_Wake_tr'], 'post_wake_dat')
save([save_loc '\' name '_pre_REM_tr'],   'pre_rem_dat')
save([save_loc '\' name '_pre_Wake_tr'],  'pre_wake_dat')

clc
end