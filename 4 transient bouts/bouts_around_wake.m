function bouts_around_wake
%% constants
epoch_len = 5;
Fs        = 1000;
min_len   = 20 * Fs;
cut_off   = 60 * Fs;
ln        = 5;

%% load data
files = dir('*.mat');
for i = 1:length(files)
    if contains(files(i).name, 'EEG') || contains(lower(files(i).name), 'labels', 'IgnoreCase', true)
        load(files(i).name)
    end
end

%% get sleep states
slp  = parse_states (labels, epoch_len, Fs);
wake = slp.wake_loc * Fs; wake_dur = wake(:, 2) - wake(:, 1);
wake = wake(wake_dur > min_len, :); wake_dur = wake_dur(wake_dur > min_len);
nrem = slp.nrem_loc * Fs;

wake_sh = wake(wake_dur <  cut_off, :);
wake_ln = wake(wake_dur >= cut_off, :);

%% load transient data
[name, folder] = uigetfile('E:\2 transitions\0 NAcc\Analysis\1 transients - all bouts');
loc = fullfile(folder, name);
load(loc); %#ok<LOAD> 

%% get the NREM that follows short wake
nrem_tr_sh_a = nan(size(wake_sh, 1), ln);
nrem_tr_sh_b = nan(size(wake_sh, 1), ln);

for i = 1:size(wake_sh, 1)
    wake_end    = wake_sh(i, 2);
    nrem_idx_a = find(nrem(:, 1) >= wake_end);
    nrem_idx_b = find(nrem(:, 1) <  wake_end);

    chk_a = length(nrem_idx_a);
    if chk_a > ln; nrem_idx_a = nrem_idx_a(1:ln); end

    chk_b = length(nrem_idx_b);
    if chk_b > ln; nrem_idx_b = nrem_idx_b(end - ln + 1:end); end

    nrem_tr_sh_a(i, 1:length(nrem_idx_a)) = rate.nrem_tot_rate(nrem_idx_a);
    nrem_tr_sh_b(i, end - length(nrem_idx_b) + 1:end) = rate.nrem_tot_rate(nrem_idx_b);
end

%% get the NREM that follows long wake
nrem_tr_ln_a = nan(size(wake_ln, 1), ln);
nrem_tr_ln_b = nan(size(wake_ln, 1), ln);

for i = 1:size(wake_ln, 1)
    wake_end    = wake_ln(i, 2);
    nrem_idx_a = find(nrem(:, 1) >= wake_end);
    nrem_idx_b = find(nrem(:, 1) <  wake_end);

    chk_a = length(nrem_idx_a);
    if chk_a > ln; nrem_idx_a = nrem_idx_a(1:ln); end

    chk_b = length(nrem_idx_b);
    if chk_b > ln; nrem_idx_b = nrem_idx_b(end - ln + 1:end); end

    nrem_tr_ln_a(i, 1:length(nrem_idx_a)) = rate.nrem_tot_rate(nrem_idx_a);
    nrem_tr_ln_b(i, end - length(nrem_idx_b) + 1:end) = rate.nrem_tot_rate(nrem_idx_b);
end

%% save data
save_loc     = 'E:\2 transitions\0 NAcc\Analysis\12 nrem around wake\0 wt light';
[~, name, ~] = fileparts(pwd);
save([save_loc '\' name '_nrem_tr_sh_a'], 'nrem_tr_sh_a')
save([save_loc '\' name '_nrem_tr_sh_b'], 'nrem_tr_sh_b')
save([save_loc '\' name '_nrem_tr_ln_a'], 'nrem_tr_ln_a')
save([save_loc '\' name '_nrem_tr_ln_b'], 'nrem_tr_ln_b')

clc
end