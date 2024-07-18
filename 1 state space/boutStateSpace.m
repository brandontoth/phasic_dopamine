function boutStateSpace
%% constant declaration
clear, clc, close all

epoch_len  = 5;    % epoch length in s
Fs         = 1000; % sampling rate in Hz
down_Fs    = 200;  % downsampled rate in Hz
wake_label = 2;    % AccuSleep wake label
nrem_label = 3;    % AccuSleep NREM label  
rem_label  = 1;    % AccuSleep REM label
cat_label  = 4;    % AccuSleep cataplexy label

%% load data
fprintf('Loading data. \n');

files = dir('*.mat');
for i = 1:length(files)
    if contains(files(i).name, 'EEG', 'IgnoreCase', true) || contains(files(i).name, 'Signal') || ...
            contains(lower(files(i).name), 'labels', 'IgnoreCase', true)
        load(files(i).name)
    end
end
 
if exist('zPhotoSyncRight', 'var')
    photo_signal = zPhotoSyncRight;
elseif exist('zPhotoSyncLeft', 'var')
    photo_signal = zPhotoSyncLeft;
elseif exist('photoSignal', 'var')
    photo_signal = photoSignal;
end

if exist('adjustedLabels', 'var')
    labels = adjustedLabels;
end

%% preprocess data and get transients
res        = remove_drift(photo_signal, Fs); % baseline correct signal
pk_array   = new_transient_detection(res, Fs);           % detect and make struct with transients
transients = quantify_transients(labels, pk_array, Fs);

% downsample EEG
signal = resample(EEG, down_Fs, Fs);

%% sort all sleep bouts
% get sleep states
slp_str = parse_states(labels, epoch_len, down_Fs);

% remove bouts < 20 s
wake_loc      = slp_str.wake_loc * down_Fs;
[~, ~, brief] = intersect(slp_str.brief_wake_loc, slp_str.wake_loc, 'rows');
wake_loc (brief, :) = [];

nrem_loc = slp_str.nrem_loc * down_Fs;
nrem_len = nrem_loc(:, 2) - nrem_loc(:, 1);
nrem_loc = nrem_loc(20 * down_Fs <= nrem_len, :);

rem_loc  = slp_str.rem_loc  * down_Fs;
rem_len  = rem_loc(:, 2) - rem_loc(:, 1);
rem_loc  = rem_loc(20 * down_Fs <= rem_len, :);

% concatenate bout LOCs, labels, and transient rates
wake_loc = [wake_loc zeros(length(wake_loc), 1) + wake_label transients.wake_tot_rate'];
nrem_loc = [nrem_loc zeros(length(nrem_loc), 1) + nrem_label transients.nrem_tot_rate'];
rem_loc  = [rem_loc  zeros(length(rem_loc),  1) + rem_label  transients.rem_tot_rate'];

if ~isempty(slp_str.cat_loc)
    cat_loc  = [slp_str.cat_loc * down_Fs zeros(length(slp_str.cat_loc), 1) + cat_label ...
    transients.cat_tot_rate'];
end

% sort data to iterate through
if ~exist('cat_loc', 'var')
    all_bouts = [wake_loc; nrem_loc; rem_loc];
else
    all_bouts = [wake_loc; nrem_loc; rem_loc; cat_loc];
end
all_bouts = sortrows(all_bouts);
bout_dur  = all_bouts(:, 2) - all_bouts(:, 1);

%% filter data from 1-60 Hz
order      = 4;  % IIR filter order
fcutlow    = 1;  % lower bandpass freq
fcuthigh   = 60; % upper bandpass freq
[b, a]     = butter(order, [fcutlow fcuthigh] / (down_Fs / 2), 'bandpass');
filt_sig   = filtfilt(b, a, signal');

%% run power analysis
ratio1 = nan(length(all_bouts), 1); % preallocate arrays
ratio2 = nan(length(all_bouts), 1);

for i    = 1:length(all_bouts)
    xdft = fft(filt_sig(all_bouts(i, 1):all_bouts(i, 2)));
    xdft = xdft(1:(bout_dur(i) / 2) + 1);
    psdx = (1 / (down_Fs * bout_dur(i))) * abs(xdft) .^ 2;
    psdx(2:end - 1) = 2 * psdx(2:end - 1);
    
    oneHz     = floor(length(psdx) / (down_Fs / 2));
    ratio1(i) = sum(psdx(oneHz * 6:oneHz * 10)) / sum(psdx(oneHz * 1:oneHz * 10));
    ratio2(i) = sum(psdx(oneHz * 1:oneHz * 16)) / sum(psdx(oneHz * 1:oneHz * 55));
end

%% plot everything
if ~exist('cat_loc', 'var')
    color = [0 1 0; 0 0 1; 1 0 0];
else
    color = [0 1 0; 0 0 1; 1 0 0; 1 0 1];
end

figure;
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 .6 .5])
ax1 = subplot(1, 2, 1);
scatter(ratio1, ratio2, 15, 'filled', 'CData', all_bouts(:, 3));
colormap(ax1, color);
shading interp  
colorbar();
xlabel('Ratio 1 (6-10/1-10 Hz)')
ylabel('Ratio 2 (1-16/1-55 Hz)')

ax2 = subplot(1, 2, 2);
scatter(ratio1, ratio2, 15, 'filled', 'CData', all_bouts(:, 4));
colormap(ax2, "jet"); caxis([0 0.35])
shading interp  
colorbar();
xlabel('Ratio 1 (6-10/1-10 Hz)')
ylabel('Ratio 2 (1-16/1-55 Hz)')

%% choose save location
save_loc  = uigetdir('', 'Choose save location');

%% save results
[~, name, ~] = fileparts(pwd);
save([save_loc, '\tr_state_space_all_bouts_light_', name, '.mat'], 'ratio1', 'ratio2', 'all_bouts')

clc
fprintf(['Done running analysis for ', name, '. \n']);
cd('..')

end