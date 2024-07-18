function ph = ph_main
%% clear workspace
clear, clc, close all
Fs = 1000;

%% choose save location
save_loc  = uigetdir('', 'Choose save location');

%% navigate to desired dir
selpath = uigetdir();
filedir = selpath;
cd(filedir);

D = dir;
D = D(~ismember({D.name}, {'.', '..'}));
D = D([D(:).isdir]);

%% iterate through folders
for k = 1:numel(D)
    currD = D(k).name; % get the current subdirectory name
    cd(currD)          % change the directory

    fprintf(1, 'Now reading %s\n', currD);

    %% load data preprocess data
    % load files in dir
    [photo_signal, raw, labels, food, lick, rw] = load_data;
    
    res = remove_drift(photo_signal, Fs); % get rid of REM drop

    %% get phase amplitudes
    beh_state = get_beh_state_fp (res, labels, food, lick, rw);
    amp       = state_phase_amp  (beh_state);
    state_fft = get_beh_state_fft(raw, labels, food, lick, rw);
    
    %% save everything
    [~, name, ~] = fileparts(pwd);
    save([save_loc '\' name '_beh_phase_fp'],  'beh_state')
    save([save_loc '\' name '_beh_phase_amp'], 'amp')
    save([save_loc '\' name '_beh_state_fft'], 'state_fft')

    clear photo_signal labels food lick rw beh_state amp res state_fft

    %% return to the parent directory
    cd('..')
end

%% concatenate across mice
ph = concat_phase(save_loc);

clc, fprintf('Done running analysis. \n');
end