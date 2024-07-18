function new_analysis_main
%% housekeeping
% clear workspace
clear, clc, close all

% declare some variables
Fs  = 1000;   % sampling rate

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
    [photo_signal, ~, labels, ~, ~, ~] = load_data;
    
    res = remove_drift(photo_signal, Fs); % get rid of REM drop

    %% transient correlational stuff
%     transient_correlation(res, labels, save_loc);
    [locs, widths, proms] = new_transient_detection(res, Fs);
    [w, pr] = quantify_width_prom(labels, widths, proms, Fs);
    rate    = quantify_transients(labels, locs, Fs);
%     sub     = subsample(labels, res, 5, Fs);
%     slp     = parse_states(labels, 5, Fs);
%     dff     = plotStateDFF(res, labels, slp.brief_wake_loc, Fs, 5);
%     early_tr_rate = tr_rate_early_late(labels, locs, Fs, "early");
%     late_tr_rate  = tr_rate_early_late(labels, locs, Fs, "late");

    [~, name, ~] = fileparts(pwd);
    save([save_loc '\' name '_tr_width'], 'w')
    save([save_loc '\' name '_tr_proms'], 'pr')
    save([save_loc '\' name '_tr_rate'] , 'rate')
%     save([save_loc '\' name '_sub']     , 'sub')
%     save([save_loc '\' name '_dff']     , 'dff')
%     save([save_loc '\' name '_early_tr_rate'], 'early_tr_rate')
%     save([save_loc '\' name '_late_tr_rate'],  'late_tr_rate')

    %% return to the parent directory
    cd('..')
end

clc, fprintf('Done running analysis. \n');
end