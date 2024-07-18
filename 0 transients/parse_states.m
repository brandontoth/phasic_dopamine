function results = parse_states(labels, epoch_len, Fs)
% Programmer: BA Toth

% Date: June 2023

% Purpose: identify start and stop times for each bout of WAKE, NREM,
% REM, and cataplexy and return arrays with this information

%% constants
NREM_LABEL = 3; % labels for Accusleep purposes
REM_LABEL  = 1;
WAKE_LABEL = 2;
CAT_LABEL  = 4;
results    = struct; % structure to hold stuff

%% if there are cataplexy episodes, locate them
if sum(ismember(labels, 4)) > 0
    cat_on  = epoch_len * (find(labels(1:end - 1) ~= CAT_LABEL & labels(2:end) == CAT_LABEL));
    cat_off = epoch_len * (find(labels(1:end - 1) == CAT_LABEL & labels(2:end) ~= CAT_LABEL));

    if labels(end) == CAT_LABEL
        cat_on = cat_on(1:end - 1);
    end

    cat_loc = [cat_on cat_off];
else
    cat_on  = [];
    cat_off = [];
    cat_loc = [];
end

%% find REM locations
rem_on  = epoch_len * (find(labels(1:end - 1) ~= REM_LABEL & labels(2:end) == REM_LABEL));
rem_off = epoch_len * (find(labels(1:end - 1) == REM_LABEL & labels(2:end) ~= REM_LABEL));

if labels(1) == REM_LABEL
    rem_off = rem_off(2:end);
end

if labels(end) == REM_LABEL
    rem_on = rem_on(1:end - 1);
end

rem_loc = [rem_on rem_off];

%% find NREM locations
nrem_on  = epoch_len * (find(labels(1:end - 1) ~= NREM_LABEL & labels(2:end) == NREM_LABEL));
nrem_off = epoch_len * (find(labels(1:end - 1) == NREM_LABEL & labels(2:end) ~= NREM_LABEL));

if labels(1) == NREM_LABEL
    beginning = 1;
    nrem_on = [beginning; nrem_on];
end

if labels(end) == NREM_LABEL
    nrem_off = [nrem_off; epoch_len * length(labels)];
end
nrem_loc = [nrem_on nrem_off];

%% find wake locations
wake_on  = epoch_len * (find(labels(1:end - 1) ~= WAKE_LABEL & labels(2:end) == WAKE_LABEL));
wake_off = epoch_len * (find(labels(1:end - 1) == WAKE_LABEL & labels(2:end) ~= WAKE_LABEL));

if labels(1) == WAKE_LABEL
    beginning = 1;
    wake_on = [beginning; wake_on];
end

if labels(end) == WAKE_LABEL
    wake_off = [wake_off; epoch_len * length(labels)];
end

wake_loc = [wake_on wake_off];
wake_len = wake_off - wake_on;

% remove all wake bouts shorter than 20 s
long_wake_loc  = wake_loc(wake_len >  20, :);
brief_wake_loc = wake_loc(wake_len <= 20, :);

% long state transitions
long_wtn = intersect(long_wake_loc(:, 2), nrem_on);
long_ntw = intersect(long_wake_loc(:, 1), nrem_off);

%% save everything to the results structure for convenience
% brief arousals
results.num_brief_wake = length(brief_wake_loc);
results.brief_wake_loc = brief_wake_loc;

% state locations
results.cat_loc  = cat_loc;
results.rem_loc  = rem_loc;
results.nrem_loc = nrem_loc;
results.wake_loc = wake_loc;

% cataplexy info
if ~isempty(cat_loc)
    results.num_cat_bouts = size(cat_loc, 1);
    results.tot_cat_time  = sum((cat_loc(:, 2) - cat_loc(:, 1)));
    results.avg_cat_len   = results.tot_cat_time / results.num_cat_bouts;
else
    results.num_cat_bouts = 0;
    results.tot_cat_time  = 0;
    results.avg_cat_len   = 0;
end

% rem info
if ~isempty(rem_loc)
    results.num_rem_bouts = size(rem_loc, 1);
    results.tot_rem_time  = sum((rem_loc(:, 2) - rem_loc(:, 1)));
    results.avg_rem_len   = results.tot_rem_time / results.num_rem_bouts;
else
    results.num_rem_bouts = 0;
    results.tot_rem_time  = 0;
    results.avg_rem_len   = 0;
end

% nrem info
results.num_nrem_bouts = size(nrem_loc, 1);
results.tot_nrem_time  = sum((nrem_loc(:, 2) - nrem_loc(:, 1)));
results.avg_nrem_len   = results.tot_nrem_time / results.num_nrem_bouts;

% wake info
results.num_wake_bouts = size(long_wake_loc, 1);
results.tot_wake_time  = sum((long_wake_loc(:, 2) - long_wake_loc(:, 1))) +...
    sum(brief_wake_loc(:, 2) - brief_wake_loc(:, 1));
results.avg_wake_len   = results.tot_wake_time / results.num_wake_bouts;

% state transitions 
results.wake_to_cat  = cat_on   * Fs;
results.cat_to_wake  = cat_off  * Fs;
results.nrem_to_rem  = rem_on   * Fs;
results.rem_to_wake  = rem_off  * Fs;
results.wake_to_nrem = long_wtn * Fs;
results.nrem_to_wake = long_ntw * Fs;

end