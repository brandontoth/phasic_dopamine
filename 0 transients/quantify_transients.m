function results = quantify_transients(labels, transients, Fs)
%% intialize struct
results = struct;

%% extend labels file to count transients per state
extend  = extendLabels(labels, transients, 5, Fs);

peaks_w = transients(extend == 2);
peaks_w(isnan(peaks_w)) = [];

peaks_n = transients(extend == 3);
peaks_n(isnan(peaks_n)) = [];

peaks_r = transients(extend == 1);
peaks_r(isnan(peaks_r)) = [];

%% get sleep states
sl_strct = parse_states(labels, 5, Fs);

%% wake stuff
% get state locs
wake_loc = sl_strct.wake_loc * Fs;

% remove brief arousals
[~, ~, brief] = intersect(sl_strct.brief_wake_loc, sl_strct.wake_loc, 'rows');
wake_loc(brief, :) = [];

% iterate through bouts
for i = 1:size(wake_loc, 1)
    % current bout
    rng = wake_loc(i, 1):wake_loc(i, 2);

    % step size to divide bouts into quartiles
    step = floor(length(rng) / 4);
    
    % num of transients in bout
    cur_bout = transients(rng);
    cur_bout(isnan(cur_bout)) = [];
    cur_num  = numel(cur_bout);
    
    % rate of transients in bout
    rate_w(i) = numel(cur_bout) / (numel(rng) / Fs);
    num_w (i) = cur_num;

    % step through wake bouts to get quartiles
    str = rng(1); fin = rng(1) + step;
    for k = 1:4
        q_num     = sum(~isnan(transients(str:fin)));
        q_rate(k) = q_num / (step / Fs); 
        
        q_amp(k)  = mean(rmmissing(transients(str:fin)));

        str = str + step; fin = fin + step;
    end

    % save to array
    wake_trans_quart(i ,:) = q_rate;
    wake_amp_quart  (i, :) = q_amp;
    clear q_rate q_amp
end

% remove nans
wake_trans_quart(isnan(wake_trans_quart)) = 0;
wake_amp_quart(isnan(wake_amp_quart)) = 0;

% add to struct
results.wake_dff        = peaks_w;
results.wake_num        = num_w;
results.wake_tot_rate   = rate_w;
results.wake_avg_rate   = mean(rate_w);
results.wake_avg_dff    = mean(peaks_w);
results.wake_skew       = skewness(peaks_w);
results.wake_kurt       = kurtosis(peaks_w);
results.wake_rate_quart = wake_trans_quart;
results.wake_amp_quart  = wake_amp_quart;

%% NREM stuff
if ~isempty(sl_strct.nrem_loc)
    % get states locs
    nrem_loc = sl_strct.nrem_loc * Fs;
%     nrem_loc = [nrem_loc(:, 1) nrem_loc(:, 2) - (5 * Fs)]; % cut off last 5 s of each bout
%     nrem_len = nrem_loc(:, 2) - nrem_loc(:, 1);
%     nrem_loc = nrem_loc(20 * Fs <= nrem_len, :);

    % iterate through bouts
    for i = 1:size(nrem_loc, 1)
        % current bout
        rng = nrem_loc(i, 1):nrem_loc(i, 2);

        % step size to divide bouts into quartiles
        step = floor(length(rng) / 4);
        
        % num of transients in bout
        cur_bout = transients(rng);
        cur_bout(isnan(cur_bout)) = [];
        cur_num  = numel(cur_bout);
        
        % rate of transients in bout
        rate_n(i) = numel(cur_bout) / (numel(rng) / Fs);
        num_n (i) = cur_num;
        
        % step through NREM bouts to get quartiles
        str = rng(1); fin = rng(1) + step;
        for k = 1:4
            q_num     = sum(~isnan(transients(str:fin)));
            q_rate(k) = q_num / (step / Fs); 

            q_amp(k)  = mean(rmmissing(transients(str:fin)));

            str = str + step; fin = fin + step;
        end

        % save to array
        nrem_trans_quart(i ,:) = q_rate;
        nrem_amp_quart  (i, :) = q_amp;
        clear q_rate q_amp
    end
    
    % remove nans
    nrem_trans_quart(isnan(nrem_trans_quart)) = 0;
    nrem_amp_quart(isnan(nrem_amp_quart)) = 0;
    
    % add to struct
    results.nrem_dff        = peaks_n;
    results.nrem_num        = num_n;
    results.nrem_tot_rate   = rate_n;
    results.nrem_avg_rate   = mean(rate_n);
    results.nrem_avg_dff    = mean(peaks_n);
    results.nrem_skew       = skewness(peaks_n);
    results.nrem_kurt       = kurtosis(peaks_n);
    results.nrem_rate_quart = nrem_trans_quart;
    results.nrem_amp_quart  = nrem_amp_quart;
end

%% REM stuff
% get state locs
rem_loc = sl_strct.rem_loc * Fs;
% rem_len = rem_loc(:, 2) - rem_loc(:, 1);
% rem_loc = rem_loc(20 * Fs <= rem_len, :);

if ~isempty(rem_loc)
    % iterate through bouts
    for i = 1:size(rem_loc, 1)
        % current bout
        rng = rem_loc(i, 1):rem_loc(i, 2);

        % step size to divide bouts into quartiles
        step = floor(length(rng) / 4);
        
        % num of transients in bout
        cur_bout = transients(rng);
        cur_bout(isnan(cur_bout)) = [];
        cur_num  = numel(cur_bout);
        
        % rate of transients in bout
        rate_r(i) = numel(cur_bout) / (numel(rng) / Fs);
        num_r (i) = cur_num;

        % step through REM bouts to get quartiles
        str = rng(1); fin = rng(1) + step;
        for k = 1:4
            q_num     = sum(~isnan(transients(str:fin)));
            q_rate(k) = q_num / (step / Fs); 

            q_amp(k)  = mean(rmmissing(transients(str:fin)));

            str = str + step; fin = fin + step;
        end

        % save to array
        rem_trans_quart(i ,:) = q_rate;
        rem_amp_quart  (i, :) = q_amp;
        clear q_rate q_amp
    end

    % remove nans
    rem_trans_quart(isnan(rem_trans_quart)) = 0;
    rem_amp_quart(isnan(rem_amp_quart)) = 0;

    % add to struct
    results.rem_dff        = peaks_r;
    results.rem_num        = num_r;
    results.rem_tot_rate   = rate_r;
    results.rem_avg_rate   = mean(rate_r);
    results.rem_avg_dff    = mean(peaks_r);
    results.rem_skew       = skewness(peaks_r);
    results.rem_kurt       = kurtosis(peaks_r);
    results.rem_rate_quart = rem_trans_quart;
    results.rem_amp_quart  = rem_amp_quart;
end

%% NREM to REM
if ~isempty(rem_loc)
    % get state locs
    [~, ~, common] = intersect(rem_loc(:, 1), nrem_loc(:, 2), 'rows');
    ntr_loc = nrem_loc;
    ntr_loc = ntr_loc(common, :);
else
    ntr_loc = [];
end

if ~isempty(ntr_loc)
    % iterate through bouts
    for i = 1:size(ntr_loc, 1)
        % current bout
        rng = ntr_loc(i, 1):ntr_loc(i, 2);

        % mean within bout
        ntr_rate(i) = sum(~isnan(transients(rng))) / (numel(rng) / Fs);

        % step size to divide bouts into quartiles
        step = floor(length(rng) / 4);
        
        % step through REM bouts to get quartiles
        str = rng(1); fin = rng(1) + step;
        for k = 1:4
            q_num     = sum(~isnan(transients(str:fin)));
            q_rate(k) = q_num / (step / Fs); 

            q_amp(k)  = mean(rmmissing(transients(str:fin)));

            str = str + step; fin = fin + step;
        end

        % save to array
        ntr_trans_quart(i, :) = q_rate;
        ntr_amp_quart  (i, :) = q_amp;
        clear q_rate q_amp
    end

    % remove nans
    ntr_trans_quart(isnan(ntr_trans_quart)) = 0;
    ntr_amp_quart(isnan(ntr_amp_quart)) = 0;

    % add to struct
    results.ntr_rate       = ntr_rate;
    results.ntr_rate_quart = ntr_trans_quart;
    results.ntr_amp_quart  = ntr_amp_quart;
end

%% NREM to wake
% get state locs
[~, ~, common] = intersect(wake_loc(:, 1), nrem_loc(:, 2), 'rows');
ntw_loc = nrem_loc;
ntw_loc = ntw_loc(common, :);

if ~isempty(ntw_loc)
    % iterate through bouts
    for i = 1:size(ntw_loc, 1)
        % current bout
        rng = ntw_loc(i, 1):ntw_loc(i, 2);

        % mean within bout
        ntw_rate(i) = sum(~isnan(transients(rng))) / (numel(rng) / Fs);

        % step size to divide bouts into quartiles
        step = floor(length(rng) / 4);
        
        % step through REM bouts to get quartiles
        str = rng(1); fin = rng(1) + step;
        for k = 1:4
            q_num     = sum(~isnan(transients(str:fin)));
            q_rate(k) = q_num / (step / Fs); 

            q_amp(k)  = mean(rmmissing(transients(str:fin)));

            str = str + step; fin = fin + step;
        end

        % save to array
        ntw_trans_quart(i, :) = q_rate;
        ntw_amp_quart  (i, :) = q_amp;
        clear q_rate q_amp
    end

    % remove nans
    ntw_trans_quart(isnan(ntw_trans_quart)) = 0;
    ntw_amp_quart(isnan(ntw_amp_quart)) = 0;

    % add to struct
    results.ntw_rate       = ntw_rate;
    results.ntw_rate_quart = ntw_trans_quart;
    results.ntw_amp_quart  = ntw_amp_quart;
end

%% CAT stuff
if ~isempty(sl_strct.cat_loc)
    cat_loc = sl_strct.cat_loc * Fs;
    for i = 1:size(cat_loc, 1)
        rng = cat_loc(i, 1):cat_loc(i, 2);
    
        cur_bout = transients(rng);
        cur_bout(isnan(cur_bout)) = [];
        cur_num  = numel(cur_bout);
    
        rate_c(i) = numel(cur_bout) / (numel(rng) / Fs);
        num_c (i) = cur_num;
    end

    peaks_c = transients(extend == 4);
    peaks_c(isnan(peaks_c)) = [];

    results.cat_dff      = peaks_c;
    results.cat_num      = num_c;
    results.cat_tot_rate = rate_c;
    results.cat_avg_rate = mean(rate_c);
    results.cat_avg_dff  = mean(peaks_c);
    results.cat_skew     = skewness(peaks_c);
    results.cat_kurt     = kurtosis(peaks_c);
end

end