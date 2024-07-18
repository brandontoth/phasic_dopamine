function [w, pr] = quantify_width_prom(labels, widths, proms, Fs)
%% initialize structs
w  = struct;
pr = struct;

%% get sleep states
sl_strct = parse_states(labels, 5, Fs);

%% REM sleep stuff
rem_loc = sl_strct.rem_loc * Fs;
% rem_len = rem_loc(:, 2) - rem_loc(:, 1);
% rem_loc = rem_loc(20 * Fs <= rem_len, :);

if ~isempty(rem_loc)
    % get quartiles for widths and proms
    for i = 1:size(rem_loc, 1)
        % current bout
        rng  = rem_loc(i, 1):rem_loc(i, 2);
    
        % step size to divide bouts into quartiles
        step = floor(length(rng) / 4);
        
        % average width and prominence in that bout
        widths_r(i) = mean(rmmissing(widths(rng)));
        proms_r (i) = mean(rmmissing(proms (rng)));
        
        % step through REM bouts to get quartiles
        str = rng(1); fin = rng(1) + step;
        for k = 1:4
            cur_width_quart(k) = mean(rmmissing(widths(str:fin)));
            cur_proms_quart(k) = mean(rmmissing(proms (str:fin)));
            str = str + step; fin = fin + step;
        end
        
        % save to array
        rem_width_quart(i, :) = cur_width_quart;
        rem_proms_quart(i, :) = cur_proms_quart;
        clear cur_width_quart cur_proms_quart
    end
    
    % convert nan to zero
    rem_width_quart(isnan(rem_width_quart)) = 0;
    rem_proms_quart(isnan(rem_proms_quart)) = 0;
    
    % save to struct 
    w.rem_widths        = widths_r;
    w.rem_width_quart   = rem_width_quart;
    pr.rem_proms        = proms_r;
    pr.rem_proms_quart  = rem_proms_quart;
end

%% NREM sleep stuff
nrem_loc = sl_strct.nrem_loc * Fs;
% nrem_loc = [nrem_loc(:, 1) nrem_loc(:, 2) - (5 * Fs)]; % cut off last 5 s of each bout
% nrem_len = nrem_loc(:, 2) - nrem_loc(:, 1);
% nrem_loc = nrem_loc(20 * Fs <= nrem_len, :);

% get quartiles for widths and proms
for i = 1:size(nrem_loc, 1)
    % current bout
    rng  = nrem_loc(i, 1):nrem_loc(i, 2);

    % step size to divide bouts into quartiles
    step = floor(length(rng) / 4);
    
    % average width and prominence in that bout
    widths_n(i) = mean(rmmissing(widths(rng)));
    proms_n (i) = mean(rmmissing(proms (rng)));
    
    % step through NREM bouts to get quartiles
    str = rng(1); fin = rng(1) + step;
    for k = 1:4
        cur_width_quart(k) = mean(rmmissing(widths(str:fin)));
        cur_proms_quart(k) = mean(rmmissing(proms (str:fin)));
        str = str + step; fin = fin + step;
    end
    
    % save to array
    nrem_width_quart(i, :) = cur_width_quart;
    nrem_proms_quart(i, :) = cur_proms_quart;
    clear cur_width_quart cur_proms_quart
end

% convert nan to zero
nrem_width_quart(isnan(nrem_width_quart)) = 0;
nrem_proms_quart(isnan(nrem_proms_quart)) = 0;

% save to struct
w.nrem_widths       = widths_n;
w.nrem_width_quart  = nrem_width_quart;
pr.nrem_proms       = proms_n;
pr.nrem_proms_quart = nrem_proms_quart;

%% only NREM to REM
[~, ~, common] = intersect(rem_loc(:, 1), nrem_loc(:, 2), 'rows');
ntr_loc = nrem_loc;
ntr_loc = ntr_loc(common, :);

% get quartiles for widths and proms
if ~isempty(ntr_loc)
    for i = 1:size(ntr_loc, 1)
        % current bout
        rng  = ntr_loc(i, 1):ntr_loc(i, 2);
    
        % step size to divide bouts into quartiles
        step = floor(length(rng) / 4);
        
        % average width and prominence in that bout
        widths_ntr(i) = mean(rmmissing(widths(rng)));
        proms_ntr (i) = mean(rmmissing(proms (rng)));
        
        % step through NREM bouts to get quartiles
        str = rng(1); fin = rng(1) + step;
        for k = 1:4
            cur_width_quart(k) = mean(rmmissing(widths(str:fin)));
            cur_proms_quart(k) = mean(rmmissing(proms (str:fin)));
            str = str + step; fin = fin + step;
        end
        
        % save to array
        ntr_width_quart(i, :) = cur_width_quart;
        ntr_proms_quart(i, :) = cur_proms_quart;
        clear cur_width_quart cur_proms_quart
    end
    
    % convert nan to zero
    ntr_width_quart(isnan(ntr_width_quart)) = 0;
    ntr_proms_quart(isnan(ntr_proms_quart)) = 0;

    % save to struct
    w.ntr_widths        = widths_ntr;
    w.ntr_width_quart   = ntr_width_quart;
    pr.ntr_proms        = proms_ntr;
    pr.ntr_proms_quart  = ntr_proms_quart;
end

%% only NREM to wake
% get state locs
wake_loc = sl_strct.wake_loc * Fs;
[~, ~, brief] = intersect(sl_strct.brief_wake_loc, sl_strct.wake_loc, 'rows');
wake_loc(brief, :) = [];

[~, ~, common] = intersect(wake_loc(:, 1), nrem_loc(:, 2), 'rows');
ntw_loc = nrem_loc;
ntw_loc = ntw_loc(common, :);

if ~isempty(ntw_loc)
    % get quartiles for widths and proms
    for i = 1:size(ntw_loc, 1)
        % current bout
        rng  = ntw_loc(i, 1):ntw_loc(i, 2);
    
        % step size to divide bouts into quartiles
        step = floor(length(rng) / 4);
        
        % average width and prominence in that bout
        widths_ntw(i) = mean(rmmissing(widths(rng)));
        proms_ntw (i) = mean(rmmissing(proms (rng)));
        
        % step through NREM bouts to get quartiles
        str = rng(1); fin = rng(1) + step;
        for k = 1:4
            cur_width_quart(k) = mean(rmmissing(widths(str:fin)));
            cur_proms_quart(k) = mean(rmmissing(proms (str:fin)));
            str = str + step; fin = fin + step;
        end
        
        % save to array
        ntw_width_quart(i, :) = cur_width_quart;
        ntw_proms_quart(i, :) = cur_proms_quart;
        clear cur_width_quart cur_proms_quart
    end
    
    % convert nan to zero
    ntw_width_quart(isnan(ntw_width_quart)) = 0;
    ntw_proms_quart(isnan(ntw_proms_quart)) = 0;
    
    % save to struct
    w.ntw_widths        = widths_ntw;
    w.ntw_width_quart   = ntw_width_quart;
    pr.ntw_proms        = proms_ntw;
    pr.ntw_proms_quart  = ntw_proms_quart;
end

end