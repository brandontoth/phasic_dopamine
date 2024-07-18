function beh_state = get_beh_state_fp(signal, labels, food, lick, rw)
%% initialize struct, constant declaration
beh_state = struct;
epoch_len = 5;
Fs        = 1000;
down_fs   = 60;
len       = 30;

%% get fp phase
fp_ph  = get_fp_phase(signal', down_fs, Fs);
signal = resample(signal, down_fs, Fs);

%% get fp within sleep states
% find states locs
slp_str  = parse_states(labels, epoch_len, down_fs);

% upsample to fp
wake_loc = slp_str.wake_loc * down_fs;
nrem_loc = slp_str.nrem_loc * down_fs;
rem_loc  = slp_str.rem_loc  * down_fs;

% fp at wake
wake_len = wake_loc(:,2) - wake_loc(:,1);
wake_fp  = nan(size(wake_loc, 1), max(wake_len));
wake_ph  = nan(size(wake_loc, 1), max(wake_len));
for i = 1:size(wake_loc, 1)
    cur_ph = fp_ph (wake_loc(i, 1):wake_loc(i, 2));
    cur_fp = signal(wake_loc(i, 1):wake_loc(i, 2));
    
    wake_ph(i, 1:size(cur_ph, 2)) = cur_ph;
    wake_fp(i, 1:size(cur_fp, 2)) = cur_fp;
end
wake_ph = wake_ph(:, 1:end - 1);
wake_fp = wake_fp(:, 1:end - 1);

% fp at nrem
nrem_len = nrem_loc(:,2) - nrem_loc(:,1);
nrem_fp  = nan(size(nrem_loc, 1), max(nrem_len));
nrem_ph  = nan(size(nrem_loc, 1), max(nrem_len));
for i = 1:size(nrem_loc, 1)
    cur_ph = fp_ph (nrem_loc(i, 1):nrem_loc(i, 2));
    cur_fp = signal(nrem_loc(i, 1):nrem_loc(i, 2));
    
    nrem_ph(i, 1:size(cur_ph, 2)) = cur_ph;
    nrem_fp(i, 1:size(cur_fp, 2)) = cur_fp;
end
nrem_ph = nrem_ph(:, 1:end - 1);
nrem_fp = nrem_fp(:, 1:end - 1);

% fp at rem
rem_len = rem_loc(:,2) - rem_loc(:,1);
rem_fp  = nan(size(rem_loc, 1), max(rem_len));
rem_ph  = nan(size(rem_loc, 1), max(rem_len));
for i = 1:size(rem_loc, 1)
    cur_ph = fp_ph (rem_loc(i, 1):rem_loc(i, 2));
    cur_fp = signal(rem_loc(i, 1):rem_loc(i, 2));
    
    rem_ph(i, 1:size(cur_ph, 2)) = cur_ph;
    rem_fp(i, 1:size(cur_fp, 2)) = cur_fp;
end
rem_ph = rem_ph(:, 1:end - 1);
rem_fp = rem_fp(:, 1:end - 1);

% fp at cat (if applicable)
if ~isempty(slp_str.cat_loc)
    cat_loc = slp_str.cat_loc * down_fs;

    cat_len = cat_loc(:,2) - cat_loc(:,1);
    cat_fp  = nan(size(cat_loc, 1), max(cat_len));
    cat_ph  = nan(size(cat_loc, 1), max(cat_len));
    for i = 1:size(cat_loc, 1)
        cur_ph = fp_ph (cat_loc(i, 1):cat_loc(i, 2));
        cur_fp = signal(cat_loc(i, 1):cat_loc(i, 2));
        
        cat_ph(i, 1:size(cur_ph, 2)) = cur_ph;
        cat_fp(i, 1:size(cur_fp, 2)) = cur_fp;
    end
    cat_ph = cat_ph(:, 1:end - 1);
    cat_fp = cat_fp(:, 1:end - 1);
end

%% get fp for diff behaviors
% find behavior onset/offset times
beh = get_beh_times(food, lick, rw, down_fs, len);

% fp at feeding
food_len = beh.pellet_times(:,2) - beh.pellet_times(:,1);
food_fp  = nan(size(beh.pellet_times, 2), max(food_len) + 1);
food_ph  = nan(size(beh.pellet_times, 2), max(food_len) + 1);
if length(food_len) > 2
    for i = 1:length(beh.pellet_times)
        food_ph(i, :) = fp_ph (beh.pellet_times(i, 1):beh.pellet_times(i, 2));
        food_fp(i, :) = signal(beh.pellet_times(i, 1):beh.pellet_times(i, 2));
    end
end

% fp at licking
lick_len = beh.lick_times(:,2) - beh.lick_times(:,1);
lick_fp  = nan(size(beh.lick_times, 1), max(lick_len));
lick_ph  = nan(size(beh.lick_times, 1), max(lick_len));
if length(lick_len) > 2
    for i = 1:length(beh.lick_times)
        lick_ph(i, 1:lick_len(i) + 1) = fp_ph (beh.lick_times(i, 1):beh.lick_times(i, 2));
        lick_fp(i, 1:lick_len(i) + 1) = signal(beh.lick_times(i, 1):beh.lick_times(i, 2));
    end
end
lick_ph = lick_ph(:, 1:end - 1);
lick_fp = lick_fp(:, 1:end - 1);

% fp at rw
rw_len = beh.rw_times(:,2) - beh.rw_times(:,1);
rw_fp  = nan(size(beh.rw_times, 1), max(rw_len));
rw_ph  = nan(size(beh.rw_times, 1), max(rw_len));
if length(rw_len) > 2
    for i = 1:length(beh.rw_times)
        rw_ph(i, 1:rw_len(i) + 1) = fp_ph (beh.rw_times(i, 1):beh.rw_times(i, 2));
        rw_fp(i, 1:rw_len(i) + 1) = signal(beh.rw_times(i, 1):beh.rw_times(i, 2));
    end
end
rw_ph = rw_ph(:, 1:end - 1);
rw_fp = rw_fp(:, 1:end - 1);

%% save everything
beh_state.wake_ph = wake_ph;
beh_state.wake_fp = wake_fp;
beh_state.nrem_ph = nrem_ph;
beh_state.nrem_fp = nrem_fp;
beh_state.rem_ph  = rem_ph;
beh_state.rem_fp  = rem_fp;
beh_state.food_ph = food_ph;
beh_state.food_fp = food_fp;
beh_state.lick_ph = lick_ph;
beh_state.lick_fp = lick_fp;
beh_state.rw_ph   = rw_ph;
beh_state.rw_fp   = rw_fp;

if exist("cat_fp", "var"); beh_state.cat_fp = cat_fp;
    beh_state.cat_ph = cat_ph; end

end