function fin = get_beh_times(food, lick, rw, down_fs, len)
%% constants
fin      = struct; % save data

%% downsample arrays
food = interp1(1:length(food), food, linspace(1, length(food), down_fs * 3600 * 12))';
lick = interp1(1:length(lick), lick, linspace(1, length(lick), down_fs * 3600 * 12))';
rw   = interp1(1:length(rw),   rw,   linspace(1, length(rw),   down_fs * 3600 * 12))';

% interpolation creates some edge effects
% zero out anything that isn't the max value to clean up TTLs
food(food ~= 1) = 0;
lick(lick ~= 1) = 0;
rw  (rw   ~= 5) = 0;

%% feeding
pellet_times  = find(diff(food) > 0.8);
pellet_onset  = pellet_times - 10 * down_fs;
pellet_offset = pellet_times + 10 * down_fs;

%% licking
lick_total  = find(diff(lick) == 1);
lick_onsets = lick_total;
for i = 2:length(lick_total)
    if lick_onsets(i) < lick_total(i - 1) + len * down_fs
        lick_onsets(i) = NaN;
    else
        lick_onsets(i) = lick_total(i);
    end
end
lick_onsets(isnan(lick_onsets)) = [];

lick_offsets = lick_total;
for i = 1:length(lick_total) - 1
    if lick_offsets(i) + len * down_fs < lick_total(i + 1)
        lick_offsets(i) = lick_total(i);       
    else
        lick_offsets(i) = NaN;
    end
end
lick_offsets(isnan(lick_offsets)) = [];

lick_len = lick_offsets - lick_onsets;
lick_onsets (lick_len < down_fs * len) = [];
lick_offsets(lick_len < down_fs * len) = [];

%% RW
rw_total  = find(diff(rw > 1));
rw_onsets = rw_total;

for i = 2:length(rw_total)
    if rw_onsets(i) < rw_total(i - 1) + len * down_fs
        rw_onsets(i) = NaN;
    else
        rw_onsets(i) = rw_total(i);
    end
end
rw_onsets(isnan(rw_onsets)) = [];

rw_offsets = rw_total;
for i = 1:length(rw_total) - 1
    if rw_offsets(i) + len * down_fs < rw_total(i + 1)
        rw_offsets(i) = rw_total(i);       
    else
        rw_offsets(i) = NaN;
    end
end
rw_offsets(isnan(rw_offsets)) = [];

rw_len = rw_offsets - rw_onsets;
rw_onsets (rw_len < down_fs * len) = [];
rw_offsets(rw_len < down_fs * len) = [];

%% save to struct
fin.pellet_times  = [pellet_onset pellet_offset];
fin.lick_times    = [lick_onsets lick_offsets];
fin.rw_times      = [rw_onsets rw_offsets];

end