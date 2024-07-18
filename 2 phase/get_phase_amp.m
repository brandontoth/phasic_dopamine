function amp = get_phase_amp(fp, ph)

% 10 degree phase bins
bins = -180:10:180;

% get amplitude in each bin, per bout
amp = nan(size(fp, 1), length(bins) - 1);
for i = 1:size(fp, 1)
    cur_fp = fp(i, :);
    cur_ph = ph(i, :);

    for j = 1:length(bins) - 1
        rng    = bins(j):bins(j + 1);
        
        amp(i, j) = mean(cur_fp((rng(1) <= cur_ph) & (cur_ph <= rng(end))));
    end
end

end