function [loc_array, width_array, prom_array] = new_transient_detection(photo_signal, Fs)

[pks, locs, widths, proms] = findpeaks(photo_signal, "MinPeakHeight", mean(photo_signal) ...
    + std(photo_signal) * 2, "MinPeakDistance", Fs / 2, ...
    "MinPeakProminence", mean(photo_signal) + std(photo_signal) * 2);

% findpeaks(photo_signal, "MinPeakHeight", mean(photo_signal) ...
%     + std(photo_signal) * 2, "MinPeakDistance", Fs / 2, ...
%     "MinPeakProminence", mean(photo_signal) + std(photo_signal) * 2, ...
%     "Annotate", "extents")

loc_array = nan(1, length(photo_signal));
for i = 1:length(locs)
    cur = locs(i);
    loc_array(cur) = pks(i);
end

width_array = nan(1, length(photo_signal));
for i = 1:length(locs)
    cur = locs(i);
    width_array(cur) = widths(i);
end

prom_array = nan(1, length(photo_signal));
for i = 1:length(locs)
    cur = locs(i);
    prom_array(cur) = proms(i);
end

end