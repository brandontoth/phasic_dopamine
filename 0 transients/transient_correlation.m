function transient_correlation(photo_signal, labels, save_loc)

%% constant declaration
Fs = 1000;

%% get sleep state info
sl_str   = parse_states(labels, 5, Fs);
rem_len  = sl_str.rem_loc (:, 2) - sl_str.rem_loc (:, 1);
nrem_len = sl_str.nrem_loc(:, 2) - sl_str.nrem_loc(:, 1);

if ~isempty(sl_str.cat_loc); cat_len = sl_str.cat_loc(:, 2) - sl_str.cat_loc(:, 1); end

%% get transients
peak   = new_transient_detection(photo_signal, Fs);
tr     = quantify_transients(labels, peak, Fs);

%% put all data in a struct
corr.rem_len   = rem_len;
corr.nrem_len  = nrem_len;
corr.rem_rate  = tr.RemTransientTotRate;
corr.nrem_rate = tr.NremTransientTotRate;
corr.rem_num   = tr.RemCountPerBout;
corr.nrem_num  = tr.NremCountPerBout;

if exist('cat_len', 'var')
    corr.cat_len  = cat_len;
    corr.cat_rate = tr.CatTransientTotRate;
    corr.cat_num  = tr.CatCountPerBout;
end

%% save data
[~, name, ~] = fileparts(pwd);
save([save_loc '\' name '_tr_corr'], 'corr')

end