function amp = state_phase_amp(beh_state)
%% define struct
amp  = struct;

%% phase amplitude across state
% wake
wake_amp = get_phase_amp(beh_state.wake_fp, beh_state.wake_ph);

% nrem
nrem_amp = get_phase_amp(beh_state.nrem_fp, beh_state.nrem_ph);

% rem
rem_amp  = get_phase_amp(beh_state.rem_fp, beh_state.rem_ph);

% cataplexy (if applicable)
if isfield(beh_state, 'cat_fp')
    cat_amp = get_phase_amp(beh_state.cat_fp, beh_state.cat_ph);
end

%% phase amplitude across behavior
% food
food_amp = get_phase_amp(beh_state.food_fp, beh_state.food_ph);

% lick
lick_amp = get_phase_amp(beh_state.lick_fp, beh_state.lick_ph);

% rw
rw_amp   = get_phase_amp(beh_state.rem_fp, beh_state.rem_ph);

%% save everything
amp.wake = wake_amp;
amp.nrem = nrem_amp;
amp.rem  = rem_amp;
amp.food = food_amp;
amp.lick = lick_amp;
amp.rw   = rw_amp;

if exist('cat_amp', 'var'); amp.cat = cat_amp; end

end