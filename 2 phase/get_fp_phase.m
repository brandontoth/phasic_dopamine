function fp_ph = get_fp_phase(signal, downFs, Fs)
%% downsample fp
signal = resample(signal, downFs, Fs);

%% filter between .5 - 4 Hz
order    = 4;   % IIR filter order
fcutlow  = .5;  % lower bandpass freq
fcuthigh = 4;   % upper bandpass freq
[b, a]   = butter(order, [fcutlow fcuthigh] / (downFs / 2), 'bandpass');
filt_sig = filtfilt(b, a, signal');

%% extract phase
fp_ph = rad2deg(angle(hilbert(filt_sig)));

end