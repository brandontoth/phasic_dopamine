function res = remove_drift(signal, Fs)
    fprintf('Processing data. \n');
    
    % Load data and add noise.
    y = resample(signal', 10, Fs);
            
    % Filter parameters
    fc = 0.0035;     % fc : cut-off frequency (cycles/sample)
    d  = 2;          % d  : filter order parameter (d = 1 or 2)
    
    % Positivity bias (peaks are positive)
    r = 5;          % r : asymmetry parameter
    
    % Regularization parameters
    amp  = 0.8;
    lam0 = 0.5 * amp;
    lam1 = 5 * amp;
    lam2 = 4 * amp;
    
    tic
    [~, f1, ~] = beads(y, d, fc, r, lam0, lam1, lam2);
    toc
    fprintf('\n');
    
    res = y - f1;
    res = resample(res', Fs, 10);
end