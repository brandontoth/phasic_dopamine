function corr_plots(data_struct1, varargin)

color = ['k', 'b', 'r', 'g', 'c', 'y', 'm'];

if isempty(varargin)
    %% REM
    figure
    tiledlayout(2, 2)
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 .3 .5])

    nexttile
    scatter(data_struct1.rem_len, data_struct1.rem_num, 'filled', color(1), 'SizeData', 5)
    hold on
    p = polyfit(data_struct1.rem_len, data_struct1.rem_num, 1);
    f = polyval(p, data_struct1.rem_len);
    plot(data_struct1.rem_len, f, color(1))
    [r, val] = corrcoef(data_struct1.rem_len, data_struct1.rem_num);
    title(['r^2 = ', num2str(r(2, 1)), ' p = ', num2str(val(2, 1))], 'FontWeight', 'normal')
    ylabel('Total transients (#)')
    xlabel('REM bout length (s)')
    pbaspect([1 1 1])

    nexttile
    scatter(data_struct1.rem_len, data_struct1.rem_rate, 'filled', color(1), 'SizeData', 5)
    hold on
    p = polyfit(data_struct1.rem_len, data_struct1.rem_rate, 1);
    f = polyval(p, data_struct1.rem_len);
    plot(data_struct1.rem_len, f, color(1))
    [r, val] = corrcoef(data_struct1.rem_len, data_struct1.rem_rate);
    title(['r^2 = ', num2str(r(2, 1)), ' p = ', num2str(val(2, 1))], 'FontWeight', 'normal')
    ylabel('Transient rate (Hz)')
    xlabel('REM bout length (s)')
    pbaspect([1 1 1])

    %% NREM
    nexttile
    scatter(data_struct1.nrem_len, data_struct1.nrem_num, 'filled', color(1), 'SizeData', 5)
    hold on
    p = polyfit(data_struct1.nrem_len, data_struct1.nrem_num, 1);
    f = polyval(p, data_struct1.nrem_len);
    plot(data_struct1.nrem_len, f, color(1))
    [r, val] = corrcoef(data_struct1.nrem_len, data_struct1.nrem_num);
    title(['r^2 = ', num2str(r(2, 1)), ' p = ', num2str(val(2, 1))], 'FontWeight', 'normal')
    ylabel('Total transients (#)')
    xlabel('NREM bout length (s)')
    pbaspect([1 1 1])

    nexttile
    scatter(data_struct1.nrem_len, data_struct1.nrem_rate, 'filled', color(1), 'SizeData', 5)
    hold on
    p = polyfit(data_struct1.nrem_len, data_struct1.nrem_rate, 1);
    f = polyval(p, data_struct1.nrem_len);
    plot(data_struct1.nrem_len, f, color(1))
    [r, val] = corrcoef(data_struct1.nrem_len, data_struct1.nrem_rate);
    title(['r^2 = ', num2str(r(2, 1)), ' p = ', num2str(val(2, 1))], 'FontWeight', 'normal')
    ylabel('Transient rate (Hz)')
    xlabel('NREM bout length (s)')
    pbaspect([1 1 1])

    set(gcf, 'renderer', 'Painters', 'PaperOrientation', 'landscape')
    print('corr', '-dpdf')

else
    %% REM
    figure
    scatter(data_struct1.rem_len, data_struct1.rem_num, color(1))
    hold on
    p = polyfit(data_struct1.rem_len, data_struct1.rem_num, 1);
    f = polyval(p, data_struct1.rem_len);
    plot(data_struct1.rem_len, f, color(1))
    
    for i = 1:size(varargin, 2)
        scatter(varargin{i}.rem_len, varargin{i}.rem_num, color(i + 1))
        p = polyfit(varargin{i}.rem_len, varargin{i}.rem_num, 1);
        f = polyval(p, varargin{i}.rem_len);
        plot(varargin{i}.rem_len, f, color(i + 1))
    end
    ylabel('Total transients (#)')
    xlabel('REM bout length (s)')

    figure
    scatter(data_struct1.rem_len, data_struct1.rem_rate, color(1))
    hold on
    p = polyfit(data_struct1.rem_len, data_struct1.rem_rate, 1);
    f = polyval(p, data_struct1.rem_len);
    plot(data_struct1.rem_len, f, color(1))
    
    for i = 1:size(varargin, 2)
        scatter(varargin{i}.rem_len, varargin{i}.rem_rate, color(i + 1))
        p = polyfit(varargin{i}.rem_len, varargin{i}.rem_rate, 1);
        f = polyval(p, varargin{i}.rem_len);
        plot(varargin{i}.rem_len, f, color(i + 1))
    end
    ylabel('Transient rate (Hz)')
    xlabel('REM bout length (s)')

    %% NREM
    figure
    scatter(data_struct1.nrem_len, data_struct1.nrem_num, color(1))
    hold on
    p = polyfit(data_struct1.nrem_len, data_struct1.nrem_num, 1);
    f = polyval(p, data_struct1.nrem_len);
    plot(data_struct1.nrem_len, f, color(1))
    
    for i = 1:size(varargin, 2)
        scatter(varargin{i}.nrem_len, varargin{i}.nrem_num, color(i + 1))
        p = polyfit(varargin{i}.nrem_len, varargin{i}.nrem_num, 1);
        f = polyval(p, varargin{i}.nrem_len);
        plot(varargin{i}.nrem_len, f, color(i + 1))
    end
    ylabel('Total transients (#)')
    xlabel('NREM bout length (s)')

    figure
    scatter(data_struct1.nrem_len, data_struct1.nrem_rate, color(1))
    hold on
    p = polyfit(data_struct1.nrem_len, data_struct1.nrem_rate, 1);
    f = polyval(p, data_struct1.nrem_len);
    plot(data_struct1.nrem_len, f, color(1))
    
    for i = 1:size(varargin, 2)
        scatter(varargin{i}.nrem_len, varargin{i}.nrem_rate, color(i + 1))
        p = polyfit(varargin{i}.nrem_len, varargin{i}.nrem_rate, 1);
        f = polyval(p, varargin{i}.nrem_len);
        plot(varargin{i}.nrem_len, f, color(i + 1))
    end
    ylabel('Transient rate (Hz)')
    xlabel('NREM bout length (s)')
end

end