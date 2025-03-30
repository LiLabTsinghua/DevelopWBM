function [y_param, stats_param] = plotCumDist_v2(E_parameter, legends, colors)
    figure;
    hold on;
    
    % set plot size
    set(gcf, 'Units', 'Inches', 'Position', [0, 0, 4, 3.5]);

    for i = 1:length(E_parameter)
        [y_param(i), stats_param(i)] = cdfplot(cell2mat(E_parameter{i}));
        set(y_param(i), 'LineWidth', 2, 'Color', colors{i});
    end

    % set FontName, FontSize and Ticks
    set(gca, 'FontName', 'Arial', 'FontSize', 8, ...
             'XScale', 'log', 'XGrid', 'off', 'YGrid', 'off', ...
             'LineWidth', 0.5, 'TickLength', [0.015, 0.015], ...
             'XTick', [10e-11, 10e-9, 10e-7, 10e-5, 10e-3, 10e-1, 10e1, 10e3], ...
             'XMinorTick', 'off', 'YTick', [0, 0.2, 0.4, 0.6, 0.8, 1]);

    title('');

    % labels
    ylabel('Cumulative distribution', 'FontName', 'Arial', 'FontSize', 8, 'FontWeight', 'normal', 'Color', 'k');
    xlabel('Variability range [mmol/gDw/h]', 'FontName', 'Arial', 'FontSize', 8, 'FontWeight', 'normal', 'Color', 'k');

    % set legend
    legend(y_param, legends, 'FontName', 'Arial', 'FontSize', 8, 'Box', 'off', 'TextColor', 'k');

    ax = gca;
    ax.Box = 'on';
    ax.TickLength = [0.015 0.015];
    set(ax.XAxis, 'Color', 'k');
    set(ax.YAxis, 'Color', 'k');

    hold off;
end