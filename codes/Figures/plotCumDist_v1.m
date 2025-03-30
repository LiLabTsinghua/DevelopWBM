function [y_param, stats_param] = plotCumDist_v1(E_parameter,legends,titlestr,colors, median)
for i=1:length(E_parameter)
    
    if nargin<5
        median = false;
    end
      
    if median
        str{i} = horzcat(legends{i},' (',num2str(length(E_parameter{i})),' / ',num2str(median(E_parameter{i})), ')');
    else
        str{i} = legends{i};
    end
    [y_param(i),stats_param(i)] = cdfplot(cell2mat(E_parameter{i}));
    set(y_param(i), 'LineWidth', 2, 'Color', colors{i});
    
    title(titlestr,'FontName','Arial','FontSize',8,'FontWeight','bold')
    ylabel('Cumulative distribution','FontName','Arial','FontSize',8,'FontWeight','normal');
    xlabel('Variability range [mmol/gDw/h]','FontName','Arial','FontSize',8,'FontWeight','normal');
    set(gca,'FontSize',8, 'XScale', 'log','XGrid','off','YGrid','off', 'LineWidth',1, 'TickLength',[0.01 0.01],...
        'XTick',[10e-11, 10e-9, 10e-7, 10e-5, 10e-3, 10e-1, 10e1, 10e3], 'XMinorTick', 'off', 'YTick',[0, 0.2, 0.4, 0.6, 0.8, 1])
    
    hold on
    
end
legend(y_param,str,'FontName','Arial','FontSize',7);
end