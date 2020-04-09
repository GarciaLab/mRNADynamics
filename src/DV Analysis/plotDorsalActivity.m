function [fit, model] = plotDorsalActivity(x, y,activity, nc, DataType, ymean, se, plotScatter, varargin)
   
    xx = repmat(x, size(y,2), 1)';

    fitOpts = varargin;
    
    %p(1)=rate coefficient, p(2)=kd, p(3)=hill coefficient p(4) y offset
    [fit, model] = fitDorsalActivity(xx, y, fitOpts{:});
    idx = ~any(isnan(ymean),2);
    x4 = xx(idx);
%     xxx = min(x(:)) : 1:max(x4(end)*1.3,fit(2)*2.5);
    xxx = min(x(:)) : 1: x4(end)*1.3;
    
    figure('Units', 'points', 'Position', [0, 0, 200, 200]);
    clr = 'r';
    plot(xxx, model(fit,xxx), '-', 'DisplayName',['fit: ',num2str(round(fit))], 'LineWidth', 2, 'Color', clr);
    set(gca,'Units', 'normalized', 'OuterPosition', [0, 0, 1, 1]);
    axis square
    hold on
    if plotScatter
     errorbar(xx(idx), ymean(idx), se(idx), 'o', 'DisplayName', DataType, 'MarkerSize', 4, 'MarkerFaceColor', clr, 'CapSize', 0);
    end
    xlabel('dorsal concentration (au)');
    ylabel(activity);
    title([DataType, ' nc',num2str(nc+11)]);
    %p(1)=rate coefficient, p(2)=kd, p(3)=hill coefficient p(4) y offset
%     legend(DataType, {['fit: ',num2str(round(fit))], 'amplitude, KD(au), n, y offset'})
    legend(DataType, {[num2str(round(fit))], 'amplitude, KD(au), n, y offset'})

    leg = get(gca, 'Legend'); w=.02;h=.01;set(leg, 'Units', 'normalized', 'Position', [0.5,0.8,w, h], 'Box','off');

end