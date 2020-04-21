function ax1 = plotInLoop(plotIndex, cmap, varargin)

xRange = [];
yRange = [];
ax1= [];
legendVisible = 'off';

%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end

if plotIndex == 1 && isempty(ax1)
    
    ax1 =gca;
    for k = 1:length(ax1.Children)
        if isprop(ax1.Children(k), 'Color')
            set(ax1.Children(k), 'Color', cmap(1, :));
        end
        if isprop(ax1.Children(k), 'MarkerFaceColor')
            set(ax1.Children(k), 'MarkerFaceColor', cmap(1, :));
        end
    end
    ax = gca;
    ax.Color = [253 249 207] /255;

    
else
    
    ax = gca;
    copyPlot(ax, ax1, 'close', 'colorMap', cmap(plotIndex, :));
    
    if isempty(xRange)
        xmax = get(gca, 'XLim');
        xlim([0, xmax(2)]);
    else
        xlim(xRange);
    end
    
    if isempty(yRange)
        ymax = get(gca, 'YLim');
        ylim([0, ymax(2)]);
    else
        ylim(yRange);
    end
    
    leg = get(gca, 'Legend'); w=.02;h=.01;set(leg, 'Units', 'normalized', 'Position', [1-w, 1-h,w, h], 'Box','off');
    set(leg, 'Visible', legendVisible);
    ax = gca;
     ax.Color = [253 249 207] /255;
    
end

end