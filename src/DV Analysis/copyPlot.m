function copyPlot(ax1, ax2, varargin)

    closeFig = false;
    clrmpFlag = false;
    for i = 1:length(varargin)
        if strcmpi(varargin{i}, 'close')
            closeFig = true;
        elseif strcmpi(varargin{i}, 'colorMap')
            clrmpFlag = true;
            clr = varargin{i+1};
        end
    end
    
    if clrmpFlag 
        set(ax1.Children, 'Color', clr, 'MarkerFaceColor', clr);
    end
    for i = 1:length(ax1.Children)
        copyobj(ax1.Children(i), ax2)
    end

    figure(ax2.Parent);
    if closeFig
        close(ax1.Parent);
    end
    
end
