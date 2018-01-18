function standardizeFigure(ax, legend, varargin)
    
    color(1,:) = [0 0 0];
    axesLineWidth = 5;
    fig = gcf;
    dataObj = get(ax, 'Children');
    dataType = get(dataObj, 'Type');
    if ~iscell(dataType)
        dataType = {dataType};
    end
    legendSize = 20;
    
    for i = 1:length(varargin)
       if strcmpi(varargin{i}, 'axeslinewidth')
            axesLineWidth = varargin{i+1};            
        elseif strcmpi(varargin{i}, 'red')
            color(i,:) = [213,108,85]/255;
        elseif strcmpi(varargin{i}, 'yellow')
            color(i,:) = [234,194,100]/255;
        elseif strcmpi(varargin{i}, 'cyan')
            color(i,:) = [108,188,233]/255;
        elseif strcmpi(varargin{i}, 'magenta')
            color(i,:) = [208,109,171]/255;
        elseif strcmpi(varargin{i}, 'lightblue')
            color(i,:) = [115,142,193]/255;
       elseif strcmpi(varargin{i}, 'legendFontSize')
            legendSize = varargin{i+1};
        end
    end
    
    if ~isempty(legend)
        legend.FontSize = legendSize;
        legend.Box = 'off';
    end
 
    for i = 1:length(dataObj)
        if strcmpi(dataType{i}, 'scatter')
            dataObj(i).MarkerFaceColor = color(i,:);
            dataObj(i).MarkerEdgeColor = color(i,:);
        elseif strcmpi(dataType{i}, 'bar') || strcmpi(dataType{i}, 'histogram')
            dataObj(i).LineStyle = 'none';
            dataObj(i).FaceColor = color(i,:);
        elseif strcmpi(dataType{i}, 'line') || strcmpi(dataType{i}, 'errorbar')
            dataObj(i).LineWidth = 5;
            dataObj(i).Color = color(i,:);
            dataObj(i).Marker = '.';
            dataObj(i).MarkerFaceColor = color(i,:);
            dataObj(i).MarkerEdgeColor = color(i,:);
            dataObj(i).MarkerSize = 30;
            if strcmpi(dataType{i}, 'errorbar')
                %insert errorbar specific things here.
            end
        end
    end
    
    set(ax, 'TickLength',[0.01 0.01],...
        'FontSize', 20, 'FontName', 'Myriad Pro', 'FontWeight', 'bold');
    ax.TickDir = 'out';
    ax.LineWidth = axesLineWidth;
    faceColor = [255,251,206]/255; %yellow axis face.
    ax.Color = faceColor;
    fig.Color = [255,255,255]/255; %white figure background
    
end