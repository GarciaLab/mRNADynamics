function standardizeFigure(ax, leg, varargin)

try
    colorDict = struct();
    colorDict.magenta = [208,109,171]/256;
    colorDict.lightBlue = [115,142,193]/256;
    colorDict.yellow = [234,194,100]/256;
    colorDict.red = [213,108,85]/256;
    colorDict.brown = [207 178 147] /256;
    colorDict.cyan = [108,188,233]/256;
    colorDictFields = fields(colorDict);
    
    color(1,:) = [0 0 0];
    
    axesLineWidth = .5;
    
    if isempty(leg)
        fig = gcf;
        leg= findobj(fig, 'Type', 'Legend');
    end
    
    dataObj = get(ax, 'Children');
    dataType = get(dataObj, 'Type');
    
    if ~iscell(dataType)
        dataType = {dataType};
    end
    
    legendSize = 8;
    fontSize = 8;
    
    for i = 1:length(varargin)
        if strcmpi(varargin{i}, 'axeslinewidth')
            axesLineWidth = varargin{i+1};
        elseif strcmpi(varargin{i}, 'red')
            color(i,:) = colorDict.red;
        elseif strcmpi(varargin{i}, 'yellow')
            color(i,:) = colorDict.yellow;
        elseif strcmpi(varargin{i}, 'cyan')
            color(i,:) = colorDict.cyan;
        elseif strcmpi(varargin{i}, 'magenta')
            color(i,:) = colorDict.magenta;
        elseif strcmpi(varargin{i}, 'lightBlue')
            color(i,:) = colorDict.lightBlue;
        elseif strcmpi(varargin{i}, 'brown')
            color(i,:) = colorDict.brown;
        elseif strcmpi(varargin{i}, 'legendFontSize')
            legendSize = varargin{i+1};
        elseif strcmpi(varargin{i}, 'fontSize')
            fontSize = varargin{i+1};
        end
    end
    
    if ~isempty(leg)
        for i = 1:length(leg)
            leg(i).FontSize = legendSize;
            leg(i).Box = 'off';
        end
    end
    
    for i = 1:length(dataObj)
        if strcmpi(dataType{i}, 'scatter')
            dataObj(i).Marker = '.';
            dataObj(i).SizeData = 8;
            %Change color to physical biology colors as long as the number
            %of colors needed is less than 5.
            if i <= length(colorDictFields)
                %                 dataObj(i).Color = colorDict.(colorDictFields{i});
                dataObj(i).MarkerFaceColor = colorDict.(colorDictFields{i});
                dataObj(i).MarkerEdgeColor = 'none';
            end
        elseif strcmpi(dataType{i}, 'bar') || strcmpi(dataType{i}, 'histogram')
            %             dataObj(i).LineStyle = 'none';
            if i <= length(colorDictFields)
                dataObj(i).FaceColor = colorDict.(colorDictFields{i});
            end
        elseif strcmpi(dataType{i}, 'line') || strcmpi(dataType{i}, 'errorbar')
            dataObj(i).LineWidth = .5;
            dataObj(i).Marker = 'o';
            dataObj(i).MarkerSize = 8;
            %Change color to physical biology colors as long as the number
            %of colors needed is less than 5.
            if i <= length(colorDictFields)
                dataObj(i).Color = colorDict.(colorDictFields{i});
                dataObj(i).MarkerFaceColor = colorDict.(colorDictFields{i});
                dataObj(i).MarkerEdgeColor = 'none';
            end
            if strcmpi(dataType{i}, 'errorbar')
                %insert errorbar specific things here.
            end
        end
    end
    
    set(ax, 'TickLength',[.01 .01],...
        'FontSize', fontSize, 'FontName', 'Lucida Sans OT', 'FontWeight', 'bold');
    ax.TickDir = 'in';
    ax.LineWidth = axesLineWidth;
    faceColor = 'none'; %yellow axis face.
    ax.Color = faceColor;
    ax.Box = 'on';
    
catch
    disp('couldn''t standardize figure');
end

end