function StandardFigurePBoC(PlotHandle,AxisHandle)

%Make it such that it can identify things like '-o'.



for i=1:length(PlotHandle)

    %Line style
    
    if ~strcmpi(get(PlotHandle(i),'Type'),'scatter') && ~strcmpi(get(PlotHandle(i),'Type'),'image')
        LineStyle=get(PlotHandle(i),'LineStyle');
        if (~isempty(strmatch(LineStyle,'-')))|...
                (~isempty(strmatch(LineStyle,'--')))|...
                (~isempty(strmatch(LineStyle,'-.')))|...
                (~isempty(strmatch(LineStyle,'-o')))
            set(PlotHandle(i),'LineWidth',1)
        end
    end
    

    %Marker style. Do this only if we're not dealing with a histogram, bar
    %plot, image, or area plot
    if ~strcmp(get(PlotHandle(i),'Type'),'histogram') && ...
                        ~strcmpi(get(PlotHandle(i),'Type'),'image') && ...
                        ~strcmp(get(PlotHandle(i),'Type'),'bar') && ...
                        ~strcmpi(get(PlotHandle(i),'Type'),'area')
        if ~strcmp(get(PlotHandle(i),'Marker'),'none')
            if get(PlotHandle(i),'Marker')=='.'
                set(PlotHandle(i),'MarkerSize',15)
                set(PlotHandle(i),'LineWidth',1)
            elseif (get(PlotHandle(i),'Marker')=='d')|...
                (strcmp(get(PlotHandle(i),'Marker'),'square'))|...
                (get(PlotHandle(i),'Marker')=='o')
                if strcmpi(get(PlotHandle(i),'Type'),'scatter')
                    set(PlotHandle(i),'SizeData',5)
                else
                    set(PlotHandle(i),'MarkerSize',5)
                end
                set(PlotHandle(i),'LineWidth',1)
                if ~isempty(strmatch(get(PlotHandle(i),'MarkerFaceColor'),'none'))
                    set(PlotHandle(i),'MarkerFaceColor','w')
                end
            end
        end
    end

    %Color section
    if strcmp(get(PlotHandle(i),'Type'),'hggroup')
        try ChangeColorPBoC2(PlotHandle(i),'FaceColor')
        end
        try ChangeColorPBoC2(PlotHandle(i),'Color')
        end
    elseif strcmp(get(PlotHandle(i),'Type'),'errorbarseries')
    
    %New Matlab histogram
    elseif strcmp(get(PlotHandle(i),'Type'),'histogram')
        %Change the colors of the bars. Note that this will only work if
        %we've specified the colors. If we let Matlab choose colors, then
        %this function won't change them.
        ChangeColorPBoC2(PlotHandle(i),'FaceColor')
    elseif strcmp(get(PlotHandle(i),'Type'),'bar')
        %Change the colors of the bars. Note that this will only work if
        %we've specified the colors. If we let Matlab choose colors, then
        %this function won't change them.
        ChangeColorPBoC2(PlotHandle(i),'FaceColor')
    else
        if ~strcmpi(get(PlotHandle(i),'Type'),'image')
            %Why do I have this auto thingy?
            if isempty(strmatch(get(PlotHandle(i),'Color'),'auto'))
                ChangeColorPBoC2(PlotHandle(i),'Color')
            end

            if isempty(strmatch(get(PlotHandle(i),'MarkerEdgeColor'),'auto'))
                ChangeColorPBoC2(PlotHandle(i),'MarkerEdgeColor')
            end

            if isempty(strmatch(get(PlotHandle(i),'MarkerFaceColor'),'auto'))&...
                    isempty(strmatch(get(PlotHandle(i),'MarkerFaceColor'),'none'))
                ChangeColorPBoC2(PlotHandle(i),'MarkerFaceColor')
            end
        end
    end
    
end

%This was for the first edition:
% set(AxisHandle,...
%     'Box','off',...
%     'Color',[0.8510,0.8510,0.8510],...
%     'XColor','w','YColor','w','ZColor','w',...
%     'TickLength',[0.02,0.05])
%Second edition
set(AxisHandle,...
    'Box','off',...
    'Color',[228,221,209]/255,...
    'XColor','w','YColor','w','ZColor','w',...
    'TickLength',[0.02,0.05])

%Added this so that the log tick marks wouldn't go away.
set(gca,'XTickMode','manual')
set(gca,'YTickMode','manual')


set(get(AxisHandle,'XLabel'),'FontSize',15,'FontName','Lucida Sans')
set(get(AxisHandle,'YLabel'),'FontSize',15,'FontName','Lucida Sans')
set(get(AxisHandle,'ZLabel'),'FontSize',15,'FontName','Lucida Sans')
set(AxisHandle,'FontSize',15,'FontName','Lucida Sans')

box on