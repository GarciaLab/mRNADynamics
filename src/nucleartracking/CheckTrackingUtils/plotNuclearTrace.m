function [traceFig, traceFigAxes,traceHandles] = ...
    plotNuclearTrace(traceFig, traceFigAxes,traceHandles, cntState,...
    anaphaseInMins, ElapsedTime,...
    ncFrames,Prefix)
%PLOTTRACE
%plot traces in checkparticletracking
hasInputChannels = ~isempty(cntState.liveExperiment.inputChannels);
ncPresent =  cntState.schnitzcells(cntState.CurrentNucleus).cycle;
priorAnaphaseInMins = anaphaseInMins(ncPresent(1)-8); %min  subtracts 8 because the first element corresponds to nc 9
priorAnaphase = ncFrames(ncPresent(1)-8); %frame

if ncPresent < 14
    nextAnaphase =  ncFrames(ncPresent(1)-7);
else
    nextAnaphase = length(cntState.Ellipses);
end

numFrames = length(cntState.Ellipses);
nucleusFirstTimePoint = ElapsedTime(...
    cntState.schnitzcells(cntState.CurrentNucleus).frames(1)); %min

%traceFigTimeAxis = ElapsedTime(cntState.schnitzcells(cntState.CurrentNucleus).frames) - nucleusFirstTimePoint; %min
traceFigTimeAxis = cntState.schnitzcells(cntState.CurrentNucleus).frames.';
approvedNucleiFrames = logical(cntState.schnitzcells(cntState.CurrentNucleus).FrameApproved);
if hasInputChannels
    if isempty(traceHandles{1})
        if ~isempty(cntState.Frames(approvedNucleiFrames))
            set(traceHandles{1},'Visible','on')
            traceHandles{1} = plot(traceFigAxes,traceFigTimeAxis(approvedNucleiFrames),...
                cntState.MidMedFluo(approvedNucleiFrames), 'k.-');
        else
            set(traceHandles{1},'Visible','off') 
        end
    elseif ~isempty(cntState.Frames(approvedNucleiFrames))
        set(traceHandles{1},'Visible','on')
        set(traceHandles{1}, 'XData', cntState.Frames(approvedNucleiFrames),...
                'YData', cntState.MidMedFluo(approvedNucleiFrames),...
                'Color', 'k', 'LineStyle', '-');

    else
       set(traceHandles{1},'Visible','off') 
    end

    hold(traceFigAxes, 'on')

    if isempty(traceHandles{2})
        if ~isempty(cntState.Frames(approvedNucleiFrames))
            set(traceHandles{2},'Visible','on')
            traceHandles{2} = plot(traceFigAxes,cntState.Frames(approvedNucleiFrames),...
                cntState.MaxFluo(approvedNucleiFrames), 'b.-');
        else
            set(traceHandles{2},'Visible','off') 
        end
    elseif ~isempty(cntState.Frames(approvedNucleiFrames))
        set(traceHandles{2},'Visible','on')
        set(traceHandles{2}, 'XData', cntState.Frames(approvedNucleiFrames),...
                'YData', cntState.MaxFluo(approvedNucleiFrames),...
                'Color', 'b', 'LineStyle', '-');
    else
        set(traceHandles{2},'Visible','off') 

    end

    if isempty(traceHandles{3})
        if ~isempty(cntState.Frames(approvedNucleiFrames))
            set(traceHandles{3},'Visible','on')
            traceHandles{3} = plot(traceFigAxes,cntState.Frames(approvedNucleiFrames),...
                cntState.MedFluo(approvedNucleiFrames), 'g.-');
        else
            set(traceHandles{3},'Visible','off')
        end
    elseif ~isempty(cntState.Frames(approvedNucleiFrames))
        set(traceHandles{3},'Visible','on')
        set(traceHandles{3}, 'XData', cntState.Frames(approvedNucleiFrames),...
                'YData', cntState.MedFluo(approvedNucleiFrames),...
            'Color', 'g', 'LineStyle', '-');
    else
        set(traceHandles{3},'Visible','off')
    end
    if isempty(traceHandles{4})
        if ~isempty(traceFigTimeAxis(~approvedNucleiFrames))
            set(traceHandles{4},'Visible','on')
            traceHandles{4} = plot(traceFigAxes,traceFigTimeAxis(~approvedNucleiFrames),...
                cntState.MidMedFluo(~approvedNucleiFrames), 'r.');
        else
            set(traceHandles{4},'Visible','off')
        end
    elseif ~isempty(traceFigTimeAxis(~approvedNucleiFrames))
        set(traceHandles{4},'Visible','on')
        set(traceHandles{4}, 'XData', cntState.Frames(~approvedNucleiFrames),...
                'YData', cntState.MidMedFluo(~approvedNucleiFrames),...
            'Color', 'r', 'Marker', '.');  
    else
        set(traceHandles{4},'Visible','off')
    end
    set(get(get(traceHandles{4}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');


    if isempty(traceHandles{5})
        if ~isempty(traceFigTimeAxis(cntState.Frames==cntState.CurrentFrame))
            set(traceHandles{5},'Visible','on')
            traceHandles{5} =  plot(traceFigAxes,...
                traceFigTimeAxis(cntState.Frames==cntState.CurrentFrame),...
            cntState.MidMedFluo(cntState.Frames==cntState.CurrentFrame),'ok');
        else
            set(traceHandles{5},'Visible','off')
        end
    elseif ~isempty(traceFigTimeAxis(cntState.Frames==cntState.CurrentFrame))
        set(traceHandles{5},'Visible','on')
        set(traceHandles{5}, 'XData', traceFigTimeAxis(cntState.Frames==cntState.CurrentFrame),...
                'YData',cntState.MidMedFluo(cntState.Frames==cntState.CurrentFrame));    
    else
        set(traceHandles{5},'Visible','off')
    end
    set(get(get(traceHandles{5}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');


    if isempty(traceHandles{6})
        if ~isempty(traceFigTimeAxis(~approvedNucleiFrames))
            set(traceHandles{6},'Visible','on')
            traceHandles{6} =  plot(traceFigAxes,...
                traceFigTimeAxis(~approvedNucleiFrames),...
            cntState.MaxFluo(~approvedNucleiFrames),'.r');
        else
           set(traceHandles{6},'Visible','off') 
        end
    elseif ~isempty(traceFigTimeAxis(~approvedNucleiFrames))
        set(traceHandles{6},'Visible','on')
        set(traceHandles{6}, 'XData',traceFigTimeAxis(~approvedNucleiFrames),...
                'YData',cntState.MaxFluo(~approvedNucleiFrames));  
    else
        set(traceHandles{6},'Visible','off')
    end

    set(get(get(traceHandles{6}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');


    if isempty(traceHandles{7})
        if ~isempty(traceFigTimeAxis(cntState.Frames==cntState.CurrentFrame))
            set(traceHandles{7},'Visible','on')
            traceHandles{7} =  plot(traceFigAxes,...
                traceFigTimeAxis(cntState.Frames==cntState.CurrentFrame),...
                cntState.MaxFluo(cntState.Frames==cntState.CurrentFrame),'ob');
        else
            set(traceHandles{7},'Visible','off')
        end
    elseif ~isempty(traceFigTimeAxis(cntState.Frames==cntState.CurrentFrame))
        set(traceHandles{7},'Visible','on')
        set(traceHandles{7}, 'XData',traceFigTimeAxis(cntState.Frames==cntState.CurrentFrame),...
                'YData',cntState.MaxFluo(cntState.Frames==cntState.CurrentFrame));  
    else
        set(traceHandles{7},'Visible','off')
    end

    set(get(get(traceHandles{7}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');

    if isempty(traceHandles{8})
        if ~isempty(traceFigTimeAxis(~approvedNucleiFrames))
            set(traceHandles{8},'Visible','on')
            traceHandles{8} = plot(traceFigAxes,...
                traceFigTimeAxis(~approvedNucleiFrames),...
                cntState.MedFluo(~approvedNucleiFrames),'.r');
        else
            set(traceHandles{8},'Visible','off')
        end
    elseif ~isempty(traceFigTimeAxis(~approvedNucleiFrames))
        set(traceHandles{8},'Visible','on')
        set(traceHandles{8}, 'XData',traceFigTimeAxis(~approvedNucleiFrames),...
                'YData',cntState.MedFluo(~approvedNucleiFrames));  
    else
        set(traceHandles{8},'Visible','off')
    end
    set(get(get(traceHandles{8}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');


    if isempty(traceHandles{9})
        if ~isempty(traceFigTimeAxis(cntState.Frames==cntState.CurrentFrame))
            set(traceHandles{9},'Visible','on')
            traceHandles{9} = plot(traceFigAxes,...
                traceFigTimeAxis(cntState.Frames==cntState.CurrentFrame),...
                cntState.MedFluo(cntState.Frames==cntState.CurrentFrame),...
                'og');
        else
            set(traceHandles{9},'Visible','off')
        end
    elseif ~isempty(traceFigTimeAxis(cntState.Frames==cntState.CurrentFrame))
        set(traceHandles{9},'Visible','on')
        set(traceHandles{9}, 'XData',traceFigTimeAxis(cntState.Frames==cntState.CurrentFrame),...
                'YData',cntState.MedFluo(cntState.Frames==cntState.CurrentFrame));
    else
        set(traceHandles{9},'Visible','off')
    end


    set(get(get(traceHandles{9}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');

    xlabel(traceFigAxes,'frame')
    title(traceFigAxes, '', 'Interpreter', 'none');
    %     traceFigAxes.Title.Interpreter = 'none';
    % creating legend

    str1 = 'Mid Median';
    str2 = 'Max';
    str3 = 'Median';


    %initialize curves

    ylabel(traceFigAxes,'integrated intensity (a.u.)')
    traceLeg = legend(traceFigAxes,[traceHandles{1}, traceHandles{2}, traceHandles{3}], str1,str2,str3, 'AutoUpdate', 'off', 'HandleVisibility', 'off');
    try
        setPlotsInvisible(traceFigAxes);
        xlim(traceFigAxes,double([min([priorAnaphase, min(traceFigTimeAxis)]),max([max(traceFigTimeAxis), nextAnaphase])])+[-1,1]);
        setPlotsVisible(traceFigAxes);
    end

    traceFigYLimits = [0, max(cntState.MaxFluo)*1.2];
    tr_counter = 10;
else
    if isempty(traceHandles{1})
        if ~isempty(cntState.Frames(approvedNucleiFrames))
            set(traceHandles{1},'Visible','on')
            traceHandles{1} = plot(traceFigAxes,cntState.Frames(approvedNucleiFrames),...
                ones(1,length(cntState.Frames(approvedNucleiFrames))), 'k.-');
        else
            set(traceHandles{1},'Visible','off') 
        end
    elseif ~isempty(traceFigTimeAxis(approvedNucleiFrames))
        set(traceHandles{1},'Visible','on')
        set(traceHandles{1}, 'XData', cntState.Frames(approvedNucleiFrames),...
                'YData', ones(1,length(cntState.Frames(approvedNucleiFrames))),...
                'Color', 'k', 'LineStyle', '-');

    else
       set(traceHandles{1},'Visible','off') 
    end

    hold(traceFigAxes, 'on')
    if isempty(traceHandles{2})
        if ~isempty(traceFigTimeAxis(~approvedNucleiFrames))
            set(traceHandles{2},'Visible','on')
            traceHandles{2} = plot(traceFigAxes,traceFigTimeAxis(~approvedNucleiFrames),...
                ones(1,length(cntState.Frames(~approvedNucleiFrames))), 'r.');
        else
            set(traceHandles{2},'Visible','off')
        end
    elseif ~isempty(traceFigTimeAxis(~approvedNucleiFrames))
        set(traceHandles{2},'Visible','on')
        set(traceHandles{2}, 'XData', cntState.Frames(~approvedNucleiFrames),...
                'YData', ones(1,length(cntState.Frames(~approvedNucleiFrames))),...
            'Color', 'r', 'Marker', '.');  
    else
        set(traceHandles{2},'Visible','off')
    end
    set(get(get(traceHandles{2}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');


    if isempty(traceHandles{3})
        if ~isempty(traceFigTimeAxis(cntState.Frames==cntState.CurrentFrame))
            set(traceHandles{3},'Visible','on')
            traceHandles{3} =  plot(traceFigAxes,...
                traceFigTimeAxis(cntState.Frames==cntState.CurrentFrame),...
            ones(1,1),'ok');
        else
            set(traceHandles{3},'Visible','off')
        end
    elseif ~isempty(traceFigTimeAxis(cntState.Frames==cntState.CurrentFrame))
        set(traceHandles{3},'Visible','on')
        set(traceHandles{3}, 'XData', traceFigTimeAxis(cntState.Frames==cntState.CurrentFrame),...
                'YData', ones(1,1));    
    else
        set(traceHandles{3},'Visible','off')
    end
    set(get(get(traceHandles{3}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    try
        setPlotsInvisible(traceFigAxes);
        xlim(traceFigAxes,double([min([priorAnaphase, min(traceFigTimeAxis)]),max([max(traceFigTimeAxis), nextAnaphase])])+[-1,1]);
        setPlotsVisible(traceFigAxes);
    end
    traceFigYLimits = [0, 2];
    tr_counter = 4;
end

for i = 1:length(ncFrames)

    currentAnaphaseBoundary = ncFrames(i);
    if isempty(traceHandles{tr_counter})
        if ~isempty(ones(1,2).*double(currentAnaphaseBoundary))
            set(traceHandles{tr_counter},'Visible','on')    
            traceHandles{tr_counter} = plot(traceFigAxes,ones(1,2).*double(currentAnaphaseBoundary),...
                traceFigYLimits,...
                'LineWidth',2,'Color','black', 'Marker', 'none');
        else
          set(traceHandles{tr_counter},'Visible','off')     
        end
    elseif ~isempty(ones(1,2).*double(currentAnaphaseBoundary))
        set(traceHandles{tr_counter},'Visible','on')   
        set(traceHandles{tr_counter}, 'XData',ones(1,2).*double(currentAnaphaseBoundary),...
                'YData',traceFigYLimits);
    else
        set(traceHandles{tr_counter},'Visible','off')   
    end
    

    set(get(get(traceHandles{tr_counter}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    tr_counter = tr_counter + 1;
    
end

% plotting the lines and traces

hold(traceFigAxes, 'off')

try
    setPlotsInvisible(traceFigAxes);
    ylim(traceFigAxes, [0, traceFigYLimits(2)]);
    setPlotsVisible(traceFigAxes);
end



hold(traceFigAxes,'off')
% 
% idata1 = amp1(approvedParticleFrames);
% idata2 = amp2(approvedParticleFrames);
% 
% 

% creating axis title
frame_idx = find(cntState.Frames == cntState.CurrentFrame);
numNuclei = length(cntState.schnitzcells);
firstLine = [Prefix];
% STOPPED HERE
secondLine = ['Nucleus: ',num2str(cntState.CurrentNucleus),'/',num2str(numNuclei), '   Frame: ',num2str(cntState.CurrentFrame),'/',num2str(numFrames),'    ',num2str(round(ElapsedTime(cntState.CurrentFrame))), 'm'];
if hasInputChannels
    thirdLine = ['Z: ',num2str(cntState.MidMedZ(frame_idx)),'/',num2str(cntState.ZSlices)];
end
    
if isfield(cntState.FrameInfo, 'nc')
    if hasInputChannels
        axisTitle={firstLine,...
            [secondLine,'    (nc',num2str(cntState.FrameInfo(cntState.CurrentFrame).nc),')'],...
            thirdLine};
    else
        axisTitle={firstLine,...
            [secondLine,'    (nc',num2str(cntState.FrameInfo(cntState.CurrentFrame).nc),')']};
    end
else
    if hasInputChannels
        axisTitle={firstLine,secondLine,thirdLine};
    else
        axisTitle={firstLine,secondLine};
    end
end

if cntState.HideApprovedFlag == 1
    axisTitle=[axisTitle,', Showing non-flagged particles'];
elseif cntState.HideApprovedFlag == 2
    axisTitle=[axisTitle,', Showing disapproved particles'];
end
hold(traceFigAxes, 'off');
setPlotsInvisible(traceFigAxes);
set(traceFigAxes.Title,'String', axisTitle);
setPlotsVisible(traceFigAxes);

set(traceFig,'Name',['Current Nucleus Trace'])

end

