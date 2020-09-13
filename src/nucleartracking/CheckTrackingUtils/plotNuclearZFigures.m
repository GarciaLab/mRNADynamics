function [zFig, zFigAxes, zProfileHandles, zTraceHandles] = plotNuclearZFigures(...
    zFig, zFigAxes, zProfileHandles, zTraceHandles,...
    cntState)
zProfileFigAxes = zFigAxes{1};
zTraceAxes = zFigAxes{2};
z = 1:cntState.ZSlices;
fr_idx = find(cntState.Frames == cntState.CurrentFrame);
ZProfile = cntState.schnitzcells(cntState.CurrentNucleus).Fluo(fr_idx,:);

if isempty(zProfileHandles{1})
    if ~isempty(z)
        zProfileHandles{1} = plot(zProfileFigAxes, z, ZProfile, 'r.-');
    end
else
    set(zProfileHandles{1}, 'XData',z,...
            'YData',ZProfile, 'Color', 'r', 'LineStyle', '-', 'Marker', '.');    
end
hold(zProfileFigAxes, 'on')
if isempty(zProfileHandles{2})
    if ~isempty(z(cntState.MidMedZ(fr_idx)))
        zProfileHandles{2} = plot(zProfileFigAxes,z(cntState.MidMedZ(fr_idx)),...
    ZProfile(cntState.MidMedZ(fr_idx)), 'ko');
    end
else
    set(zProfileHandles{2}, 'XData',z(cntState.MidMedZ(fr_idx)),...
            'YData',ZProfile(cntState.MidMedZ(fr_idx)),...
            'Color', 'k', 'Marker', 'o');    
end
if isempty(zProfileHandles{3})
    if ~isempty(z(cntState.MazZ(fr_idx)))
        zProfileHandles{3} = plot(zProfileFigAxes,z(cntState.MaxZ(fr_idx)),...
    ZProfile(cntState.MaxZ(fr_idx)), 'bo');
    end
else
    set(zProfileHandles{3}, 'XData',z(cntState.MaxZ(fr_idx)),...
            'YData',ZProfile(cntState.MaxZ(fr_idx)),...
            'Color', 'b',  'Marker', 'o');    
end
if isempty(zProfileHandles{4})
    if ~isempty(z(cntState.MedZ(fr_idx)))
        zProfileHandles{4} = plot(zProfileFigAxes,z(cntState.MedZ(fr_idx)),...
    ZProfile(cntState.MedZ(fr_idx)), 'go');
    end
else
    set(zProfileHandles{4}, 'XData',z(cntState.MedZ(fr_idx)),...
            'YData',ZProfile(cntState.MedZ(fr_idx)),...
            'Color', 'g', 'Marker', 'o');    
end

zProfileFigYLimits = [0, max(ZProfile)*1.2];
try
    setPlotsInvisible(zProfileFigAxes);
    ylim(zProfileFigAxes, [0, zProfileFigYLimits(2)]);
    setPlotsVisible(zProfileFigAxes);
end


if isempty(zTraceHandles{1})
    if ~isempty(cntState.Frames)
        zTraceHandles{1} = plot(zTraceAxes,cntState.Frames,...
            cntState.MidMedZ, 'k.-');
    end
else
    set(zTraceHandles{1}, 'XData',cntState.Frames,...
            'YData',cntState.MidMedZ,...
            'Color', 'k', 'LineStyle', '-', 'Marker', '.');  
end
hold(zTraceAxes,'on')
if isempty(zTraceHandles{2})
    if ~isempty(cntState.Frames)
        zTraceHandles{2} = plot(zTraceAxes, cntState.Frames,...
            cntState.MaxZ,'.-b');
    end
else
    set(zTraceHandles{2}, 'XData',cntState.Frames,...
            'YData',cntState.MaxZ,...
            'Color', 'b', 'LineStyle', '-', 'Marker', '.');  
end
if isempty(zTraceHandles{3})
    if ~isempty(cntState.Frames)
        zTraceHandles{3} = plot(zTraceAxes, cntState.Frames,...
            cntState.MedZ,'.-g');
    end
else
    set(zTraceHandles{3}, 'XData',cntState.Frames,...
            'YData',cntState.MedZ,...
            'Color', 'g', 'LineStyle', '-', 'Marker', '.');  
end
if isempty(zTraceHandles{4})
    if ~isempty(cntState.Frames(fr_idx))
        zTraceHandles{4} = xline(cntState.Frames(fr_idx),...
            'Color', 'r', 'LineStyle', '--');
    end
else
    set(zTraceHandles{4}, 'value',cntState.Frames(fr_idx),...
            'Color', 'r', 'LineStyle', '--');  
end
hold(zTraceAxes, 'off')

zTraceFigYLimits = [0, cntState.ZSlices*1.3];
try
    setPlotsInvisible(zTraceAxes);
    ylim(zTraceAxes, [0, zTraceFigYLimits(2)]);
    setPlotsVisible(zTraceAxes);
end
zFigAxes = {zProfileFigAxes, zTraceAxes};



end
