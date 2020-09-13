function [zFig, zFigAxes, zProfileHandles, zTraceHandles] = plotNuclearZFigures(...
    zFig, zFigAxes, zProfileHandles, zTraceHandles,...
    cntState)
zProfileFigAxes = zFigAxes{1};
zTraceAxes = zFigAxes{2};
z = 1:cntState.ZSlices;
fr_idx = find(cntState.Frames == cntState.CurrentFrame);
ZProfile = cntState.schnitzcells(cntState.CurrentNucleus).Fluo(fr_idx,:);

if isempty(zProfileHandles{1})
    if ~isempty(ZProfile)
        set(zProfileHandles{1},'Visible','on')
        zProfileHandles{1} = plot(zProfileFigAxes, z, ZProfile, 'r.-');
    else
        zProfileHandles{1} = plot(zProfileFigAxes, [0, 1], [0, 1], 'r.');
        set(zProfileHandles{1},'Visible','off')
    end
elseif ~isempty(ZProfile)
    set(zProfileHandles{1},'Visible','on')
    set(zProfileHandles{1}, 'XData',z,...
            'YData',ZProfile, 'Color', 'r', 'LineStyle', '-', 'Marker', '.'); 
else
    zProfileHandles{1} = plot(zProfileFigAxes, [0, 1], [0, 1], 'r.');
    set(zProfileHandles{1},'Visible','off') 
end
hold(zProfileFigAxes, 'on')
if isempty(zProfileHandles{2})
    if ~isempty(ZProfile(cntState.MidMedZ(fr_idx)))
        set(zProfileHandles{2},'Visible','on')
        zProfileHandles{2} = plot(zProfileFigAxes,z(cntState.MidMedZ(fr_idx)),...
            ZProfile(cntState.MidMedZ(fr_idx)), 'ko');
    else
        zProfileHandles{2} = plot(zProfileFigAxes, [0, 1], [0, 1], 'ko');
        set(zProfileHandles{2},'Visible','off')
    end
elseif ~isempty(ZProfile(cntState.MidMedZ(fr_idx)))
    set(zProfileHandles{2},'Visible','on') 
    set(zProfileHandles{2}, 'XData',z(cntState.MidMedZ(fr_idx)),...
            'YData',ZProfile(cntState.MidMedZ(fr_idx)),...
            'Color', 'k', 'Marker', 'o');    
else
    zProfileHandles{2} = plot(zProfileFigAxes, [0, 1], [0, 1], 'ko');
    set(zProfileHandles{2},'Visible','off')
end
if isempty(zProfileHandles{3})
    if ~isempty(ZProfile(cntState.MaxZ(fr_idx)))
        set(zProfileHandles{3},'Visible','on')
        zProfileHandles{3} = plot(zProfileFigAxes,z(cntState.MaxZ(fr_idx)),...
            ZProfile(cntState.MaxZ(fr_idx)), 'bo');
    else
        zProfileHandles{3} = plot(zProfileFigAxes, [0, 1], [0, 1], 'bo');
        set(zProfileHandles{3},'Visible','off')
    end
elseif ~isempty(ZProfile(cntState.MaxZ(fr_idx)))
    set(zProfileHandles{3},'Visible','on')
    set(zProfileHandles{3}, 'XData',z(cntState.MaxZ(fr_idx)),...
            'YData',ZProfile(cntState.MaxZ(fr_idx)),...
            'Color', 'b',  'Marker', 'o');    
else
    zProfileHandles{3} = plot(zProfileFigAxes, [0, 1], [0, 1], 'bo');
    set(zProfileHandles{3},'Visible','off') 
end
if isempty(zProfileHandles{4})
    if ~isempty(ZProfile(cntState.MedZ(fr_idx)))
        set(zProfileHandles{4},'Visible','on')
        zProfileHandles{4} = plot(zProfileFigAxes,z(cntState.MedZ(fr_idx)),...
            ZProfile(cntState.MedZ(fr_idx)), 'go');
    else
        zProfileHandles{4} = plot(zProfileFigAxes, [0, 1], [0, 1], 'go');
        set(zProfileHandles{4},'Visible','off')
    end
elseif ~isempty(ZProfile(cntState.MedZ(fr_idx)))
    set(zProfileHandles{4},'Visible','on')
    set(zProfileHandles{4}, 'XData',z(cntState.MedZ(fr_idx)),...
            'YData',ZProfile(cntState.MedZ(fr_idx)),...
            'Color', 'g', 'Marker', 'o');   
else
    zProfileHandles{4} = plot(zProfileFigAxes, [0, 1], [0, 1], 'go');
    set(zProfileHandles{4},'Visible','off')

end
if ~isempty(fr_idx)
    zProfileFigYLimits = [0, max(ZProfile)*1.2];
end
if exist('zProfileFigYLimits', 'var')
    setPlotsInvisible(zProfileFigAxes);
    ylim(zProfileFigAxes, [0, max([1, zProfileFigYLimits(2)])]);
    setPlotsVisible(zProfileFigAxes);
end


if isempty(zTraceHandles{1})
    if ~isempty(cntState.Frames)
        set(zTraceHandles{1},'Visible','on')
        zTraceHandles{1} = plot(zTraceAxes,cntState.Frames,...
            cntState.MidMedZ, 'k.-');
    else
        zTraceHandles{1} = plot(zTraceAxes, [0, 1], [0, 1], 'k.-');
        set(zTraceHandles{1},'Visible','off')
    end
elseif ~isempty(cntState.Frames)
    set(zTraceHandles{1},'Visible','on')
    set(zTraceHandles{1}, 'XData',cntState.Frames,...
            'YData',cntState.MidMedZ,...
            'Color', 'k', 'LineStyle', '-', 'Marker', '.'); 
else
    zTraceHandles{1} = plot(zTraceAxes, [0, 1], [0, 1], 'k.-');
    set(zTraceHandles{1},'Visible','off')
end
hold(zTraceAxes,'on')
if isempty(zTraceHandles{2})
    if ~isempty(cntState.Frames)
        set(zTraceHandles{2},'Visible','on')
        zTraceHandles{2} = plot(zTraceAxes, cntState.Frames,...
            cntState.MaxZ,'.-b');
    else
        zTraceHandles{2} = plot(zTraceAxes, [0, 1], [0, 1], 'b.-');
        set(zTraceHandles{2},'Visible','off')
    end
elseif ~isempty(cntState.Frames)
    set(zTraceHandles{2},'Visible','on')
    set(zTraceHandles{2}, 'XData',cntState.Frames,...
            'YData',cntState.MaxZ,...
            'Color', 'b', 'LineStyle', '-', 'Marker', '.');  
else
    zTraceHandles{2} = plot(zTraceAxes, [0, 1], [0, 1], 'b.-');
    set(zTraceHandles{2},'Visible','off')
end
if isempty(zTraceHandles{3})
    if ~isempty(cntState.Frames)
        set(zTraceHandles{3},'Visible','on')
        zTraceHandles{3} = plot(zTraceAxes, cntState.Frames,...
            cntState.MedZ,'.-g');
    else
        zTraceHandles{3} = plot(zTraceAxes, [0, 1], [0, 1], 'g.-');
        set(zTraceHandles{3},'Visible','off')
    end
elseif ~isempty(cntState.Frames)
    set(zTraceHandles{3},'Visible','on')
    set(zTraceHandles{3}, 'XData',cntState.Frames,...
            'YData',cntState.MedZ,...
            'Color', 'g', 'LineStyle', '-', 'Marker', '.');  
else
    zTraceHandles{3} = plot(zTraceAxes, [0, 1], [0, 1], 'g.-');
    set(zTraceHandles{3},'Visible','off')
end

if isempty(zTraceHandles{4})
    if ~isempty(cntState.Frames(fr_idx))
        set(zTraceHandles{4},'Visible','on')
        zTraceHandles{4} = xline(cntState.Frames(fr_idx),...
            'Color', 'r', 'LineStyle', '--');
    else
        zTraceHandles{4} = xline(zTraceAxes, 0, 'Color', 'r', 'LineStyle', '--');
        set(zTraceHandles{4},'Visible','off')
    end
elseif ~isempty(cntState.Frames(fr_idx))
    set(zTraceHandles{4},'Visible','on')
    set(zTraceHandles{4}, 'value',cntState.Frames(fr_idx),...
            'Color', 'r', 'LineStyle', '--');  
else
    zTraceHandles{4} = xline(zTraceAxes, 0, 'Color', 'r', 'LineStyle', '--');
    set(zTraceHandles{4},'Visible','off')
end
hold(zTraceAxes, 'off')

zTraceFigYLimits = [0, cntState.ZSlices*1.3];
if exist('zTraceFigYLimits', 'var')
    setPlotsInvisible(zTraceAxes);
    ylim(zTraceAxes, [0, max([zTraceFigYLimits(2), 1])]);
    setPlotsVisible(zTraceAxes);
end
zFigAxes = {zProfileFigAxes, zTraceAxes};

disp('here')

end
