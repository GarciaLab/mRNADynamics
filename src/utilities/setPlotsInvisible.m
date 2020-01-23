function setPlotsInvisible(ax)

axChild   = get(ax,'Children');
allCurves = findobj(axChild,'type','errorbar');
set(allCurves,'visible','off'); 

end