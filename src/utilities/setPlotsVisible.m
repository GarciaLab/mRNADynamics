function setPlotsVisible(ax)

axChild   = get(ax,'Children');
allCurves = findobj(axChild,'type','errorbar');
set(allCurves,'visible','on'); 

end