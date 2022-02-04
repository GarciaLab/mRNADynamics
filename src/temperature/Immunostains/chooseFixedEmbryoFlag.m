% chooseFixedEmbryoFlag.m
% author: Gabriella Martini
% date created: 9/2/21
% date last modified: 9/2/21
function [flag, choice] = chooseFixedEmbryoFlag

    d = dialog('Position',[300 300 250 150],'Name','Select a flag');
    txt = uicontrol('Parent',d,'Style','text',...
           'Position',[20 80 210 40],'String','Select a flag');
       
    popup = uicontrol('Parent',d,'Style','popup','Position',[75 70 100 25],...
           'String',{'Flag 0 (None)';'Flag 1 (Gastrulated)';'Flag 2 (Deformed or cut) ';...
           'Flag 3 (too young or unfertilized)';'Flag 4 (Embryos clumped)';...
           'Flag 5 (DV axis oriented incorrectly)'; 'Flag 6 (Bubble in the way)'; 'Flag 7 (Mounted with uneven anterior and posterior)';...
           'Flag 8 (Embryo out of focus)'; 'Flag 9 (other)'}, ...
           'Callback',@popup_callback);
       
    btn = uicontrol('Parent',d,'Position',[89 20 70 25],'String','Okay',...
           'Callback','delete(gcf)');
       
    choice = 'None';
       
    % Wait for d to close before running to completion
    uiwait(d);
       function popup_callback(popup,event)
          idx = popup.Value;
          flag = idx-1;
          popup_items = popup.String;
          % This code uses dot notation to get properties.
          choice = char(popup_items(idx,:));
       end
end