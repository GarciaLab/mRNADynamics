% chooseFlag.m
% author: Gabriella Martini
% date created: 9/13/20
% date last modified: 9/13/20
function [flag, choice] = chooseFlag

    d = dialog('Position',[300 300 250 150],'Name','Select a flag');
    txt = uicontrol('Parent',d,'Style','text',...
           'Position',[20 80 210 40],'String','Select a flag');
       
    popup = uicontrol('Parent',d,'Style','popup','Position',[75 70 100 25],...
           'String',{'Flag 0 (None)';'Flag 1 (auto-flagged)';'Flag 2 (nucleus disappears early in cycle) ';...
           'Flag 3 (nuclear appears late in cycle)';'Flag 4 (Missing chunk of frames mid-way)';...
           'Flag 5 (Nucleus disappers from z-stack)'; 'Flag 6 (Sick nucleus)'; 'Flag 7 (other)'},...
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