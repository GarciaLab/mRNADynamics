function choice = chooseProjection

    d = dialog('Position',[300 300 250 150],'Name','Select a Projection');
    txt = uicontrol('Parent',d,'Style','text',...
           'Position',[20 80 210 40],'String','Select a projection');
       
    popup = uicontrol('Parent',d,'Style','popup','Position',[75 70 100 25],...
           'String',{'None (Default)';'Max Z';'Median Z'; 'Max Z and Time'},...
           'Callback',@popup_callback);
       
    btn = uicontrol('Parent',d,'Position',[89 20 70 25],'String','Okay',...
           'Callback','delete(gcf)');
       
    choice = 'None';
       
    % Wait for d to close before running to completion
    uiwait(d);
   
       function popup_callback(popup,event)
          idx = popup.Value;
          popup_items = popup.String;
          % This code uses dot notation to get properties.
          choice = char(popup_items(idx,:));
       end
end