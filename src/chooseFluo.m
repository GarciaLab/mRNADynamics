% chooseFlag.m
% author: Gabriella Martini
% date created: 9/13/20
% date last modified: 9/13/20
function FluoLabel= chooseFluo(CNfieldnames)
    FluoStrings = {};

    for i=1:length(CNfieldnames)
        fn = CNfieldnames{i};
        if regexp(fn, '^Fluo[A-za-z_0-9]*_TimeTrace')
            FluoStrings{length(FluoStrings)+1} = strrep(fn, '_TimeTrace', '');
        end
    end
    
    FluoStrings = reshape(FluoStrings, [length(FluoStrings), 1]);
    
    

    d = dialog('Position',[300 300 250 150],'Name','Select a Fluo Field');
    txt = uicontrol('Parent',d,'Style','text',...
           'Position',[20 80 210 40],'String','Select a Fluo Field');

    popup = uicontrol('Parent',d,'Style','popup','Position',[75 70 100 25],...
           'String',FluoStrings,...
           'Callback',@popup_callback);

    btn = uicontrol('Parent',d,'Position',[89 20 70 25],'String','Okay',...
           'Callback','delete(gcf)');

    FluoLabel = 'None';

    % Wait for d to close before running to completion
    uiwait(d);
    
       function popup_callback(popup,event)
          idx = popup.Value;
          fluofield = idx;
          popup_items = popup.String;
          % This code uses dot notation to get properties.
          FluoLabel = char(popup_items(idx,:));
       end



end