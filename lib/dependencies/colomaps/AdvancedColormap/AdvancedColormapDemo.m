
FH = 600;
fh = figure('Units','Pixels','Position',[0 0 FH FH], 'ToolBar','none', 'Color','w', 'UserData',1.0);
movegui(gcf,'center');
imagesc(peaks(256));
mm = 32;
ah = gca;
AW = FH - 2*mm - 150;
AH = FH - 2*mm;
set(gca,'Units','Pixels','Position',[150+mm 0+mm AW AH ]);
axis equal;
axis tight;
colorbar('east');
colormap jet;

PauseTime = 0.5;
drawnow;                                        fprintf('%s: Figure displyed\n',mfilename);
AdvancedColormap('AddColormapControls');        fprintf('%s: ColormapControls added\n',mfilename);
drawnow; pause(PauseTime);

button = questdlg({'Would you like to see full demo (click NO for short one)?','','This may take some time'});
switch lower(button)
    case 'yes'
        fprintf('%s: Proceeding with full demo mode\n',mfilename);
        h1 = uicontrol(gcf, ...
            'Style','pushbutton', ...
            'String','Faster', ...
            'Unit','pixels', ...
            'Position',[05 FH-20 55 17], ...
            'CallBack','set(gcf,''UserData'',get(gcf,''UserData'')/2);');
        h2 = uicontrol(gcf, ...
            'Style','pushbutton', ...
            'String','Slower', ...
            'Unit','pixels', ...
            'Position',[65 FH-20 55 17], ...
            'CallBack','set(gcf,''UserData'',get(gcf,''UserData'')*2);');

        cmaps = AdvancedColormap('GetSupportedColormaps');
        Nmaps = length(cmaps);
        cmap  = AdvancedColormap;
        [r,c] = size(cmap);
        wbh = waitbar(0,'', 'Name',mfilename);
        movegui(wbh,'west');
        ttt = 0;
        for k=1:Nmaps
            cmapName = cmaps{k};
            if any(strcmpi(cmapName,{'colorcube','flag','lines','prism','white'}))
                continue;
            end
            
            try
                waitbar(k/Nmaps,wbh,sprintf('Colormap #%d of %d: "%s"',k,Nmaps,strrep(cmapName,'_','\_')));
            catch ME
                delete(h1);
                delete(h2);
                return;
            end

            if ttt>0
                set(wbh, 'Name',sprintf('%s : ~%.1f seconds left',mfilename,ttt*(Nmaps-k+1)));
            end

            tic;
            AdvancedColormap(cmapName);        title(cmapName, 'Interpreter','none');               drawnow; pause(get(gcf,'UserData'));
            AdvancedColormap('Invert');        title('inverted', 'Interpreter','none');             drawnow; pause(get(gcf,'UserData'));
            AdvancedColormap('Reverse');       title('inverted + reversed', 'Interpreter','none');  drawnow; pause(get(gcf,'UserData'));
            AdvancedColormap('Invert');        title('reversed', 'Interpreter','none');             drawnow; pause(get(gcf,'UserData'));
            AdvancedColormap('Brighten',0.25); title('brightened', 'Interpreter','none');           drawnow; pause(get(gcf,'UserData'));
            AdvancedColormap('Darken',0.5) ;   title('darkened', 'Interpreter','none');             drawnow; pause(get(gcf,'UserData'));
            AdvancedColormap(cmapName);        title(cmapName, 'Interpreter','none');               drawnow; pause(get(gcf,'UserData'));
            
            if ~ishandle(h1),    break; end
            switch mod(k,7)
                case 0
                    title('shifting up');
                    for p=1:r
                        if ~ishandle(h1),    break; end
                        AdvancedColormap('ShiftUp'); drawnow;
                    end
                    pause(get(gcf,'UserData'));
                case 1
                    title('shifting down');
                    for p=1:r
                        if ~ishandle(h1),    break; end
                        AdvancedColormap('ShiftDown'); drawnow;
                    end
                    pause(get(gcf,'UserData'));
                case 2
                    title('RGB>RBG');    AdvancedColormap('RGB>RBG'); drawnow; pause(get(gcf,'UserData'));
                case 3
                    title('RGB>BGR');    AdvancedColormap('RGB>BGR'); drawnow; pause(get(gcf,'UserData'));
                case 4
                    title('RGB>GRB');    AdvancedColormap('RGB>GRB'); drawnow; pause(get(gcf,'UserData'));
                case 5
                    title('RGB>BRG');    AdvancedColormap('RGB>BRG'); drawnow; pause(get(gcf,'UserData'));
                    title('RGB>BRG x2'); AdvancedColormap('RGB>BRG'); drawnow; pause(get(gcf,'UserData'));
                case 6
                    title('RGB>GBR');    AdvancedColormap('RGB>GBR'); drawnow; pause(get(gcf,'UserData'));
                    title('RGB>GBR x2'); AdvancedColormap('RGB>GBR'); drawnow; pause(get(gcf,'UserData'));
            end
            if ~ishandle(h1),    break; end
            ttt = toc;
        end
        if ishandle(h1),    delete(h1); end
        if ishandle(h2),    delete(h2); end
        if ishandle(wbh),   close(wbh); end
    case 'no'
        fprintf('%s: Proceeding with brief demo mode\n',mfilename);
        AdvancedColormap('hot');                        fprintf('%s: hot colormap applied\n',mfilename);
        drawnow; pause(PauseTime);
        AdvancedColormap('togglecolormapcontrols');     fprintf('%s: ColormapControls toggled\n',mfilename);
        drawnow; pause(PauseTime);
        AdvancedColormap('jet');                        fprintf('%s: jet colormap applied\n',mfilename);
        drawnow; pause(PauseTime);
        AdvancedColormap('ShowColormapControls');       fprintf('%s: ColormapControls displayed\n',mfilename);
        drawnow; pause(PauseTime);
        AdvancedColormap('Topographic');                fprintf('%s: Topographic colormap applied\n',mfilename);
        drawnow; pause(PauseTime);
        AdvancedColormap('hot');                        fprintf('%s: hot colormap applied\n',mfilename);
        drawnow; pause(PauseTime);
        AdvancedColormap('Thermal');                    fprintf('%s: Thermal colormap applied\n',mfilename);
        drawnow; pause(PauseTime);
        cmap  = colormap; %AdvancedColormap('GetCurrentColormap');
        [r,c] = size(cmap);
        for k=1:r
            AdvancedColormap('ShiftUp');
            drawnow;
        end
        msgbox('Now you can play with ADVANCEDCOLORMAP controls by yourself ;o)','Play','none','modal');
        %button = questdlg({'Would you like to see full demo (click NO for short one)?','','This may take some time'});
    case 'cancel'
        fprintf('%s: Canceled by user\n',mfilename);
end
