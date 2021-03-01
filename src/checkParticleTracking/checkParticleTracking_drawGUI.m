function [OverlayFig, overlayAxes, snippetFigAxes, rawDataAxes,...
    gaussianAxes, traceFig, traceFigAxes, zProfileFigAxes,...
    zTraceAxes, HisOverlayFig,HisOverlayFigAxes, multiFig]...
    ...
    = checkParticleTracking_drawGUI(UseHistoneOverlay,...
    fish, plot3DGauss, ExperimentType, multiView, xSize, ySize, PlotInputChannel)
%% Description
% This script is for W

OverlayFig = figure;
traceFig = figure;
multiFig = [];
HisOverlayFig = [];
HisOverlayFigAxes = [];
traceFigAxes = [];
zTraceAxes = [];
rawDataAxes = [];
gaussianAxes = [];

if ~exist('PlotInputChannel', 'var')
    PlotInputChannel = false;
end



if UseHistoneOverlay
    HisOverlayFig = figure;
    HisOverlayFigAxes = axes(HisOverlayFig);
end

% Define the overlayAxes
overlayAxes = axes(OverlayFig);
%tb = axtoolbar(overlayAxes);

% Separate the traceFig
if ~fish
    traceFigAxes = axes(traceFig);
    xlabel(traceFigAxes,'frame')
    title(traceFigAxes, '', 'Interpreter', 'none');
    %     traceFigAxes.Title.Interpreter = 'none';
    yyaxis(traceFigAxes,'left')
    % creating legend
    if plot3DGauss
        str1 = '3-slice mRNA';
        str2 = '3D-Gaussian fit mRNA';
    else
        str1 = '1-slice mRNA';
        str2 = 'multi-slice mRNA';
    end
    %initialize curves
    e1 = errorbar(traceFigAxes,[0, 1], [0, 1], [1, 1], 'k.-');
    hold(traceFigAxes, 'on')
    e2 = errorbar(traceFigAxes,[0, 1], [0, 1], [1, 1], 'b.-');
    ylabel(traceFigAxes,'integrated intensity (a.u.)')
    if strcmpi(ExperimentType, 'inputoutput')
        yyaxis(traceFigAxes,'right')
        e3 = plot(traceFigAxes,[0, 1], [0, 1], 'r.-', 'DisplayName', 'protein');
        ylabel(traceFigAxes,'input protein intensity (a.u.)');
        traceLeg = legend(traceFigAxes,[e1, e2, e3], str1,str2, 'protein', 'AutoUpdate', 'off', 'HandleVisibility', 'off');
    elseif strcmpi(ExperimentType, '2spot')
        if plot3DGauss
            str3 = 'Twin 3-slice mRNA';
            str4 = 'Twin 3D-Gaussian fit mRNA';
        else
            str3 = 'Twin 1-slice mRNA';
            str4 = 'Twin multi-slice mRNA';
        end
        e3 = errorbar(traceFigAxes,[0, 1], [0, 1], [1, 1], 'g.-');
        e4 = errorbar(traceFigAxes,[0, 1], [0, 1], [1, 1], 'm.-');
        
        if PlotInputChannel
            str5 = 'protein';
            yyaxis(traceFigAxes,'right')
            e5 = plot(traceFigAxes,[0, 1], [0, 1], 'r.-', 'DisplayName', 'protein');
            ylabel(traceFigAxes,'input protein intensity (a.u.)');
            traceFigAxes.YAxis(2).Visible = 'off';
            traceLeg = legend(traceFigAxes,[e1, e2, e3, e4, e5], str1,str2,str3, str4, str5, 'AutoUpdate', 'off', 'HandleVisibility', 'off', 'Location', 'eastoutside');
        else
            traceFigAxes.YAxis(2).Visible = 'off';
            traceLeg = legend(traceFigAxes,[e1, e2, e3, e4], str1,str2,str3, str4, 'AutoUpdate', 'off', 'HandleVisibility', 'off', 'Location', 'eastoutside');
        end
    else
        traceFigAxes.YAxis(2).Visible = 'off';
        traceLeg = legend(traceFigAxes,[e1, e2], str1,str2, 'AutoUpdate', 'off', 'HandleVisibility', 'off', 'Location', 'eastoutside');
    end
    
    hold(traceFigAxes, 'off')
end

zFig = figure;
if ~fish & ~strcmpi(ExperimentType, '2spot')
    zProfileFigAxes = subplot(1, 2, 1, 'Parent', zFig);
    zTraceAxes = subplot(1, 2, 2, 'Parent', zFig);
    ylabel(zProfileFigAxes,'intensity(au)', 'FontSize',12);
    xlabel(zProfileFigAxes,'z-slice', 'FontSize',12);
    xlabel(zTraceAxes,'frame')
    ylabel(zTraceAxes,'z-slice')
    title(zTraceAxes,'brightest Z trace')
elseif ~fish & strcmpi(ExperimentType, '2spot')
    zProfileFigAxes{1} = subplot(1, 3, 1, 'Parent', zFig);
    p1 = plot(zProfileFigAxes{1},[0, 1], [0, 1],'b.-');
    hold(zProfileFigAxes{1}, 'on')
    p2 = plot(zProfileFigAxes{1},[0, 1], [0, 1],  'm.-');
    c1 = plot(zProfileFigAxes{1},[0], [0],  'bo');
    c2 = plot(zProfileFigAxes{1},[0], [0],  'mo');
    zProfileFigAxes{2} = subplot(1, 3, 2, 'Parent', zFig);
    p3 = plot(zProfileFigAxes{2},[0, 1], [0, 1], 'b.-');
    hold(zProfileFigAxes{2}, 'on')
    p4 = plot(zProfileFigAxes{2},[0, 1], [0, 1],  'm.-');
    c3 = plot(zProfileFigAxes{2},[0], [0],  'bo');
    c4 = plot(zProfileFigAxes{2},[0], [0],  'mo');
    zTraceAxes = subplot(1, 3, 3, 'Parent', zFig);
    p5 = plot(zTraceAxes,[0, 1], [0, 1],  'b.-');
    hold(zTraceAxes, 'on')
    p6 = plot(zTraceAxes,[0, 1], [0, 1], 'm.-');
    c5 = plot(zTraceAxes,[0], [0],  'bo');
    c6 = plot(zTraceAxes,[0], [0],  'mo');
    ylabel(zProfileFigAxes{1},'x (pixels)', 'FontSize',12);
    xlabel(zProfileFigAxes{1},'frame', 'FontSize',12);
    ylabel(zProfileFigAxes{2},'y (pixels)', 'FontSize',12);
    xlabel(zProfileFigAxes{2},'frame', 'FontSize',12);
    ylabel(zTraceAxes,'z-slice', 'FontSize',12);
    xlabel(zTraceAxes,'frame', 'FontSize',12);
else
    zProfileFigAxes = axes(zFig);
end

snipFig = figure();
% snippetFigAxes = subplot(1, 3, 1, 'Parent', snipFig);
snippetFigAxes = axes(snipFig);
%   rawDataAxes = subplot(1, 3, 2, 'Parent', snipFig);
%   gaussianAxes = subplot(1, 3, 3, 'Parent', snipFig);

if UseHistoneOverlay
    set(HisOverlayFig, 'units', 'normalized', 'position', [0.01, 0.1, .33, .33]);
end

% Define the size of the figures/subplots
if ~fish
    %overlayDim = [.82, .45];
    % Define the size of the overlayFig using the format of the
    % movie/image (xSize, ySize)
    % The figure looks inconvenient if the height is smaller than 0.4, so
    % we will fix that height, and adjust the width depending on the ratio
    % of xSize and ySize.
    yDim = 0.45;
    ratio = xSize/ySize;
    xDim = min([0.6 yDim*ratio]);
    % maxSize = max(xSize, ySize);
    overlayDim = [xDim yDim];
    set(OverlayFig, 'units', 'normalized', 'OuterPosition', [0,1-overlayDim(2), overlayDim(1), overlayDim(2)]);
    %set(OverlayFig, 'units', 'normalized', 'position', [0,1-overlayDim(2), overlayDim(1), overlayDim(2)]);
    set(traceFig, 'units', 'normalized', 'position', [overlayDim(1)+0.05, 0.6, .4 .3])
    %set(overlayAxes, 'units', 'normalized', 'position', [-.25 .06 .9 .9])
    %set(traceFigAxes, 'units', 'normalized', 'position', [.48 .17 .48 .63])
    set(snipFig, 'units', 'normalized', 'position', [0.355, 0.15, 3 * (.2 / 2), .33 / 2]);
    if ~strcmpi(ExperimentType, '2spot')
        set(zFig, 'units', 'normalized', 'position', [0.67, 0.15, .2, .33 / 2]);
    else
        set(zFig, 'units', 'normalized', 'position', [0.45, 0.15, .55, .3]);
    end
else
    set(snipFig, 'units', 'normalized', 'position', [0.355, 0.05, 3 * (.2 / 2), .33 / 2]);
    set(zFig, 'units', 'normalized', 'position', [0.67, 0.05, .2, .33 / 2]);
end


if multiView
    
    multiFig = figure('units', 'normalized', 'Position', [.6,.08, .8*overlayDim(2), overlayDim(2)]);
end

end