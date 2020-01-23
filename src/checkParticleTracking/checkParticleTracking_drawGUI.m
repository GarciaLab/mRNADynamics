function [Overlay, overlayAxes, snippetFigAxes, rawDataAxes, gaussianAxes, traceFigAxes, zProfileFigAxes,...
    zTraceAxes, HisOverlayFig,HisOverlayFigAxes] = checkParticleTracking_drawGUI(UseHistoneOverlay, fish, plot3DGauss, ExperimentType)


Overlay = figure;
HisOverlayFig = [];
HisOverlayFigAxes = [];
traceFigAxes = [];
zTraceAxes = [];
rawDataAxes = [];
gaussianAxes = [];



if UseHistoneOverlay
    HisOverlayFig = figure;
    HisOverlayFigAxes = axes(HisOverlayFig);
end


if ~fish
    overlayAxes = subplot(1, 2, 1, 'Parent', Overlay);
    tb = axtoolbar(overlayAxes);
    traceFigAxes = subplot(1, 2, 2, 'Parent', Overlay);
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
    e1 = errorbar(traceFigAxes,1, 1, 1, 'k');
    hold(traceFigAxes, 'on')
    e2 = errorbar(traceFigAxes,1, 1, 1, 'b');
    ylabel(traceFigAxes,'integrated intensity (a.u.)')
    if strcmpi(ExperimentType, 'inputoutput')
        yyaxis(traceFigAxes,'right')
        e3 = errorbar(traceFigAxes,1, 1, 1, 'r.-', 'DisplayName', 'protein');
        ylabel(traceFigAxes,'input protein intensity (a.u.)');
        traceLeg = legend(traceFigAxes,[e1, e2, e3], str1,str2, 'protein', 'AutoUpdate', 'off');
    else
        traceFigAxes.YAxis(2).Visible = 'off';
        traceLeg = legend(traceFigAxes,[e1, e2], str1,str2, 'AutoUpdate', 'off');
    end
else
    overlayAxes = axes(Overlay);
end

zFig = figure;
if ~fish
    zProfileFigAxes = subplot(1, 2, 1, 'Parent', zFig);
    zTraceAxes = subplot(1, 2, 2, 'Parent', zFig);
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

if ~fish
    set(Overlay, 'units', 'normalized', 'position', [0.01, .45, .82, .33]);
    set(overlayAxes, 'units', 'normalized', 'position', [-.25 .06 .9 .9])
    set(traceFigAxes, 'units', 'normalized', 'position', [.48 .17 .48 .63])
    set(snipFig, 'units', 'normalized', 'position', [0.355, 0.15, 3 * (.2 / 2), .33 / 2]);
    set(zFig, 'units', 'normalized', 'position', [0.67, 0.15, .2, .33 / 2]);
else
    set(snipFig, 'units', 'normalized', 'position', [0.355, 0.05, 3 * (.2 / 2), .33 / 2]);
    set(zFig, 'units', 'normalized', 'position', [0.67, 0.05, .2, .33 / 2]);
end

end