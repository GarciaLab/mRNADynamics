% checkNuclearTracking_drawGUI.m
%author: Gabriella Martini
% date created: 9/7/20
% date last modified: 9/7/20
function [OverlayFig, overlayAxes, snippetFigAxes, rawDataAxes,...
    gaussianAxes, traceFig, traceFigAxes, zProfileFigAxes,...
    zTraceAxes, HisOverlayFig,HisOverlayFigAxes, multiFig]...
    ...
    = checkNuclearTracking_drawGUI(UseHistoneOverlay,ExperimentType,xSize, ySize)
%% Description
% This script is for 

OverlayFig = figure;
traceFig = figure;
multiFig = [];
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

% Define the overlayAxes
overlayAxes = axes(OverlayFig);
%tb = axtoolbar(overlayAxes);

% Separate the traceFig
traceFigAxes = axes(traceFig);
xlabel(traceFigAxes,'frame')
title(traceFigAxes, '', 'Interpreter', 'none');
%     traceFigAxes.Title.Interpreter = 'none';
yyaxis(traceFigAxes,'left')
% creating legend

str1 = 'Median Protein Fluo';
str2 = 'Max Protein Fluo';
str3 = 'Mid Median Protein Fluo';

%initialize curves
e1 = errorbar(traceFigAxes,[0, 1], [0, 1], [1, 1], 'k.-');
hold(traceFigAxes, 'on')
e2 = errorbar(traceFigAxes,[0, 1], [0, 1], [1, 1], 'b.-');
e3 = errorbar(traceFigAxes,[0, 1], [0, 1], [1, 1], 'g.-');
ylabel(traceFigAxes,'integrated intensity (a.u.)')
traceFigAxes.YAxis(2).Visible = 'off';
traceLeg = legend(traceFigAxes,[e1, e2, e3], str1,str2,str3, 'AutoUpdate', 'off', 'HandleVisibility', 'off');



zFig = figure;

zProfileFigAxes = subplot(1, 2, 1, 'Parent', zFig);
zTraceAxes = subplot(1, 2, 2, 'Parent', zFig);
ylabel(zProfileFigAxes,'intensity(au)', 'FontSize',12);
xlabel(zProfileFigAxes,'z-slice', 'FontSize',12);
xlabel(zTraceAxes,'frame')
ylabel(zTraceAxes,'z-slice')
title(zTraceAxes,'brightest Z trace')


snipFig = figure();
% snippetFigAxes = subplot(1, 3, 1, 'Parent', snipFig);
snippetFigAxes = axes(snipFig);
%   rawDataAxes = subplot(1, 3, 2, 'Parent', snipFig);
%   gaussianAxes = subplot(1, 3, 3, 'Parent', snipFig);

if UseHistoneOverlay
    set(HisOverlayFig, 'units', 'normalized', 'position', [0.01, 0.1, .33, .33]);
end

% Define the size of the figures/subplots
set(snipFig, 'units', 'normalized', 'position', [0.355, 0.05, 3 * (.2 / 2), .33 / 2]);
set(zFig, 'units', 'normalized', 'position', [0.67, 0.05, .2, .33 / 2]);



end