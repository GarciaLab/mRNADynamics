% checkNuclearTracking_drawGUI.m
%author: Gabriella Martini
% date created: 9/7/20
% date last modified: 10/13/20
function [OverlayFig, overlayAxes, ImageHandle,...
    snipFig, snipAxes, snipHandles,...
    traceFig, traceFigAxes,traceHandles,...
    zFig, zFigAxes, zProfileHandles, zTraceHandles]...
    ...
    = checkNuclearTracking_drawGUI(cntState,...
    Prefix, xSize, ySize, ImageHisMat,ElapsedTime,...
    anaphaseInMins, ncFrames)
%% Description
% This script is for
hasInputChannels = ~isempty(cntState.liveExperiment.inputChannels);


close all
OverlayFig = figure(1);
snipFig = figure(2);
traceFig = figure(3);
if hasInputChannels
    zFig = figure(4);
else
    zFig = [];
end

multiFig = [];
HisOverlayFig = [];
HisOverlayFigAxes = [];
traceFigAxes = [];
zTraceAxes = [];
rawDataAxes = [];
gaussianAxes = [];

%
%% Setup Figure 1
hisImage = ImageHisMat(:, :, cntState.CurrentFrame);
numFrames = length(cntState.Ellipses);
ApprovedNuclei = [cntState.schnitzcells.Approved];
CheckedNuclei = [cntState.schnitzcells.Approved];
% Define the overlayAxes
overlayAxes = axes(OverlayFig);
if cntState.UseHistoneOverlay
    if hasInputChannels
        if isempty(cntState.DisplayRange)
            HisOverlayImageMat=cat(3,mat2gray(hisImage),mat2gray(cntState.MaxImageMat),zeros(size(cntState.MaxImageMat)));
        else
            HisOverlayImageMat=cat(3,mat2gray(hisImage,double(cntState.DisplayRange)),mat2gray(cntState.MaxImageMat),zeros(size(cntState.MaxImageMat)));
        end
    else
        if isempty(cntState.DisplayRange)
            HisOverlayImageMat=cat(3,mat2gray(hisImage),zeros(size(hisImage)),zeros(size(hisImage)));
        else
            HisOverlayImageMat=cat(3,mat2gray(hisImage,double(cntState.DisplayRange)),zeros(size(hisImage)),zeros(size(hisImage)));
        end
    end
    ImageHandle = imshow(HisOverlayImageMat,[],'Border','Tight','Parent',overlayAxes, 'InitialMagnification', 'fit');
else
    if hasInputChannels
        ImageHandle = imshow(cntState.MaxImageMat,...
            cntState.DisplayRangeSpot, 'Border', 'Tight', 'Parent', overlayAxes, ...
            'InitialMagnification', 'fit');
    else
        error('Must use histone overlay or protein channel inputs.')
    end
    
end
overlayAxes.XAxis.Visible = 'off';
overlayAxes.YAxis.Visible = 'off';
%%





hold(overlayAxes, 'on')


EllipseHandle = [];
EllipseHandleBlue = [];
EllipseHandleRed = [];
EllipseHandleMagenta = [];
EllipseHandleOrange = [];
EllipseHandleGreen = [];
EllipseHandleWhite = [];


schnitzCellNo_ApprovedUnchecked=[];
schnitzCellNo_ApprovedUncheckedPassed=[];
schnitzCellApproved=[];
schnitzCellCheckedRejected=[];
schnitzCellUncheckedRejected=[];
schnitzCellCheckedRejectedUnflagged = [];


for i=1:cntState.numNuclei()
    if ~isempty(cntState.CurrentNucleusCellNo) & (i == cntState.CurrentNucleus)
        EllipseHandleGreen = ellipseCellCNT(cntState, cntState.CurrentNucleusCellNo, 'g', 10, overlayAxes);
    elseif (cntState.schnitzcells(i).Checked ==0) & (cntState.schnitzcells(i).Approved ==1) & (cntState.schnitzcells(i).FirstPass == 0)
        frame_idx = find(cntState.schnitzcells(i).frames == cntState.CurrentFrame, 1);
        if ~isempty(frame_idx)
            schnitzCellNo_ApprovedUnchecked= [schnitzCellNo_ApprovedUnchecked,cntState.schnitzcells(i).cellno(frame_idx)];
        end
    elseif (cntState.schnitzcells(i).Checked ==0) & (cntState.schnitzcells(i).Approved ==1) & (cntState.schnitzcells(i).FirstPass == 1)
        frame_idx = find(cntState.schnitzcells(i).frames == cntState.CurrentFrame, 1);
        if ~isempty(frame_idx)
            schnitzCellNo_ApprovedUncheckedPassed= [schnitzCellNo_ApprovedUncheckedPassed,cntState.schnitzcells(i).cellno(frame_idx)];
        end
    elseif (cntState.schnitzcells(i).Approved ==1) & (cntState.schnitzcells(i).Checked ==1)
        frame_idx = find(cntState.schnitzcells(i).frames == cntState.CurrentFrame, 1);
        if ~isempty(frame_idx)
            schnitzCellApproved = [schnitzCellApproved,cntState.schnitzcells(i).cellno(frame_idx)];
        end
    elseif (cntState.schnitzcells(i).Approved ~=1) & (cntState.schnitzcells(i).Checked ==1) & (cntState.schnitzcells(i).Flag ~= 0)
        frame_idx = find(cntState.schnitzcells(i).frames == cntState.CurrentFrame, 1);
        if ~isempty(frame_idx)
            schnitzCellCheckedRejected = [schnitzCellCheckedRejected,cntState.schnitzcells(i).cellno(frame_idx)];
        end
    elseif (cntState.schnitzcells(i).Approved ~=1) & (cntState.schnitzcells(i).Checked ==1) & (cntState.schnitzcells(i).Flag ==0)
        frame_idx = find(cntState.schnitzcells(i).frames == cntState.CurrentFrame, 1);
        if ~isempty(frame_idx)
            schnitzCellCheckedRejectedUnflagged = [schnitzCellCheckedRejectedUnflagged,cntState.schnitzcells(i).cellno(frame_idx)];
        end
    elseif (cntState.schnitzcells(i).Approved ~=1) & (cntState.schnitzcells(i).Checked ==1)
        frame_idx = find(cntState.schnitzcells(i).frames == cntState.CurrentFrame, 1);
        if ~isempty(frame_idx)
            schnitzCellCheckedRejected = [schnitzCellCheckedRejected,cntState.schnitzcells(i).cellno(frame_idx)];
        end
    elseif (cntState.schnitzcells(i).Approved ~=1) & (cntState.schnitzcells(i).Checked ==0)
        frame_idx = find(cntState.schnitzcells(i).frames == cntState.CurrentFrame, 1);
        if ~isempty(frame_idx)
            schnitzCellUncheckedRejected = [schnitzCellUncheckedRejected,cntState.schnitzcells(i).cellno(frame_idx)];
        end
    end
end

EllipseHandle = notEllipseCellCNT(cntState, schnitzCellNo_ApprovedUnchecked, 'y', 10, overlayAxes);
EllipseHandleMagenta = notEllipseCellCNT(cntState, schnitzCellNo_ApprovedUncheckedPassed, 'm', 10, overlayAxes);
EllipseHandleBlue = notEllipseCellCNT(cntState, schnitzCellApproved, 'c', 10, overlayAxes);
EllipseHandleRed = notEllipseCellCNT(cntState, schnitzCellCheckedRejected, 'r', 10, overlayAxes);
EllipseHandleOrange = notEllipseCellCNT(cntState, schnitzCellUncheckedRejected, [1, 0.5, 0], 10, overlayAxes);
EllipseHandleWhite = notEllipseCellCNT(cntState, schnitzCellCheckedRejectedUnflagged, 'w', 10, overlayAxes);







if (ApprovedNuclei(cntState.CurrentNucleus) == 1) & (CheckedNuclei(cntState.CurrentNucleus) == 1)
    set(OverlayFig,'Color','g')
elseif ApprovedNuclei(cntState.CurrentNucleus) <= 0
    set(OverlayFig,'Color','r')
elseif (ApprovedNuclei(cntState.CurrentNucleus) == 2) & (CheckedNuclei(cntState.CurrentNucleus) == 1)
    set(OverlayFig,'Color','y')
else
    set(OverlayFig,'Color','default')
end





if isfield(cntState.FrameInfo, 'nc')
    if isfield(cntState.schnitzcells, 'anaphaseFrame')
        set(OverlayFig,'Name',['Schnitz Cell: ',num2str(cntState.CurrentNucleus),'/',num2str(cntState.numNuclei()),...
            ', Frame: ',num2str(cntState.CurrentFrame),'/',num2str(numFrames),...
            ', Ana: ',num2str(cntState.schnitzcells(cntState.CurrentNucleus).anaphaseFrame),...
            ' nc: ', num2str(cntState.FrameInfo(cntState.CurrentFrame).nc)]);
        
    elseif hasInputChannels
        set(OverlayFig,'Name',['Schnitz Cell: ',num2str(cntState.CurrentNucleus),'/',num2str(cntState.numNuclei()),...
            ', Frame: ',num2str(cntState.CurrentFrame),'/',num2str(numFrames),...
            ', Z: ',num2str(cntState.CurrentZ),'/',num2str(cntState.ZSlices),...
            ' nc: ', num2str(cntState.FrameInfo(cntState.CurrentFrame).nc)]);
    else
        set(OverlayFig,'Name',['Schnitz Cell: ',num2str(cntState.CurrentNucleus),'/',num2str(cntState.numNuclei()),...
            ', Frame: ',num2str(cntState.CurrentFrame),'/',num2str(numFrames),...
            ' nc: ', num2str(cntState.FrameInfo(cntState.CurrentFrame).nc)]);
    end
end

%% Figure 2

% snippetFigAxes = subplot(1, 3, 1, 'Parent', snipFig);
if hasInputChannels
    snipHistoneAxes = subplot(1, 4, 1, 'Parent', snipFig);
    snipMidMedInputAxes = subplot(1, 4, 2, 'Parent', snipFig);
    snipInputAxes = subplot(1, 4, 3, 'Parent', snipFig);
    snipMedInputAxes = subplot(1, 4, 4, 'Parent', snipFig);
    snipAxes = {snipHistoneAxes, snipMidMedInputAxes, snipInputAxes, snipMedInputAxes};
else
    snipHistoneAxes = axes(snipFig);
    snipAxes = {snipHistoneAxes};
end

%%

fr_idx = find(cntState.Frames == cntState.CurrentFrame);
if hasInputChannels
    maxz_idx = cntState.MaxZ(fr_idx);
    medz_idx = cntState.MedZ(fr_idx);
    midmedz_idx = cntState.MidMedZ(fr_idx);
end
pixelSize = cntState.FrameInfo(1).PixelSize; % (um)
nucleusDiameter = getDefaultParameters(cntState.FrameInfo,['d9']);
nucleusDiameter_pixels = nucleusDiameter/pixelSize;

scale = 2; %magnification of snippet
snippet_size = double(round(nucleusDiameter_pixels));

xSchnitz = cntState.getCurrentX();
ySchnitz = cntState.getCurrentY();

imSnippet = zeros(2*snippet_size + 1, 2*snippet_size + 1, 'double');

imMedSnippet = zeros(2*snippet_size + 1, 2*snippet_size + 1, 'double');
imMidMedSnippet = zeros(2*snippet_size + 1, 2*snippet_size + 1, 'double');
imHisSnippet = zeros(2*snippet_size + 1, 2*snippet_size + 1, 'double');

if (~isempty(xSchnitz)) & (~isempty(ySchnitz))
    if ySchnitz-snippet_size < 1
        rmin_idx = snippet_size-ySchnitz+2;
    else
        rmin_idx = 1;
    end
    if xSchnitz-snippet_size < 1
        cmin_idx = snippet_size-xSchnitz+2;
    else
        cmin_idx = 1;
    end
    if ySchnitz+snippet_size > ySize
        rmax_idx = rmin_idx + 2*snippet_size - (ySchnitz + snippet_size-ySize);
        %ySchnitz + snippet_size + 1 - ySize;
    else
        rmax_idx = size(imSnippet, 1);
    end
    if xSchnitz+snippet_size > xSize
        cmax_idx = cmin_idx + 2*snippet_size - (xSchnitz + snippet_size-xSize);
    else
        cmax_idx =size(imSnippet, 2);
    end
    if hasInputChannels
        imMidMedSnippet(rmin_idx:rmax_idx, cmin_idx:cmax_idx) = ...
            mat2gray( double(cntState.ImageMat(max(1,ySchnitz-snippet_size):min(ySize,ySchnitz+snippet_size),...
            max(1,xSchnitz-snippet_size):min(xSize,xSchnitz+snippet_size))));
        
        imSnippet(rmin_idx:rmax_idx, cmin_idx:cmax_idx) = ...
            mat2gray( double(cntState.MaxImageMat(max(1,ySchnitz-snippet_size):min(ySize,ySchnitz+snippet_size),...
            max(1,xSchnitz-snippet_size):min(xSize,xSchnitz+snippet_size))));
        
        imMedSnippet(rmin_idx:rmax_idx, cmin_idx:cmax_idx) = ...
            mat2gray( double(cntState.MedImageMat(max(1,ySchnitz-snippet_size):min(ySize,ySchnitz+snippet_size),...
            max(1,xSchnitz-snippet_size):min(xSize,xSchnitz+snippet_size))));
    end
    
    imHisSnippet(rmin_idx:rmax_idx, cmin_idx:cmax_idx) = ...
        mat2gray( double(hisImage(max(1,ySchnitz-snippet_size):min(ySize,ySchnitz+snippet_size),...
        max(1,xSchnitz-snippet_size):min(xSize,xSchnitz+snippet_size))));
end
IntegrationRadius = 2/pixelSize; % 6*ceil(sqrt(212/pixelSize)); %integrate 109 pixels around the spot with 212nm pixel size
[xGrid, yGrid] = meshgrid(1:2*snippet_size+1,1:2*snippet_size+1);
rGrid = sqrt((xGrid-ceil(snippet_size)).^2 + (yGrid-ceil(snippet_size)).^2);
IntegrationArea= rGrid < IntegrationRadius & (rGrid+1) >= IntegrationRadius;


if hasInputChannels
    SnippetOverlay=cat(3,IntegrationArea/2 + ...
        +imSnippet,imSnippet,imSnippet);
    SnippetMedOverlay=cat(3,IntegrationArea/2 + ...
        +imMedSnippet,imMedSnippet,imMedSnippet);
    SnippetMidMedOverlay=cat(3,IntegrationArea/2 + ...
        +imMidMedSnippet,imMidMedSnippet,imMidMedSnippet);
end
HisSnippetOverlay=cat(3,IntegrationArea/2 + ...
    +imHisSnippet,imHisSnippet,imHisSnippet);
%%
snipHandles = {};



if hasInputChannels
    scale = 25;
    snipHandles{1} = imshow(HisSnippetOverlay,...
        [],'Border','Tight','InitialMagnification',scale*2, 'Parent', snipHistoneAxes);
    title(snipHistoneAxes, 'Histone');
    snipHandles{2} = imshow(SnippetMidMedOverlay,...
        [],'Border','Tight','InitialMagnification',scale*2, 'Parent', snipMidMedInputAxes);
    title(snipMidMedInputAxes, 'Mid Median Z')
    snipHandles{3} = imshow(SnippetOverlay,...
        [],'Border','Tight','InitialMagnification',scale*2, 'Parent', snipInputAxes);
    title(snipInputAxes, 'Max Z')
    snipHandles{4} = imshow(SnippetMedOverlay,...
        [],'Border','Tight','InitialMagnification',scale*2, 'Parent', snipMedInputAxes);
    title(snipMedInputAxes, 'Median Z');
else
    snipHandles{1} = imshow(HisSnippetOverlay,...
        [],'Border','Tight','InitialMagnification',100, 'Parent', snipHistoneAxes);
    title(snipHistoneAxes, 'Histone');
end
set(snipFig,'Name',['Image Snippets'])


%% Figure 3


ncPresent =  cntState.schnitzcells(cntState.CurrentNucleus).cycle;
priorAnaphaseInMins = anaphaseInMins(ncPresent(1)-8); %min  subtracts 8 because the first element corresponds to nc 9
priorAnaphase = ncFrames(ncPresent(1)-8); %frame
if ncPresent < 14
    nextAnaphase =  ncFrames(ncPresent(1)-7);
else
    nextAnaphase = length(cntState.Ellipses);
end

nucleusFirstTimePoint = ElapsedTime(...
    cntState.schnitzcells(cntState.CurrentNucleus).frames(1)); %min

%traceFigTimeAxis = ElapsedTime(cntState.schnitzcells(cntState.CurrentNucleus).frames) - nucleusFirstTimePoint; %min
traceFigTimeAxis = cntState.schnitzcells(cntState.CurrentNucleus).frames.';
approvedNucleiFrames = cntState.schnitzcells(cntState.CurrentNucleus).FrameApproved;


traceFigAxes = axes(traceFig);
traceHandles = {};
if hasInputChannels
    if ~isempty(traceFigTimeAxis(approvedNucleiFrames))
        traceHandles{1} = plot(traceFigAxes,traceFigTimeAxis(approvedNucleiFrames),...
            cntState.MidMedFluo(approvedNucleiFrames), 'k.-');
        hold(traceFigAxes, 'on')
        traceHandles{2}  = plot(traceFigAxes,traceFigTimeAxis(approvedNucleiFrames),...
            cntState.MaxFluo(approvedNucleiFrames), 'b.-');
        traceHandles{3} = plot(traceFigAxes,traceFigTimeAxis(approvedNucleiFrames),...
            cntState.MedFluo(approvedNucleiFrames), 'g.-');
    else
        traceHandles{1} = plot(traceFigAxes, [0, 1], [0, 1], 'k.-');
        set(traceHandles{1},'Visible','off')
        hold(traceFigAxes, 'on')
        traceHandles{2} = plot(traceFigAxes, [0, 1], [0, 1], 'b.-');
        set(traceHandles{2},'Visible','off')
        traceHandles{3} = plot(traceFigAxes, [0, 1], [0, 1], 'g.-');
        set(traceHandles{3},'Visible','off')
    end
    
    if ~isempty(traceFigTimeAxis(~approvedNucleiFrames))
        traceHandles{4} = plot(traceFigAxes,traceFigTimeAxis(~approvedNucleiFrames),cntState.MidMedFluo(~approvedNucleiFrames),'.r');
        set(get(get(traceHandles{4}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    else
        traceHandles{4} = plot(traceFigAxes, [0, 1], [0, 1], 'r.');
        set(traceHandles{4},'Visible','off')
        set(get(get(traceHandles{4}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    end
    
    if ~isempty(traceFigTimeAxis(cntState.Frames==cntState.CurrentFrame))
        traceHandles{5} = plot(traceFigAxes,traceFigTimeAxis(cntState.Frames==cntState.CurrentFrame),...
            cntState.MidMedFluo(cntState.Frames==cntState.CurrentFrame),'ok');
        set(get(get(traceHandles{5}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    else
        traceHandles{5} = plot(traceFigAxes, [0, 1], [0, 1], 'ok');
        set(get(get(traceHandles{5}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
        set(traceHandles{5},'Visible','off')
    end
    
    if ~isempty(traceFigTimeAxis(~approvedNucleiFrames))
        traceHandles{6} = plot(traceFigAxes,traceFigTimeAxis(~approvedNucleiFrames),cntState.MaxFluo(~approvedNucleiFrames),'.r');
        set(get(get(traceHandles{6}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    else
        traceHandles{6} = plot(traceFigAxes, [0, 1], [0, 1], 'r.');
        set(traceHandles{6},'Visible','off')
        set(get(get(traceHandles{6}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    end
    
    if ~isempty(traceFigTimeAxis(cntState.Frames==cntState.CurrentFrame))
        traceHandles{7} = plot(traceFigAxes,traceFigTimeAxis(cntState.Frames==cntState.CurrentFrame),...
            cntState.MaxFluo(cntState.Frames==cntState.CurrentFrame),'ob');
        set(get(get(traceHandles{7}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    else
        traceHandles{7} =plot(traceFigAxes, [0, 1], [0, 1], 'ob');
        set(get(get(traceHandles{7}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
        set(traceHandles{7},'Visible','off')
    end
    
    if ~isempty(traceFigTimeAxis(~approvedNucleiFrames))
        traceHandles{8} = plot(traceFigAxes,traceFigTimeAxis(~approvedNucleiFrames),cntState.MedFluo(~approvedNucleiFrames),'.r');
        set(get(get(traceHandles{8}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    else
        traceHandles{8} = plot(traceFigAxes, [0, 1], [0, 1], 'r.');
        set(traceHandles{8},'Visible','off')
        set(get(get(traceHandles{8}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    end
    
    if ~isempty(traceFigTimeAxis(cntState.Frames==cntState.CurrentFrame))
        traceHandles{9} = plot(traceFigAxes,traceFigTimeAxis(cntState.Frames==cntState.CurrentFrame),...
            cntState.MedFluo(cntState.Frames==cntState.CurrentFrame),'og');
        set(get(get(traceHandles{9}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    else
        traceHandles{9} =plot(traceFigAxes, [0, 1], [0, 1], 'og');
        set(get(get(traceHandles{9}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
        set(traceHandles{9},'Visible','off')
    end
    
    
    
    
    
    
    xlabel(traceFigAxes,'frame')
    title(traceFigAxes, '', 'Interpreter', 'none');
    %     traceFigAxes.Title.Interpreter = 'none';
    % creating legend
    
    str1 = 'Mid Median';
    str2 = 'Max';
    str3 = 'Median';
    
    %initialize curves
    
    ylabel(traceFigAxes,'integrated intensity (a.u.)')
    traceLeg = legend(traceFigAxes,[traceHandles{1}, traceHandles{2}, traceHandles{3}], str1,str2,str3,...
        'AutoUpdate', 'off', 'HandleVisibility', 'off', 'location', 'southeast');
    nc = cntState.schnitzcells(cntState.CurrentNucleus).cycle;
    
    
    traceFigYLimits = [0, max(cntState.MaxFluo)*1.3];
    for i = 1:length(ncFrames)
        hndl_idx = length(traceHandles)+1;
        currentAnaphaseBoundary = ncFrames(i);
        traceHandles{hndl_idx} =plot(traceFigAxes,ones(1,2).*double(currentAnaphaseBoundary),traceFigYLimits,...
            'LineWidth',2,'Color','black', 'Marker', 'none');
        set(get(get(traceHandles{hndl_idx}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    end
    try
        setPlotsInvisible(traceFigAxes);
        xlim(traceFigAxes,double([min([priorAnaphase, min(traceFigTimeAxis)]),max([max(traceFigTimeAxis), nextAnaphase])])+[-1,1]);
        setPlotsVisible(traceFigAxes);
    end
    
    
    % plotting the lines and traces
    
    hold(traceFigAxes, 'off')
    
    try
        setPlotsInvisible(traceFigAxes);
        ylim(traceFigAxes, [0, traceFigYLimits(2)]);
        setPlotsVisible(traceFigAxes);
    end
    
    hold(traceFigAxes,'off')
else
    if ~isempty(traceFigTimeAxis(approvedNucleiFrames))
        traceHandles{1} = plot(traceFigAxes,traceFigTimeAxis(approvedNucleiFrames),...
            ones(1, length(traceFigTimeAxis(approvedNucleiFrames))), 'k.-');
        hold(traceFigAxes, 'on')
    else
        traceHandles{1} = plot(traceFigAxes, [0, 1], [0, 1], 'k.-');
        set(traceHandles{1},'Visible','off')
        hold(traceFigAxes, 'on')
    end
    if ~isempty(traceFigTimeAxis(~approvedNucleiFrames))
        traceHandles{2} = plot(traceFigAxes,traceFigTimeAxis(~approvedNucleiFrames),ones(1, length(traceFigTimeAxis(~approvedNucleiFrames))),'.r');
        set(get(get(traceHandles{2}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    else
        traceHandles{2} = plot(traceFigAxes, [0, 1], [0, 1], 'r.');
        set(traceHandles{2},'Visible','off')
        set(get(get(traceHandles{2}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    end
    
    if ~isempty(traceFigTimeAxis(cntState.Frames==cntState.CurrentFrame))
        traceHandles{3} = plot(traceFigAxes,traceFigTimeAxis(cntState.Frames==cntState.CurrentFrame),...
            ones(1),'ok');
        set(get(get(traceHandles{3}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    else
        traceHandles{3} = plot(traceFigAxes, [0, 1], [0, 1], 'ok');
        set(get(get(traceHandles{3}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
        set(traceHandles{3},'Visible','off')
    end
    
    
    
    xlabel(traceFigAxes,'frame')
    title(traceFigAxes, '', 'Interpreter', 'none');
    
    nc = cntState.schnitzcells(cntState.CurrentNucleus).cycle;
    
    
    traceFigYLimits = [0, 2];
    for i = 1:length(ncFrames)
        hndl_idx = length(traceHandles)+1;
        currentAnaphaseBoundary = ncFrames(i);
        traceHandles{hndl_idx} =plot(traceFigAxes,ones(1,2).*double(currentAnaphaseBoundary),traceFigYLimits,...
            'LineWidth',2,'Color','black', 'Marker', 'none');
        set(get(get(traceHandles{hndl_idx}, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    end
    try
        setPlotsInvisible(traceFigAxes);
        xlim(traceFigAxes,double([min([priorAnaphase, min(traceFigTimeAxis)]),max([max(traceFigTimeAxis), nextAnaphase])])+[-1,1]);
        setPlotsVisible(traceFigAxes);
    end
    
    
    try
        setPlotsInvisible(traceFigAxes);
        ylim(traceFigAxes, [0, traceFigYLimits(2)]);
        setPlotsVisible(traceFigAxes);
    end
    traceFigAxes.YAxis.Visible = 'off';
    hold(traceFigAxes,'off')
end
%
% idata1 = amp1(approvedParticleFrames);
% idata2 = amp2(approvedParticleFrames);
%
%

% creating axis title
frame_idx = find(cntState.Frames == cntState.CurrentFrame);
numNuclei = length(cntState.schnitzcells);
firstLine = [Prefix];
% STOPPED HERE
secondLine = ['Nucleus: ',num2str(cntState.CurrentNucleus),'/',num2str(numNuclei), ',   Frame: ',num2str(cntState.CurrentFrame),'/',num2str(numFrames),'    ',num2str(round(ElapsedTime(cntState.CurrentFrame))), 'm'];
if hasInputChannels
    thirdLine = ['Z: ',num2str(cntState.MidMedZ(frame_idx)),'/',num2str(cntState.ZSlices)];
end

if isfield(cntState.FrameInfo, 'nc')
    if hasInputChannels
        axisTitle={firstLine,...
            [secondLine,'    (nc',num2str(cntState.FrameInfo(cntState.CurrentFrame).nc),')'],...
            thirdLine};
    else
        axisTitle={firstLine,...
            [secondLine,'    (nc',num2str(cntState.FrameInfo(cntState.CurrentFrame).nc),')']};
    end
else
    if hasInputChannels
        axisTitle={firstLine,secondLine,thirdLine};
    else
        axisTitle={firstLine,secondLine};
    end
end

if cntState.HideApprovedFlag == 1
    axisTitle=[axisTitle,', Showing non-flagged particles'];
elseif cntState.HideApprovedFlag == 2
    axisTitle=[axisTitle,', Showing disapproved particles'];
end
hold(traceFigAxes, 'off');
setPlotsInvisible(traceFigAxes);
set(traceFigAxes.Title,'String', axisTitle);
setPlotsVisible(traceFigAxes);

set(traceFig,'Name',['Current Nucleus Trace'])

%% Figure 4
zProfileHandles = {};
zTraceHandles = {};
if hasInputChannels
    z = 1:cntState.ZSlices;
    fr_idx = find(cntState.Frames == cntState.CurrentFrame);
    zProfileFigAxes = subplot(1, 2, 1, 'Parent', zFig);
    
    if ~isempty(fr_idx)
        ZProfile = cntState.schnitzcells(cntState.CurrentNucleus).Fluo(fr_idx,:);
        
        zProfileHandles{1} = plot(zProfileFigAxes, z, ZProfile, 'r.-');
        hold(zProfileFigAxes, 'on')
        zProfileHandles{2} = plot(zProfileFigAxes,z(cntState.MidMedZ(fr_idx)),...
            ZProfile(cntState.MidMedZ(fr_idx)), 'ko');
        zProfileHandles{3} = plot(zProfileFigAxes,z(cntState.MaxZ(fr_idx)),...
            ZProfile(cntState.MaxZ(fr_idx)), 'bo');
        zProfileHandles{4} = plot(zProfileFigAxes,z(cntState.MedZ(fr_idx)),...
            ZProfile(cntState.MedZ(fr_idx)), 'go');
    else
        zProfileHandles{1} = plot(zProfileFigAxes, [0, 1], [0, 1], 'r.');
        hold(zProfileFigAxes, 'on')
        set(zProfileHandles{1},'Visible','off')
        zProfileHandles{2} = plot(zProfileFigAxes, [0, 1], [0, 1], 'ko');
        set(zProfileHandles{2},'Visible','off')
        zProfileHandles{3} = plot(zProfileFigAxes, [0, 1], [0, 1], 'bo');
        set(zProfileHandles{3},'Visible','off')
        zProfileHandles{4} = plot(zProfileFigAxes, [0, 1], [0, 1], 'go');
        set(zProfileHandles{4},'Visible','off')
    end
    
    
    xlabel(zProfileFigAxes,'z-slice', 'FontSize',12)
    ylabel(zProfileFigAxes,'intensity (au)', 'FontSize',12)
    title(zProfileFigAxes, 'z-Profile', 'Interpreter', 'none');
    zProfileFigAxes.Title.FontSize = 10;
    if ~isempty(fr_idx)
        zProfileFigYLimits = [0, max(ZProfile)*1.2];
    end
    if exist('zProfileFigYLimits', 'var')
        setPlotsInvisible(zProfileFigAxes);
        ylim(zProfileFigAxes, [0, max([zProfileFigYLimits(2), 1])]);
        setPlotsVisible(zProfileFigAxes);
    end
    
    
    zTraceAxes = subplot(1, 2, 2, 'Parent', zFig);
    
    zTraceHandles{1} = plot(zTraceAxes,cntState.Frames,...
        cntState.MidMedZ, 'k.-');
    hold(zTraceAxes,'on')
    zTraceHandles{2} = plot(zTraceAxes, cntState.Frames,...
        cntState.MaxZ,'.-b');
    zTraceHandles{3} = plot(zTraceAxes, cntState.Frames,...
        cntState.MedZ,'.-g');
    if ~isempty(fr_idx)
        zTraceHandles{4} = xline(zTraceAxes, cntState.Frames(fr_idx),...
            'Color', 'r', 'LineStyle', '--');
    else
        zTraceHandles{4} = xline(zTraceAxes, 0, 'Color', 'r', 'LineStyle', '--');
        set(zTraceHandles{4},'Visible','off')
    end
    
    hold(zTraceAxes, 'off')
    xlabel(zTraceAxes,'frame')
    ylabel(zTraceAxes,'z-slice')
    title(zTraceAxes,'Z traces')
    zTraceAxes.Title.FontSize = 10;
    zTraceLeg = legend(zTraceAxes,[zTraceHandles{1}, zTraceHandles{2}, zTraceHandles{3}], str1,str2,str3, 'AutoUpdate', 'off', 'HandleVisibility', 'off');
    
    zTraceFigYLimits = [0, cntState.ZSlices*1.3];
    try
        setPlotsInvisible(zTraceAxes);
        ylim(zTraceAxes, [0, zTraceFigYLimits(2)]);
        setPlotsVisible(zTraceAxes);
    end
    zFigAxes = {zProfileFigAxes, zTraceAxes};
    
    set(zFig,'Name',['Current Nucleus Z Profile'])
else
    zFigAxes = [];
end
%%

% Define the size of the figures/subplots
set(OverlayFig, 'units', 'normalized', 'position', [0.01, 0.45, .4, .4]);
if hasInputChannels
    set(snipFig, 'units', 'normalized', 'position', [0.01, 0.1, .4, .2]);
else
    set(snipFig, 'units', 'normalized', 'position', [0.01, 0.1, .2, .2]);
end
set(traceFig, 'units', 'normalized', 'position', [0.5, 0.5, .33, .33]);
if hasInputChannels
    set(zFig, 'units', 'normalized', 'position', [0.5, 0.05, .45, .33]);
end



end

