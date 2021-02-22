% 9/12 To do list:
% - write prefix of the data set of the particle in the figure title

close all;
dataset = 'mcp_opt';

colors = [0.1961    0.8039    0.1961;... % (1) lime green
    0.2941         0    0.7059;... % (2) indigo
    0.8353    0.4235    0.3333;... % (3) red (1st of PBoC colors)
    0.9176    0.7608    0.3922;... % (4) yellow
    0.4235    0.7373    0.9137;... % (5) cyan
    0.8157    0.4275    0.6706;... % (6) magenta
    0.4510    0.5569    0.7569;]; % (7) light blue

colorTraceOneSlice = colors(1,:);
colorTraceThreeSlices = colors(2,:);
colorsForMultiple = colors(3:7,:);
transparencyFitted = 0.55;


markerSize = 25;
lineWidth = 3;
fittedFactor = 2/3;
lineWidthFitted = lineWidth*fittedFactor;
waitTime = 0.1;
figTitleFontSize  = 20;
axisFontSize = 14;

%% Loading the data ------------------------------------------------------------------------
currentParticle = 1;
currentChannel = 1;

%---- section for running with many data sets ----
[allParticles, particlesCorrespondingDataSetNumber,...
    allSpots,spotsCorrespondingDataSetNumber,...
    Prefixes, otherRelevantData] = concatenateDataSets(dataset);

numberOfParticles = size(allParticles,2);
%---- section for running with many data sets ----


%% making trace figure ------------------------------------------------------------------------
particleTraceFig = figure();
positionOfParticleTraceFig = [4 49 996 635];%[4 97 996 587];%[403 49 562 635];
set(gcf,'Position',positionOfParticleTraceFig);

% Just going to draw one particle instead of drawing all of them
drawCurrentTrace(currentParticle,allParticles,particlesCorrespondingDataSetNumber,...
    allSpots,spotsCorrespondingDataSetNumber,currentChannel,... % particle info
    otherRelevantData,... % time info
    markerSize,lineWidth,colorTraceOneSlice,colorTraceThreeSlices,figTitleFontSize,axisFontSize) % plot info

% making slider to change particle
particleSliderPosition = [0.9268    0.2000    0.0315    0.7000]; % Units : normalized
% textPosition = [0.91 0.9 0.04 0.9100]; % Units : normalized
particleSliderSmallStep = 1/(numberOfParticles-1); % clicking will change it by 1
particleSliderBigStep = particleSliderSmallStep*4; % sliding will change it by units of 4
% Create slider
particleSlider = uicontrol('Style', 'slider',...
    'Min',1,'Max',numberOfParticles,'Value',currentParticle,...
    'Units', 'normalized','Tag','particleSlider',...
    'Position', particleSliderPosition,...
    'SliderStep',[particleSliderSmallStep particleSliderBigStep],...
    'Callback', {@chooseParticleAction,...
    allParticles,particlesCorrespondingDataSetNumber,...
    allSpots,spotsCorrespondingDataSetNumber,currentChannel,... % particle info
    otherRelevantData,... % time info
    markerSize,lineWidth,colorTraceOneSlice,colorTraceThreeSlices,figTitleFontSize,axisFontSize}); % plot info

traceFigHandles = guihandles(particleTraceFig);


%% making histogram figure  ------------------------------------------------------------------------
% this figure will show the histogram of interest
ncPresentInMovie = unique(correspondingNCInfo);
timeOnHistogramAllFigure = figure();
timeOnHistogramAPFigure = figure();
firstMeasuredHistogramAllFigure = figure();
firstMeasuredHistogramAPFigure = figure();
overlayTimeOnFirstMeasured = figure();
initialSlopeHistogramAllFigure = figure();
initialSlopeHistogramAPFigure = figure();

allHistogramHandles = [timeOnHistogramAllFigure, timeOnHistogramAPFigure,...
    firstMeasuredHistogramAllFigure, firstMeasuredHistogramAPFigure,...
    overlayTimeOnFirstMeasured,...
    initialSlopeHistogramAllFigure initialSlopeHistogramAPFigure];
positionOfHistogramFig = [4 49 996 635];%[4 97 996 587];

for currentHistogramFigureIndex = 1:length(allHistogramHandles)
    currentFigure = allHistogramHandles(currentHistogramFigureIndex);
    figure(currentFigure)
    set(gcf,'Position',positionOfHistogramFig);
    
    drawHistograms(currentFigure,... % figure info
        allParticles,particlesCorrespondingDataSetNumber,... % Particle info
        allSpots,spotsCorrespondingDataSetNumber,currentChannel,... % Spots info
        ncPresentInMovie,otherRelevantData,... % Time info
        colorTraceOneSlice,colorsForMultiple,... % plot color info
        figTitleFontSize,axisFontSize) % plot info
    
end

%% making sorted trace figure ------------------------------------------------------------------
sortedTraceFig = figure();
positionOfSortedTraceFig = [4 49 996 635];%[4 97 996 587];%[403 49 562 635];
set(gcf,'Position',positionOfSortedTraceFig);

% Just going to draw one particle instead of drawing all of them
drawCurrentSortedTrace(currentParticle,allParticles,particlesCorrespondingDataSetNumber,...
    allSpots,spotsCorrespondingDataSetNumber,currentChannel,... % particle info
    otherRelevantData,... % time info
    markerSize,lineWidth,colorTraceOneSlice,colorsForMultiple,figTitleFontSize,axisFontSize) % plot info

% making slider to change particle
sortedParticleSliderPosition = [0.9268    0.2000    0.0315    0.7000]; % Units : normalized
% textPosition = [0.91 0.9 0.04 0.9100]; % Units : normalized
particleSliderSmallStep = 1/(numberOfParticles-1); % clicking will change it by 1
particleSliderBigStep = particleSliderSmallStep*4; % sliding will change it by units of 4
% Create slider
sortedParticleSlider = uicontrol('Style', 'slider',...
    'Min',1,'Max',numberOfParticles,'Value',currentParticle,...
    'Units', 'normalized','Tag','particleSlider',...
    'Position', sortedParticleSliderPosition,...
    'SliderStep',[particleSliderSmallStep particleSliderBigStep],...
    'Callback', {@chooseSortedTraceAction,...
    allParticles,particlesCorrespondingDataSetNumber,...
    allSpots,spotsCorrespondingDataSetNumber,currentChannel,... % particle info
    otherRelevantData,... % time info
    markerSize,lineWidth,colorTraceOneSlice,colorsForMultiple,figTitleFontSize,axisFontSize}); % plot info

sortedTraceFigHandles = guihandles(sortedTraceFig);

%% making option figure ------------------------------------------------------------------------
optionsFig = figure();
positionOfOptionFig = [1017 48 255 635];%[1017 97 255 586];
maxButtonYPoint = 0.8318;
firstXPoint = 0.0745;
fullButtonWidth = 0.7;
buttonHeight = 0.0341;
tabFactor = 1.1;
checkBoxOneIndentFactor = 3;
sectionSpacing = buttonHeight*2;
currentSection = 0;
buttonSpacing = buttonHeight;
buttonCounter = 0; % this will increase by 1 before each check box
defaultOption = 0;

%% Single Trace Controls -----------------------------
% Add a text for the Title
textHeight = 0.0341;
textLength = 0.9;
textPosition = [firstXPoint maxButtonYPoint+2*textHeight...
    textLength textHeight*2];
titleString = 'Control Panel';
titleFontSize = 18;
titleTxt = uicontrol('Style','text',...
    'Units', 'normalized',...
    'Position',textPosition,...
    'FontSize',titleFontSize,...
    'String',titleString);

% Add a text for the section title
textHeight = 0.0341;
textLength = 0.9;
textPosition = [firstXPoint maxButtonYPoint+textHeight ...
    textLength textHeight*1.2];
sectionFontSize = 12;
sectionString = 'Metric Graphing Options';
sectionTitleTxt = uicontrol('Style','text',...
    'Units', 'normalized',...
    'Position',textPosition,...
    'FontSize',sectionFontSize,...
    'String',sectionString);

% Fitted initial slope button
buttonCounter = buttonCounter + 1;
set(gcf,'Position',positionOfOptionFig);
initialSlopeOption = uicontrol('Style','checkbox',...
    'String',{'Initial Slope'},'Min',defaultOption,'Max',buttonCounter,...
    'Value',defaultOption,'Tag','initialSlopeButton',...
    'Callback', {@initialSlopeOptionAction,particleTraceFig,traceFigHandles,...
    allParticles,particlesCorrespondingDataSetNumber,...
    allSpots,spotsCorrespondingDataSetNumber,currentChannel,... % particle info
    otherRelevantData,... % time info
    markerSize,lineWidth,lineWidthFitted,colorTraceOneSlice,colorTraceThreeSlices,...
    transparencyFitted,figTitleFontSize,axisFontSize}); % plot info
initialSlopeOption.Units = 'normalized';
currentButtonYPoint = maxButtonYPoint- buttonSpacing*(buttonCounter-1)...
    - currentSection*sectionSpacing;
currentButtonXPoint = firstXPoint;
initialSlopeOption.Position = [currentButtonXPoint currentButtonYPoint...
    fullButtonWidth buttonHeight];

% toggling evenly spacing with inital slope option
% this checkbox is indented
buttonCounter = buttonCounter + 1;
toggleRespacingOption = uicontrol('Style','checkbox',...
    'String',{'Center in X'},'Min',defaultOption,'Max',buttonCounter,...
    'Value',defaultOption,'Tag','toggleRespacingButton',...
    'Callback', {@toggleResizingOptionAction,particleTraceFig,traceFigHandles,...
    allParticles,particlesCorrespondingDataSetNumber,...
    allSpots,spotsCorrespondingDataSetNumber,currentChannel,... % particle info
    otherRelevantData,... % time info
    markerSize,lineWidth,lineWidthFitted,colorTraceOneSlice,colorTraceThreeSlices,transparencyFitted...
    ,figTitleFontSize,axisFontSize}); % plot info
toggleRespacingOption.Units = 'normalized';
currentButtonYPoint = maxButtonYPoint- buttonSpacing*(buttonCounter-1)...
    - currentSection*sectionSpacing;
% indenting this checkbox
currentButtonXPoint = firstXPoint + buttonSpacing*checkBoxOneIndentFactor;
toggleRespacingOption.Position = [currentButtonXPoint currentButtonYPoint...
    fullButtonWidth buttonHeight];


% looking at particle traces sorted by intensity....
buttonCounter = buttonCounter + 1;
toggleRespacingOption = uicontrol('Style','checkbox',...
    'String',{'Sort by Intensity'},'Min',defaultOption,'Max',buttonCounter,...
    'Value',defaultOption,'Tag','toggleRespacingButton',...
    'Callback', {@showSortedTraceFig,particleTraceFig,traceFigHandles,...
    sortedTraceFig, sortedTraceFigHandles}); % plot info
toggleRespacingOption.Units = 'normalized';
currentButtonYPoint = maxButtonYPoint- buttonSpacing*(buttonCounter-1)...
    - currentSection*sectionSpacing;
currentButtonXPoint = firstXPoint;
toggleRespacingOption.Position = [currentButtonXPoint currentButtonYPoint...
    fullButtonWidth buttonHeight];




%% Histogram controls -----------------------------
currentSection = currentSection + 1;
histogramCounter = 0;


% Add a text for the section title
textHeight = 0.0341;
textLength = 0.9;
currentYPoint = maxButtonYPoint- buttonSpacing*(buttonCounter-1)...
    - currentSection*sectionSpacing;
textPosition = [firstXPoint currentYPoint ...
    textLength textHeight*1.1];
sectionTitleString = 'Options';
sectionFontSize = 12;
sectionString = 'Histogram Options';
sectionTitleTxt = uicontrol('Style','text',...
    'Units', 'normalized',...
    'Position',textPosition,...
    'FontSize',sectionFontSize,...
    'String',sectionString);


% Initial Time Histograms ------------------------
buttonCounter = buttonCounter + 1;
currentYPoint = currentYPoint - buttonSpacing;
textPosition = [firstXPoint currentYPoint ...
    textLength textHeight*1.1];
sectionFontSize = 8;
sectionString = 'Initial Time Histograms';
initialTimeSubsection = uicontrol('Style','text',...
    'Units', 'normalized',...
    'Position',textPosition,...
    'FontSize',sectionFontSize,...
    'String',sectionString,...
    'HorizontalAlignment','left');

% all nc
buttonCounter = buttonCounter + 1;
histogramCounter = histogramCounter + 1;
initialTimeHistogram = uicontrol('Style','checkbox',...
    'String','All Nuclear Cycles & AP','Min',defaultOption,'Max',histogramCounter,...
    'Value',defaultOption,'Tag','All Nuclear Cycles & AP',...
    'Callback', {@histogramOptionAction,particleTraceFig});
initialTimeHistogram.Units = 'normalized';
currentYPoint = currentYPoint - buttonSpacing./2;
% indenting this checkbox
currentButtonXPoint = firstXPoint + buttonSpacing*checkBoxOneIndentFactor;
initialTimeHistogram.Position = [currentButtonXPoint currentYPoint...
    fullButtonWidth buttonHeight];

% by ap position
buttonCounter = buttonCounter + 1;
histogramCounter = histogramCounter + 1;
initialMeasuredHistogram = uicontrol('Style','checkbox',...
    'String','% AP Position',...
    'Min',defaultOption,'Max',histogramCounter,...
    'Value',defaultOption,'Tag','% AP Position',...
    'Callback', {@histogramOptionAction,particleTraceFig});
initialMeasuredHistogram.Units = 'normalized';
currentYPoint = currentYPoint - buttonSpacing;
% indenting this checkbox
currentButtonXPoint = firstXPoint + buttonSpacing*checkBoxOneIndentFactor;
initialMeasuredHistogram.Position = [currentButtonXPoint currentYPoint...
    fullButtonWidth buttonHeight];

% First measured point ------------------
currentYPoint = currentYPoint - buttonSpacing*2;
textPosition = [firstXPoint currentYPoint ...
    textLength textHeight*1.1];
sectionFontSize = 8;
sectionString = 'First Measured Time Point Histograms';
firstMeasuredPointSubsection = uicontrol('Style','text',...
    'Units', 'normalized',...
    'Position',textPosition,...
    'FontSize',sectionFontSize,...
    'String',sectionString,...
    'HorizontalAlignment','left');

% all nc and ap
buttonCounter = buttonCounter + 1;
histogramCounter = histogramCounter + 1;
initialMeasuredHistogram = uicontrol('Style','checkbox',...
    'String','All Nuclear Cycles & AP',...
    'Min',defaultOption,'Max',histogramCounter,...
    'Value',defaultOption,'Tag','All Nuclear Cycles & AP',...
    'Callback', {@histogramOptionAction,particleTraceFig});
initialMeasuredHistogram.Units = 'normalized';
currentYPoint = currentYPoint - buttonSpacing./2;
% indenting this checkbox
currentButtonXPoint = firstXPoint + buttonSpacing*checkBoxOneIndentFactor;
initialMeasuredHistogram.Position = [currentButtonXPoint currentYPoint...
    fullButtonWidth buttonHeight];

% by ap
buttonCounter = buttonCounter + 1;
histogramCounter = histogramCounter + 1;
initialMeasuredHistogram = uicontrol('Style','checkbox',...
    'String','% AP',...
    'Min',defaultOption,'Max',histogramCounter,...
    'Value',defaultOption,'Tag','% AP',...
    'Callback', {@histogramOptionAction,particleTraceFig});
initialMeasuredHistogram.Units = 'normalized';
currentYPoint = currentYPoint - buttonSpacing;
% indenting this checkbox
currentButtonXPoint = firstXPoint + buttonSpacing*checkBoxOneIndentFactor;
initialMeasuredHistogram.Position = [currentButtonXPoint currentYPoint...
    fullButtonWidth buttonHeight];


% First measured and time on!!!! ------------------
currentYPoint = currentYPoint - buttonSpacing*2;
textPosition = [firstXPoint currentYPoint ...
    textLength textHeight*1.1];
sectionFontSize = 8;
sectionString = 'First Measured & Time On Overlay Histograms';
firstMeasuredAndTimeOnSubsection = uicontrol('Style','text',...
    'Units', 'normalized',...
    'Position',textPosition,...
    'FontSize',sectionFontSize,...
    'String',sectionString,...
    'HorizontalAlignment','left');

% overlay
buttonCounter = buttonCounter + 1;
histogramCounter = histogramCounter + 1;
initialMeasuredHistogram = uicontrol('Style','checkbox',...
    'String','Overlay',...
    'Min',defaultOption,'Max',histogramCounter,...
    'Value',defaultOption,'Tag','overlay',...
    'Callback', {@histogramOptionAction,particleTraceFig});
initialMeasuredHistogram.Units = 'normalized';
currentYPoint = currentYPoint - buttonSpacing./2;
% indenting this checkbox
currentButtonXPoint = firstXPoint + buttonSpacing*checkBoxOneIndentFactor;
initialMeasuredHistogram.Position = [currentButtonXPoint currentYPoint...
    fullButtonWidth buttonHeight];


% Initial Slope ------------------
currentYPoint = currentYPoint - buttonSpacing*2;
textPosition = [firstXPoint currentYPoint ...
    textLength textHeight*1.1];
sectionFontSize = 8;
sectionString = 'Initial Slope Histograms';
firstMeasuredPointSubsection = uicontrol('Style','text',...
    'Units', 'normalized',...
    'Position',textPosition,...
    'FontSize',sectionFontSize,...
    'String',sectionString,...
    'HorizontalAlignment','left');

% all nc and ap
buttonCounter = buttonCounter + 1;
histogramCounter = histogramCounter + 1;
initialSlopeHistogram = uicontrol('Style','checkbox',...
    'String','All Nuclear Cycles & AP',...
    'Min',defaultOption,'Max',histogramCounter,...
    'Value',defaultOption,'Tag','Initial Slope',...
    'Callback', {@histogramOptionAction,particleTraceFig});
initialSlopeHistogram.Units = 'normalized';
currentYPoint = currentYPoint - buttonSpacing./2;
% indenting this checkbox
currentButtonXPoint = firstXPoint + buttonSpacing*checkBoxOneIndentFactor;
initialSlopeHistogram.Position = [currentButtonXPoint currentYPoint...
    fullButtonWidth buttonHeight];

% vs % ap
buttonCounter = buttonCounter + 1;
histogramCounter = histogramCounter + 1;
initialSlopeHistogram = uicontrol('Style','checkbox',...
    'String','% AP Position',...
    'Min',defaultOption,'Max',histogramCounter,...
    'Value',defaultOption,'Tag','Initial Slope',...
    'Callback', {@histogramOptionAction,particleTraceFig});
initialSlopeHistogram.Units = 'normalized';
currentYPoint = currentYPoint - buttonSpacing;
% indenting this checkbox
currentButtonXPoint = firstXPoint + buttonSpacing*checkBoxOneIndentFactor;
initialSlopeHistogram.Position = [currentButtonXPoint currentYPoint...
    fullButtonWidth buttonHeight];


% % template button
% buttonCounter = buttonCounter + 1;
% select = uicontrol('Style','checkbox',...
%     'String',{'Testing'},'Min',defaultOption,'Max',buttonCounter,...
%     'Value',defaultOption,'Tag','selectButton',...
%     'Callback', {@optionAction,particleTraceFig,traceFigHandles,...
%     Particles,Spots,currentChannel,...
%     markerSize,lineWidth,lineWidthFitted,color1,color3,transparencyFitted,...
%     figTitleFontSize, axisFontSize});
% select.Units = 'normalized';
% currentButtonYPoint = maxButtonYPoint- buttonSpacing*(buttonCounter-1)...
%     - currentSection*sectionSpacing;
% currentButtonXPoint = firstXPoint;
% select.Position = [currentButtonXPoint currentButtonYPoint...
%     fullButtonWidth buttonHeight];





% have the first graph shown be the trace figure
figure(particleTraceFig);

% maybe have options for choosing the time scale? seconds vs minutes?

%% other functions (most will be placed in a new script)
% redrawing the trace of the particle that was choosen
function chooseParticleAction(source,event,allParticles,particlesCorrespondingDataSetNumber,...
    allSpots,spotsCorrespondingDataSetNumber,currentChannel,... % particle info
    otherRelevantData,... % time info
    markerSize,lineWidth,colorTraceOneSlice,colorTraceThreeSlices,...
    figTitleFontSize,axisFontSize) % plot info

currentParticle = round(source.Value);

drawCurrentTrace(currentParticle,allParticles,particlesCorrespondingDataSetNumber,...
    allSpots,spotsCorrespondingDataSetNumber,currentChannel,... % particle info
    otherRelevantData,... % time info
    markerSize,lineWidth,colorTraceOneSlice,colorTraceThreeSlices,...
    figTitleFontSize,axisFontSize) % plot info

end

function chooseSortedTraceAction(source,event,...
    allParticles,particlesCorrespondingDataSetNumber,...
    allSpots,spotsCorrespondingDataSetNumber,currentChannel,... % particle info
    otherRelevantData,... % time info
    markerSize,lineWidth,colorTraceOneSlice,colorTraceThreeSlices,...
    figTitleFontSize,axisFontSize)

smallFigurePosition = [1022 47 254 132];
smallFigureNumber = 11;

currentParticle = round(source.Value);

drawCurrentSortedTrace(currentParticle,allParticles,particlesCorrespondingDataSetNumber,...
    allSpots,spotsCorrespondingDataSetNumber,currentChannel,... % particle info
    otherRelevantData,... % time info
    markerSize,lineWidth,colorTraceOneSlice,colorTraceThreeSlices,...
    figTitleFontSize,axisFontSize) % plot info

if ishandle(smallFigureNumber)
    figure(smallFigureNumber)
else
    smallFigureHandle = figure();
    set(smallFigureHandle,'Position',smallFigurePosition)
end

drawCurrentTrace(currentParticle,allParticles,particlesCorrespondingDataSetNumber,...
    allSpots,spotsCorrespondingDataSetNumber,currentChannel,... % particle info
    otherRelevantData,... % time info
    markerSize,lineWidth,colorTraceOneSlice,colorTraceThreeSlices,...
    figTitleFontSize,axisFontSize,'smallFigure') % plot info
end

function initialSlopeOptionAction(source,event,particleTraceFig,traceFigHandles,...
    Particles,particlesCorrespondingDataSetNumber,...
    Spots,spotsCorrespondingDataSetNumber,currentChannel,... % particle info
    otherRelevantData,... % time info
    markerSize,lineWidth,lineWidthFitted,color1,color3,transparencyFitted...
    ,figTitleFontSize,axisFontSize)

figure(particleTraceFig);
rescaling = 0;

%find the current particle
currentParticle = round(traceFigHandles.particleSlider.Value);

%find the values of the other check boxes

drawCurrentTrace(currentParticle,Particles,particlesCorrespondingDataSetNumber,...
    Spots,spotsCorrespondingDataSetNumber,currentChannel,... % particle info
    otherRelevantData,... % time info
    markerSize,lineWidth,color1,color3,figTitleFontSize,axisFontSize) % plot info


if source.Value > 0 % resets figure if source.Value < 0
    plotFittedLine(currentParticle,Particles,particlesCorrespondingDataSetNumber,...
        Spots,spotsCorrespondingDataSetNumber,currentChannel,... % particle info
        otherRelevantData,... % time info
        lineWidthFitted,color1,color3,transparencyFitted,rescaling)
end

end

function toggleResizingOptionAction(source,event,particleTraceFig,traceFigHandles,...
    Particles,particlesCorrespondingDataSetNumber,...
    Spots,spotsCorrespondingDataSetNumber,currentChannel,... % particle info
    otherRelevantData,... % time info
    markerSize,lineWidth,lineWidthFitted,color1,color3,transparencyFitted,...
    figTitleFontSize,axisFontSize)

figure(particleTraceFig);

%find the current particle
currentParticle = round(traceFigHandles.particleSlider.Value);

drawCurrentTrace(currentParticle,Particles,particlesCorrespondingDataSetNumber,...
    Spots,spotsCorrespondingDataSetNumber,currentChannel,... % particle info
    otherRelevantData,... % time info
    markerSize,lineWidth,color1,color3,figTitleFontSize,axisFontSize) % plot info

if source.Value > 0
    rescaling = 1;
else
    rescaling = 0;
end
plotFittedLine(currentParticle,Particles,particlesCorrespondingDataSetNumber,...
    Spots,spotsCorrespondingDataSetNumber,currentChannel,... % particle info
    otherRelevantData,... % time info
    lineWidthFitted,color1,color3,transparencyFitted,rescaling)

end

function histogramOptionAction(source,event,particleTraceFig)
figureShift = particleTraceFig.Number;
if source.Value > 0
    figure(source.Value + figureShift)
else
    figure(particleTraceFig)
end

end


function showSortedTraceFig(source,event,particleTraceFig,traceFigHandles,...
    sortedTraceFig, sortedTraceFigHandles)
subFigureNumber = 11;

if source.Value > 0
%     figure(subFigureNumber)
    figure(sortedTraceFig)
else
    if ishandle(subFigureNumber)
        close(subFigureNumber);
    end
    figure(particleTraceFig)
end

end

%% Calculation functions
% calculate the initial slope of the current particle
% draw the current particle's trace


function drawCurrentSortedTrace(currentParticle,allParticles,particlesCorrespondingDataSetNumber,...
    allSpots,spotsCorrespondingDataSetNumber,currentChannel,... % particle info
    otherRelevantData,... % time info
    markerSize,lineWidth,colorTraceOneSlice,colorTraceThreeSlices,...
    figTitleFontSize,axisFontSize) % plot info


currentDataSet = particlesCorrespondingDataSetNumber(currentParticle);
correspondingNCInfo = otherRelevantData(currentDataSet).correspondingNCInfo;
framesInMinutes = otherRelevantData(currentDataSet).framesInMinutes;
nuclearCycleBoundariesInMinutes = otherRelevantData(currentDataSet).nuclearCycleBoundariesInMinutes;


[currentParticleSubset,ParticlesSubset,SpotsSubset] = ...
    getSubsetData(currentParticle,allParticles,particlesCorrespondingDataSetNumber,...
    allSpots,spotsCorrespondingDataSetNumber);

% getting trace of current particle
[frame,ampIntegral,ampIntegral3,~,~,~,errorIntegral,...
    ~,~,~,errorIntegral3, ~,~]=GetParticleTrace(currentParticleSubset,...
    ParticlesSubset{currentChannel},SpotsSubset{currentChannel});
traceLength = length(frame);

% adjusting the time where 0 is the first time point of the nc of the trace
frameRange = framesInMinutes(frame); % Units to seconds
ncPresent = unique(correspondingNCInfo(frame));
timeOfFirstNC = nuclearCycleBoundariesInMinutes(ncPresent(1)-8);
% adjusting frameRange
frameRange = frameRange - timeOfFirstNC;

% finding the max
[~, correspondingMaxIndex] = max(ampIntegral3);
indexRange = {1:correspondingMaxIndex; correspondingMaxIndex+1:traceLength};

% plotting the traces
orderedXTickMarks = 1:traceLength;
cla
hold on

% sorting
[sortedIntensities3, correspondingIndex3] = sort(ampIntegral3);
% up to max potion 
tempIndex = find(correspondingIndex3<=correspondingMaxIndex);
% upToMax = 1:currentMaxIndex;
currentReorderedSubsetIndex = tempIndex;
currentReorderedSubsetAmp = sortedIntensities3(tempIndex);
subTraceLength = length(tempIndex);

errorbar(currentReorderedSubsetIndex,currentReorderedSubsetAmp,...
    ones(1,subTraceLength)*errorIntegral3,...
    '.','MarkerSize',markerSize,'Color',colorTraceThreeSlices(2,:),...
    'DisplayName','Up to max')

%after max potion
if ~isequal(correspondingMaxIndex,length(frameRange))
    hold on 
    tempIndex = find(correspondingIndex3>correspondingMaxIndex);
%     afterMax = currentMaxIndex+1:length(correspondingIndex3);
    currentReorderedSubsetIndex = tempIndex;
    currentReorderedSubsetAmp = sortedIntensities3(tempIndex);
    subTraceLength = length(tempIndex);
    
    errorbar(currentReorderedSubsetIndex,currentReorderedSubsetAmp,...
        ones(1,subTraceLength)*errorIntegral3,...
        '.','MarkerSize',markerSize,'Color',colorTraceThreeSlices(3,:),...
        'DisplayName','After max')
end

xticks(orderedXTickMarks)
correspondingSortedXValues3 = frameRange(correspondingIndex3);
roundedCorrespondingSortedXValues3 = round(correspondingSortedXValues3,2);
xticklabels(num2cell(roundedCorrespondingSortedXValues3))

legend('show'); % 'FontSize',axisFontSize (?)


ncRange = correspondingNCInfo(frame);
ncRange = unique(ncRange); % removing repeats

basicStringInfo = ['Particle ' num2str(currentParticleSubset)];
% put in the frame boundaries of the nc?
if ~isempty(ncRange)
    titleInformation = {basicStringInfo, ['nc = ' num2str(ncRange)]};
else
    error('There is nothing to plot!')
end

title(titleInformation,'FontSize',figTitleFontSize)
xlabel('Time (minutes)','FontSize',axisFontSize)
ylabel('Integrated Intensity','FontSize',axisFontSize)

xPadding = 0.3*mean(frameRange);
if length(frameRange) > 1
    xPadding = 0.3*(frameRange(end)-frameRange(1));
end
% minX = frameRange(1)-xPadding;
% maxX = frameRange(end)+xPadding;
minX = 0;
maxX = traceLength+1;
xlim([minX maxX])

% change y padding to include the error bars
yPaddingFactor = [0.2, 0.6];
upperPadding = 1 + yPaddingFactor(1);
upperBounds = [(ampIntegral+errorIntegral) (ampIntegral3+errorIntegral3)];

maxY = max(upperBounds)*upperPadding;

minY = 0;
ylim([minY maxY])
end




% plot the fitted line
function plotFittedLine(currentParticle,Particles,particlesCorrespondingDataSetNumber,...
    Spots,spotsCorrespondingDataSetNumber,currentChannel,... % particle info
    otherRelevantData,... % time info
    lineWidthFitted,color1,color3,transparencyFitted,rescaling)

currentDataSet = particlesCorrespondingDataSetNumber(currentParticle);
correspondingNCInfo = otherRelevantData(currentDataSet).correspondingNCInfo;
framesInMinutes = otherRelevantData(currentDataSet).framesInMinutes;
nuclearCycleBoundariesInMinutes = otherRelevantData(currentDataSet).nuclearCycleBoundariesInMinutes;


[currentParticleSubset,ParticlesSubset,SpotsSubset] = ...
    getSubsetData(currentParticle,Particles,particlesCorrespondingDataSetNumber,...
    Spots,spotsCorrespondingDataSetNumber);

[frame,ampIntegral,ampIntegral3,~,~,~,~,~,~,~,~,~,~]=...
    GetParticleTrace(currentParticleSubset,...
    ParticlesSubset{currentChannel},SpotsSubset{currentChannel});

frameRange = framesInMinutes(frame); % Changing the units to seconds
ncPresent = unique(correspondingNCInfo(frame));
timeOfFirstNC = nuclearCycleBoundariesInMinutes(ncPresent(1)-8);
% adjusting frameRange
frameRange = frameRange - timeOfFirstNC;

% Initial Slopes
[fittedEQ,indexOfMax] = ...
    initialSlopeCalc(frameRange,ampIntegral); % one slice
[fittedEQ3,indexOfMax3] = ...
    initialSlopeCalc(frameRange,ampIntegral3); % three slice

% plotting on the same figure
hold on

% plotting fitted lines
if ~isnan(fittedEQ)
    xSlope = frameRange(1:indexOfMax); % add one more frame around edge?
    extrapolatedInitialTime = roots(fittedEQ);
    yFitted = polyval(fittedEQ,xSlope);
    
    xSlope3 = frameRange(1:indexOfMax3);
    extrapolatedInitialTime3 = roots(fittedEQ3);
    yFitted3 = polyval(fittedEQ3,xSlope3);
    
    % plot for 1 slice integration
    temp = plot([extrapolatedInitialTime xSlope],[0 yFitted],...
        'LineWidth',lineWidthFitted,...
        'Color',color1);
    temp.Color(4) = transparencyFitted;
    currentEquation = fittedEQ;
    
    firstPart = 'fit: ';
    m = num2str(currentEquation(1));
    if currentEquation(2) < 0
        addOn = ' - ';
    else
        addOn = ' + ';
    end
    b = [addOn num2str(abs(currentEquation(2)))];
    temp.DisplayName = [firstPart m 'x' b];
    
    % plot for 3 slice integration
    temp = plot([extrapolatedInitialTime3 xSlope3],[0 yFitted3],...
        'LineWidth',lineWidthFitted,...
        'Color',color3);
    temp.Color(4) = transparencyFitted;
    currentEquation = fittedEQ3;
    firstPart = 'fit3: ';
    m = num2str(currentEquation(1));
    if currentEquation(2) < 0
        addOn = ' - ';
    else
        addOn = ' + ';
    end
    b = [addOn num2str(abs(currentEquation(2)))];
    temp.DisplayName = [firstPart m 'x' b];
    
    
    furtherPoint = min(extrapolatedInitialTime,extrapolatedInitialTime3);
    if rescaling
        maxXLimit = frameRange(1)-furtherPoint + frameRange(end);
        xlim([furtherPoint maxXLimit])
    else
        xPadding = 0.3*mean(frameRange);
        if length(frameRange) > 1
            xPadding = 0.3*(frameRange(end)-frameRange(1));
        end
        xlim([furtherPoint-xPadding frameRange(end)+xPadding])
    end
    % ylim
    legend();
end


end

function [extrapolatedInitialTime, extrapolatedInitialTime3] = ...
    findInitialTime(Particles,particlesCorrespondingDataSetNumber,...
    Spots,spotsCorrespondingDataSetNumber,currentChannel,...
    numberOfParticles,otherRelevantData,...
    ncOfInterest,dataSetsOfInterest)
%  useAPBins = 0; % default is to use all particles not separate by apbins
%
% for i = 1:length(varargin)
%     switch varargin(i)
%         case 'useAPBins'
%             useAPBins = 1;
%
%     end
% end

% Making sure the variables are of type cell
if ~iscell(Particles)
    Particles = {Particles};
end
if ~iscell(Spots)
    Spots = {Spots};
end


extrapolatedInitialTime = [];
extrapolatedInitialTime3 = [];
counter = 0;

% only looking at the particles indicated in dataSetsOfInterest
particleIndexOfInterest = [];
for currentDataSet = dataSetsOfInterest
    particlesFound = find(particlesCorrespondingDataSetNumber == currentDataSet);
    particleIndexOfInterest = [particleIndexOfInterest particlesFound];
end

for currentParticle = particleIndexOfInterest
    currentDataSet = particlesCorrespondingDataSetNumber(currentParticle);
    correspondingNCInfo = otherRelevantData(currentDataSet).correspondingNCInfo;
    framesInMinutes = otherRelevantData(currentDataSet).framesInMinutes;
    nuclearCycleBoundariesInMinutes = otherRelevantData(currentDataSet).nuclearCycleBoundariesInMinutes;
    
    [currentParticleSubset,ParticlesSubset,SpotsSubset] = ...
        getSubsetData(currentParticle,Particles,particlesCorrespondingDataSetNumber,...
        Spots,spotsCorrespondingDataSetNumber);
    
    % insert subset into here
    [frame,ampIntegral,ampIntegral3,~,~,~,~,~,~,~,~,~,~]=...
        GetParticleTrace(currentParticleSubset,...
        ParticlesSubset{currentChannel},SpotsSubset{currentChannel});
    
    currentNCRange = correspondingNCInfo(frame);
    currentNCRange = unique(currentNCRange); % removing repeats
    if sum(ismember(currentNCRange,ncOfInterest))
        
        frameRange = framesInMinutes(frame); % Units to seconds
        ncPresent = unique(correspondingNCInfo(frame));
        timeOfFirstNC = nuclearCycleBoundariesInMinutes(ncPresent(1)-8);
        % adjusting frameRange
        frameRange = frameRange - timeOfFirstNC;
        
        % Initial Slopes
        [fittedEQ,~] = initialSlopeCalc(frameRange,ampIntegral); % one slice
        [fittedEQ3,~] = initialSlopeCalc(frameRange,ampIntegral3); % three slice
        
        if ~sum(isnan([fittedEQ fittedEQ3])) && ~sum(isinf([fittedEQ fittedEQ3]))
            counter = counter + 1;
            extrapolatedInitialTime(counter) = roots(fittedEQ);
            extrapolatedInitialTime3(counter) = roots(fittedEQ3);
        end
    end
end

end

function initialMeasured = ...
    findInitialMeasured(Particles,particlesCorrespondingDataSetNumber,...
    Spots,spotsCorrespondingDataSetNumber,currentChannel,...
    numberOfParticles,otherRelevantData,...
    ncOfInterest,dataSetsOfInterest)

% Making sure the variables are of type cell
if ~iscell(Particles)
    Particles = {Particles};
end
if ~iscell(Spots)
    Spots = {Spots};
end


initialMeasured = [];
counter = 0;

particleIndexOfInterest = [];
for currentDataSet = dataSetsOfInterest
    particlesFound = find(particlesCorrespondingDataSetNumber == currentDataSet);
    particleIndexOfInterest = [particleIndexOfInterest particlesFound];
end

for currentParticle = particleIndexOfInterest
    
    currentDataSet = particlesCorrespondingDataSetNumber(currentParticle);
    correspondingNCInfo = otherRelevantData(currentDataSet).correspondingNCInfo;
    framesInMinutes = otherRelevantData(currentDataSet).framesInMinutes;
    nuclearCycleBoundariesInMinutes = otherRelevantData(currentDataSet).nuclearCycleBoundariesInMinutes;
    
    [currentParticleSubset,ParticlesSubset,SpotsSubset] = ...
        getSubsetData(currentParticle,Particles,particlesCorrespondingDataSetNumber,...
        Spots,spotsCorrespondingDataSetNumber);
    
    [frame,~,~,~,~,~,~,~,~,~,~,~,~]=...
        GetParticleTrace(currentParticleSubset,...
        ParticlesSubset{currentChannel},SpotsSubset{currentChannel});
    
    currentNCRange = correspondingNCInfo(frame);
    currentNCRange = unique(currentNCRange); % removing repeats
    if sum(ismember(currentNCRange,ncOfInterest))
        
        frameRange = framesInMinutes(frame); % Units to seconds
        ncPresent = unique(correspondingNCInfo(frame));
        timeOfFirstNC = nuclearCycleBoundariesInMinutes(ncPresent(1)-8);
        % adjusting frameRange
        frameRange = frameRange - timeOfFirstNC;
        
        % first measured point
        counter = counter + 1;
        initialMeasured(counter) = frameRange(1);
    end
end

end


function drawHistograms(histogramFigure,...
    Particles,particlesCorrespondingDataSetNumber,...
    Spots,spotsCorrespondingDataSetNumber,currentChannel,ncOfInterest,...
    otherRelevantData,... % particle info
    color1,color3,figTitleFontSize,axisFontSize) % plot info

if iscell(Particles)
    numberOfParticles = size(Particles{currentChannel},2);
else
    numberOfParticles = size(Particles,2);
end

figure(histogramFigure)
cla
% do action
histogramHeights = [];
binLimits = [];
apBinWidth = 0.2; % 20% embryo length bins
makeSupTitle = 0; % make a super title if subplots were used

dataSetsOfInterest = [6 2];
switch histogramFigure.Number - 1 % there is a particle figure plotted first
    case 1 % Extrapolated Initial Time Histogram (for all ap)
        dataSetsOfInterest = unique(particlesCorrespondingDataSetNumber);
        
        [extrapolatedInitialTime, extrapolatedInitialTime3] = ...
            findInitialTime(Particles,particlesCorrespondingDataSetNumber,...
            Spots,spotsCorrespondingDataSetNumber,currentChannel,...
            numberOfParticles,otherRelevantData,ncOfInterest,dataSetsOfInterest);
        
        % Plotting histogram (1 slice)
        %         currentHistogram = histogram(extrapolatedInitialTime,...
        %             'FaceColor',color3(2,:),...
        %             'BinMethod','fd','DisplayName','One Slice');
        %         histogramHeights = [histogramHeights currentHistogram.Values];
        %         binLimits = [binLimits currentHistogram.BinLimits];
        %         hold on
        
        % Plotting histogram (3 slices)
        currentHistogram = histogram(extrapolatedInitialTime3,...
            'FaceColor',color3(3,:),...
            'BinMethod','fd','DisplayName','Three Slices');
        histogramHeights = [histogramHeights currentHistogram.Values];
        binLimits = [binLimits currentHistogram.BinLimits];
        
        
        % Setting X limits
        xlim([min(binLimits) max(binLimits)])
        
        % Setting Y limits and labels
        tallestBarHeight = max(histogramHeights);
        ylim([0 tallestBarHeight+1])
        desiredYRange = 0:tallestBarHeight+1;
        yAxisLabels = cellfun(@num2str, num2cell(desiredYRange),...
            'UniformOutput', false);
        set(gca, 'YTick',desiredYRange,'YTickLabel', yAxisLabels)
        xLabelText = 'Time (minutes)';
        quantity = 'Time on';
        legend('show')
        
    case 2 % Extrapolated Initial Time Histogram (with AP)
        
        quantity = 'Time on';
        colorCounter = 0;
        for i = dataSetsOfInterest
            [extrapolatedInitialTime, extrapolatedInitialTime3] = ...
                findInitialTime(Particles,particlesCorrespondingDataSetNumber,...
                Spots,spotsCorrespondingDataSetNumber,currentChannel,...
                numberOfParticles,otherRelevantData,ncOfInterest,i);
            
            % Plotting Histogram 1 slice
            %         subplot(1,2,1)
            %         currentHistogram = histogram(extrapolatedInitialTime,...
            %             'FaceColor',color1,...
            %             'BinMethod','fd','DisplayName','One Slice');
            %         histogramHeights = [histogramHeights currentHistogram.Values];
            %         binLimits = [currentHistogram.BinLimits];
            %         % Setting X limits
            %         xlim([min(binLimits) max(binLimits)])
            
            %         % Setting Y limits and labels
            %         tallestBarHeight = max(histogramHeights);
            %         ylim([0 tallestBarHeight+1])
            %         desiredYRange = 0:tallestBarHeight+1;
            %         yAxisLabels = cellfun(@num2str, num2cell(desiredYRange),...
            %             'UniformOutput', false);
            %         set(gca, 'YTick',desiredYRange,'YTickLabel', yAxisLabels)
            %         xlabel('Time (minutes)','FontSize',axisFontSize)
            %         ylabel('Counts','FontSize',axisFontSize)
            
            % Plotting Histogram 3 Slices
            %         subplot(1,2,2)
            colorCounter = colorCounter + 1;
            currentHistogram = histogram(extrapolatedInitialTime3,...
                'FaceColor',color3(colorCounter+1,:),...
                'BinMethod','fd','DisplayName',['Data set : ' num2str(i+9)]);
            hold on
            histogramHeights = [histogramHeights currentHistogram.Values];
            binLimits = [currentHistogram.BinLimits];
        end
        % Setting X limits
        xlim([min(binLimits) max(binLimits)])
        
        % Setting Y limits and labels
        tallestBarHeight = max(histogramHeights);
        ylim([0 tallestBarHeight+1])
        desiredYRange = 0:tallestBarHeight+1;
        yAxisLabels = cellfun(@num2str, num2cell(desiredYRange),...
            'UniformOutput', false);
        set(gca, 'YTick',desiredYRange,'YTickLabel', yAxisLabels)
        xLabelText = 'Time (minutes)';
        legend('show')
        %         makeSupTitle = 1;
        
    case 3 % Initial measured time (for all ap)
        dataSetsOfInterest = unique(particlesCorrespondingDataSetNumber);
        
        initialMeasured = findInitialMeasured(Particles,particlesCorrespondingDataSetNumber,...
            Spots,spotsCorrespondingDataSetNumber,currentChannel,...
            numberOfParticles,otherRelevantData,ncOfInterest,dataSetsOfInterest);
        
        % Plotting histogram
        currentHistogram = histogram(initialMeasured,...
            'FaceColor',color3(3,:),'BinMethod','fd');
        histogramHeights = [histogramHeights currentHistogram.Values];
        binLimits = [binLimits currentHistogram.BinLimits];
        
        % Setting X limits
        xlim([min(binLimits) max(binLimits)])
        
        % Setting Y limits and labels
        tallestBarHeight = max(histogramHeights);
        ylim([0 tallestBarHeight+1])
        desiredYRange = 0:tallestBarHeight+1;
        yAxisLabels = cellfun(@num2str, num2cell(desiredYRange),...
            'UniformOutput', false);
        set(gca, 'YTick',desiredYRange,'YTickLabel', yAxisLabels)
        quantity = 'First Data Point';
        xLabelText = 'Time (minutes)';
%         legend('show')
        
    case 4 % case for Initial time measured (with AP)
        colorCounter = 0;
        for i = dataSetsOfInterest
            initialMeasured = findInitialMeasured(Particles,particlesCorrespondingDataSetNumber,...
                Spots,spotsCorrespondingDataSetNumber,currentChannel,...
                numberOfParticles,otherRelevantData,ncOfInterest,i);
            colorCounter = colorCounter + 1;
            currentHistogram = histogram(initialMeasured,...
                'FaceColor',color3(colorCounter+1,:),'BinMethod','fd','DisplayName',['Data set : ' num2str(i+9)]);
            hold on
            histogramHeights = [histogramHeights currentHistogram.Values];
            binLimits = [binLimits currentHistogram.BinLimits];
        end
        % Setting X limits
        xlim([min(binLimits) max(binLimits)])
        
        % Setting Y limits and labels
        tallestBarHeight = max(histogramHeights);
        ylim([0 tallestBarHeight+1])
        desiredYRange = 0:tallestBarHeight+1;
        yAxisLabels = cellfun(@num2str, num2cell(desiredYRange),...
            'UniformOutput', false);
        set(gca, 'YTick',desiredYRange,'YTickLabel', yAxisLabels)
        quantity = 'First Data Point';
        xLabelText = 'Time (minutes)';
        legend('show')
    case 5 % overlay of histograms of initial measured time point and time on 
        dataSetsOfInterest = unique(particlesCorrespondingDataSetNumber);
        
        % --------- time on ----------
        [extrapolatedInitialTime, extrapolatedInitialTime3] = ...
            findInitialTime(Particles,particlesCorrespondingDataSetNumber,...
            Spots,spotsCorrespondingDataSetNumber,currentChannel,...
            numberOfParticles,otherRelevantData,ncOfInterest,dataSetsOfInterest);
     
        % Plotting histogram (3 slices)
        quantity = 'Time on';
        currentHistogram = histogram(extrapolatedInitialTime3,...
            'FaceColor',color3(2,:),...
            'BinMethod','fd','DisplayName',quantity);
        histogramHeights = [histogramHeights currentHistogram.Values];
        binLimits = [binLimits currentHistogram.BinLimits];
       
        hold on 
        % ------ initial measured -------
        initialMeasured = findInitialMeasured(Particles,particlesCorrespondingDataSetNumber,...
            Spots,spotsCorrespondingDataSetNumber,currentChannel,...
            numberOfParticles,otherRelevantData,ncOfInterest,dataSetsOfInterest);
        
        % Plotting histogram
        quantity = 'First Data Point';
        currentHistogram = histogram(initialMeasured,...
            'FaceColor',color3(3,:),'BinMethod','fd','DisplayName',quantity);
        histogramHeights = [histogramHeights currentHistogram.Values];
        binLimits = [binLimits currentHistogram.BinLimits];
        
        % Setting X limits
        xlim([min(binLimits) max(binLimits)])
        
        % Setting Y limits and labels
        tallestBarHeight = max(histogramHeights);
        ylim([0 tallestBarHeight+1])
        desiredYRange = 0:tallestBarHeight+1;
        yAxisLabels = cellfun(@num2str, num2cell(desiredYRange),...
            'UniformOutput', false);
        set(gca, 'YTick',desiredYRange,'YTickLabel', yAxisLabels)
        quantity = 'Overlay: Time On & First Measured';
        xLabelText = 'Time (minutes)';
        legend('show')
        
    case 6 % Histogram of fitted slopes
        dataSetsOfInterest = unique(particlesCorrespondingDataSetNumber);
        allParticleIndexes = 1:numberOfParticles;
        
        particlesToInclude = allParticleIndexes(...
            ismember(particlesCorrespondingDataSetNumber,dataSetsOfInterest));
        
        allFittedSlopes = [];
        for currentParticle = particlesToInclude
            [currentParticleSubset,ParticlesSubset,SpotsSubset] = ...
                getSubsetData(currentParticle,Particles,particlesCorrespondingDataSetNumber,...
                Spots,spotsCorrespondingDataSetNumber);
            
            [frame,ampIntegral,ampIntegral3,~,~,~,~,~,~,~,~,~,~]=...
                GetParticleTrace(currentParticleSubset,...
                ParticlesSubset{currentChannel},SpotsSubset{currentChannel});
            
            % Getting corresponding frame range
            currentDataSet = particlesCorrespondingDataSetNumber(currentParticle);
            correspondingNCInfo = otherRelevantData(currentDataSet).correspondingNCInfo;
            nuclearCycleBoundariesInMinutes = otherRelevantData(currentDataSet).nuclearCycleBoundariesInMinutes;
            framesInMinutes = otherRelevantData(currentDataSet).framesInMinutes;
            frameRange = framesInMinutes(frame); % Changing the units to seconds
            ncPresent = unique(correspondingNCInfo(frame));
            timeOfFirstNC = nuclearCycleBoundariesInMinutes(ncPresent(1)-8);
            % adjusting frameRange
            frameRange = frameRange - timeOfFirstNC;
            
            % Initial Slopes
            [fittedEQ,~] = ...
                initialSlopeCalc(frameRange,ampIntegral); % one slice
            [fittedEQ3,~] = ...
                initialSlopeCalc(frameRange,ampIntegral3); % three slice
            
            if ~isnan(fittedEQ) % not including particles with only one point
                allFittedSlopes = [allFittedSlopes;
                    fittedEQ(1) fittedEQ3(1)];
            end
        end
        
        % Plotting histogram ----------------------------------------------
        % Plotting 1 slice integrated cylinder
        %         currentHistogram = histogram(allFittedSlopes(:,1),...
        %             'FaceColor',color3(2,:),'BinMethod','fd','DisplayName','One slice');
        %         histogramHeights = [histogramHeights currentHistogram.Values];
        %         binLimits = [binLimits currentHistogram.BinLimits];
        
        hold on
        % Plotting 3 slice integrated cylinder
        currentHistogram = histogram(allFittedSlopes(:,2),...
            'FaceColor',color3(3,:),'BinMethod','fd','DisplayName','Three slice');
        histogramHeights = [histogramHeights currentHistogram.Values];
        binLimits = [binLimits currentHistogram.BinLimits];
        
        % Setting X limits
        xlim([min(binLimits) max(binLimits)])
        
        % Setting Y limits and labels
        tallestBarHeight = max(histogramHeights);
        ylim([0 tallestBarHeight+1])
        desiredYRange = 0:tallestBarHeight+1;
        yAxisLabels = cellfun(@num2str, num2cell(desiredYRange),...
            'UniformOutput', false);
        set(gca, 'YTick',desiredYRange,'YTickLabel', yAxisLabels)
        quantity = 'Initial Slope';
        xLabelText = 'Slope (a.u./minutes)';
        legend('show')
        
    case 7 % Histogram of fitted slopes
        
        allParticleIndexes = 1:numberOfParticles;
        
        
        % change this so it can handle AP...multiple movies...
        for currentDataSets = dataSetsOfInterest
            allFittedSlopes = [];
            particlesToInclude = allParticleIndexes(...
            ismember(particlesCorrespondingDataSetNumber,currentDataSets));
        
            for currentParticle = particlesToInclude
                [currentParticleSubset,ParticlesSubset,SpotsSubset] = ...
                    getSubsetData(currentParticle,Particles,particlesCorrespondingDataSetNumber,...
                    Spots,spotsCorrespondingDataSetNumber);
                
                [frame,ampIntegral,ampIntegral3,~,~,~,~,~,~,~,~,~,~]=...
                    GetParticleTrace(currentParticleSubset,...
                    ParticlesSubset{currentChannel},SpotsSubset{currentChannel});
                
                % Getting corresponding frame range
                currentDataSet = particlesCorrespondingDataSetNumber(currentParticle);
                correspondingNCInfo = otherRelevantData(currentDataSet).correspondingNCInfo;
                nuclearCycleBoundariesInMinutes = otherRelevantData(currentDataSet).nuclearCycleBoundariesInMinutes;
                framesInMinutes = otherRelevantData(currentDataSet).framesInMinutes;
                frameRange = framesInMinutes(frame); % Changing the units to seconds
                ncPresent = unique(correspondingNCInfo(frame));
                timeOfFirstNC = nuclearCycleBoundariesInMinutes(ncPresent(1)-8);
                % adjusting frameRange
                frameRange = frameRange - timeOfFirstNC;
                
                % Initial Slopes
                [fittedEQ,~] = ...
                    initialSlopeCalc(frameRange,ampIntegral); % one slice
                [fittedEQ3,~] = ...
                    initialSlopeCalc(frameRange,ampIntegral3); % three slice
                
                if ~isnan(fittedEQ) % not including particles with only one point
                    allFittedSlopes = [allFittedSlopes;
                        fittedEQ(1) fittedEQ3(1)];
                end
            end
        end
        
        % Plotting histogram ----------------------------------------------
        % Plotting 1 slice integrated cylinder
        %         currentHistogram = histogram(allFittedSlopes(:,1),...
        %             'FaceColor',color3(2,:),'BinMethod','fd','DisplayName','One slice');
        %         histogramHeights = [histogramHeights currentHistogram.Values];
        %         binLimits = [binLimits currentHistogram.BinLimits];
        
        hold on
        % Plotting 3 slice integrated cylinder
        currentHistogram = histogram(allFittedSlopes(:,2),...
            'FaceColor',color3(3,:),'BinMethod','fd','DisplayName','Three slice');
        histogramHeights = [histogramHeights currentHistogram.Values];
        binLimits = [binLimits currentHistogram.BinLimits];
        
        % Setting X limits
        xlim([min(binLimits) max(binLimits)])
        
        % Setting Y limits and labels
        tallestBarHeight = max(histogramHeights);
        ylim([0 tallestBarHeight+1])
        desiredYRange = 0:tallestBarHeight+1;
        yAxisLabels = cellfun(@num2str, num2cell(desiredYRange),...
            'UniformOutput', false);
        set(gca, 'YTick',desiredYRange,'YTickLabel', yAxisLabels)
        quantity = 'Initial Slope';
        xLabelText = 'Slope (a.u./minutes)';
        legend('show')      
        
end


titleName = [quantity ' histogram'];
if makeSupTitle
    superTitle = suptitle(titleName);
    set(superTitle,'FontSize',figTitleFontSize)
else
    title(titleName,'FontSize',figTitleFontSize)
    ylabel('Counts','FontSize',axisFontSize)
    xlabel(xLabelText,'FontSize',axisFontSize)
end

end

function maxIntensity = calcualteMaxIntensityOfAll(currentChannel,Particles,Spots)
numberOfParticles = size(Particles{:},2);

maxIntensity = -1;
for currentParticle = 1:numberOfParticles
    [~,AmpIntegral,AmpIntegral3,~,~,~,~,~,~,~,~,~]=...
        GetParticleTrace(currentParticle,...
        Particles{currentChannel},Spots{currentChannel});
    currentMax = max([AmpIntegral AmpIntegral3]);
    maxIntensity = max(maxIntensity, currentMax);
end

end




function FrameInfo = createNCInformation(FrameInfo,numberOfFrames,...
    nuclearCycleBoundaries)

nc9 = nuclearCycleBoundaries(1);
nc10 = nuclearCycleBoundaries(2);
nc11 = nuclearCycleBoundaries(3);
nc12 = nuclearCycleBoundaries(4);
nc13 = nuclearCycleBoundaries(5);
nc14 = nuclearCycleBoundaries(6);

for currentFrame=1:numberOfFrames
    if currentFrame<nc9
        FrameInfo(currentFrame).nc=8;
    elseif (currentFrame>=nc9)&&(currentFrame<nc10)
        FrameInfo(currentFrame).nc=9;
    elseif (currentFrame>=nc10)&&(currentFrame<nc11)
        FrameInfo(currentFrame).nc=10;
    elseif (currentFrame>=nc11)&&(currentFrame<=nc12)
        FrameInfo(currentFrame).nc=11;
    elseif (currentFrame>=nc12)&&(currentFrame<=nc13)
        FrameInfo(currentFrame).nc=12;
    elseif (currentFrame>=nc13)&&(currentFrame<=nc14)
        FrameInfo(currentFrame).nc=13;
    elseif currentFrame>=nc14
        FrameInfo(currentFrame).nc=14;
    end
end

end

% Template call function -------------------------------------
function optionAction(source,event,particleTraceFig,traceFigHandles,...
    Particles,Spots,currentChannel,...
    markerSize,lineWidth,lineWidthFitted,color1,color3,transparencyFitted,...
    figTitleFontSize,axisFontSize)

figure(particleTraceFig);

%find the current particle
currentParticle = round(traceFigHandles.slider.Value);

%find the values of the other check boxes

disp('This button works! Congrats!!!')

end







