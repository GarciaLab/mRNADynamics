close all;
% getting the prefix 
[prefix,~] = getPrefixAndFolder;

colors = [0.1961    0.8039    0.1961;...
    0.2941         0    0.7059];
color1 = colors(1,:);
color3 = colors(2,:);
transparencyFitted = 0.5;

markerSize = 25;
lineWidth = 3;
fittedFactor = 2/3;
lineWidthFitted = lineWidth*fittedFactor;
waitTime = 0.1;
figTitleFontSize  = 20;
axisFontSize = 14;

%% Loading the data 
[~,~,dropboxFolder,~,~]= DetermineLocalFolders(prefix);
[~,~,defaultDropboxFolder,~,PreProcPath]=...
    DetermineLocalFolders;
dataFolder=[dropboxFolder,filesep,prefix];
FilePrefix=[dataFolder(length(dropboxFolder)+2:end),'_'];

% Particle Information and Loading
particlePathName = [dataFolder,filesep,'Particles.mat'];
load(particlePathName)
if ~iscell(Particles)
    Particles = {Particles};
    SpotFilter = {SpotFilter};
end
numberOfParticles = size(Particles{:},2);
currentChannel = 1;

% Spots Information and Loading
spotsPathName = [dataFolder,filesep,'Spots.mat'];
load(spotsPathName)
if ~iscell(Spots)
    Spots = {Spots};
end

particleTraceFig = figure();
positionOfParticleTraceFig = [4 49 996 635];%[4 97 996 587];%[403 49 562 635];
set(gcf,'Position',positionOfParticleTraceFig);

for currentParticle = 1:numberOfParticles
    drawCurrentTrace(currentParticle,Particles,Spots,currentChannel,...
        markerSize,lineWidth,color1,color3,figTitleFontSize,axisFontSize)
    
    pause(waitTime)
end
minY = 0;
maxY = calcualteMaxIntensityOfAll(currentChannel,Particles,Spots);
ylim([minY maxY*1.03])

% making slider
sliderPosition = [0.9268    0.2000    0.0315    0.7000]; % Units : normalized
% textPosition = [0.91 0.9 0.04 0.9100]; % Units : normalized
sliderSmallStep = 1/(numberOfParticles-1); % clicking will change it by 1
sliderBigStep = sliderSmallStep*4; % sliding will change it by units of 4
% Create slider
slider = uicontrol('Style', 'slider',...
    'Min',1,'Max',numberOfParticles,'Value',currentParticle,...
    'Units', 'normalized','Tag','slider',...
    'Position', sliderPosition,...
    'SliderStep',[sliderSmallStep sliderBigStep],...
    'Callback', {@chooseParticle,Particles,Spots,currentChannel,...
    markerSize,lineWidth,lineWidthFitted,color1,color3,transparencyFitted,...
    figTitleFontSize, axisFontSize});

% Add a text uicontrol to label the slider.
% txt = uicontrol('Style','text',...
%     'Units', 'normalized',...
%     'Position',textPosition,...
%     'String','test');

traceFigHandles = guihandles(particleTraceFig); 


%% making option figure 
optionsFig = figure();
positionOfOptionFig = [1017 48 255 635];%[1017 97 255 586];
maxButtonYPoint = 0.8318;
firstXPoint = 0.0745;
fullButtonWidth = 0.7;
buttonHeight = 0.0341;
tabFactor = 1.1;
checkBoxOneIndentFactor = 3;
sectionSpacing = 0.1;
currentSection = 0;
buttonSpacing = buttonHeight;
defaultOption = 0;

%% First Section 
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
sectionTitleString = 'Options';
sectionFontSize = 12;
sectionString = 'Metric Graphing Options';
sectionTitleTxt = uicontrol('Style','text',...
    'Units', 'normalized',...
    'Position',textPosition,...
    'FontSize',sectionFontSize,...
    'String',sectionString);

buttonCounter = 1;
set(gcf,'Position',positionOfOptionFig);
initialSlopeOption = uicontrol('Style','checkbox',...
    'String',{'Initial Slope'},'Min',defaultOption,'Max',buttonCounter,...
    'Value',defaultOption,'Tag','initialSlopeButton',...
    'Callback', {@initialSlopeOptionAction,particleTraceFig,...
    traceFigHandles,...
    Particles,Spots,currentChannel,...
    markerSize,lineWidth,lineWidthFitted,color1,color3,transparencyFitted});
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
    'Callback', {@toggleResizingOptionAction,particleTraceFig,...
    traceFigHandles,...
    Particles,Spots,currentChannel,...
    markerSize,lineWidth,lineWidthFitted,color1,color3,transparencyFitted,...
    figTitleFontSize, axisFontSize});
toggleRespacingOption.Units = 'normalized';
currentButtonYPoint = maxButtonYPoint- buttonSpacing*(buttonCounter-1)...
    - currentSection*sectionSpacing;
% indenting this checkbox
currentButtonXPoint = firstXPoint + buttonSpacing*checkBoxOneIndentFactor;
toggleRespacingOption.Position = [currentButtonXPoint currentButtonYPoint...
    fullButtonWidth buttonHeight];

% %% New section
% currentSection = currentSection + 1;
% 
% % Add a text for the section title
% textHeight = 0.0341;
% textLength = 0.9;
% currentYPoint = maxButtonYPoint- buttonSpacing*(buttonCounter-1)...
%     - currentSection*sectionSpacing;
% textPosition = [firstXPoint currentYPoint+textHeight ...
%     textLength textHeight*1.1]; 
% sectionTitleString = 'Options';
% sectionFontSize = 12;
% sectionString = 'Metric Graphing Options';
% sectionTitleTxt = uicontrol('Style','text',...
%     'Units', 'normalized',...
%     'Position',textPosition,...
%     'FontSize',sectionFontSize,...
%     'String',sectionString);
% 
% 
% % Initial Time Histogram button
% % buttonCounter = buttonCounter + 1;
% initialTimeHistogram = uicontrol('Style','checkbox',...
%     'String','Initial Time','Min',defaultOption,'Max',buttonCounter,...
%     'Value',defaultOption,'Tag','selectButton',...
%     'Callback', {@initialTimeHistogramOptionAction,particleTraceFig,traceFigHandles,...
%     Particles,Spots,currentChannel,...
%     markerSize,lineWidth,lineWidthFitted,color1,color3,transparencyFitted,...
%     figTitleFontSize, axisFontSize});
% initialTimeHistogram.Units = 'normalized';
% currentButtonYPoint = maxButtonYPoint- buttonSpacing*(buttonCounter-1)...
%     - currentSection*sectionSpacing;
% currentButtonXPoint = firstXPoint;
% initialTimeHistogram.Position = [currentButtonXPoint currentButtonYPoint...
%     fullButtonWidth buttonHeight];
% 
% % Initial Slope Histogram button
% buttonCounter = buttonCounter + 1;
% initialSlopeHistogram = uicontrol('Style','checkbox',...
%     'String','Initial Slope','Min',defaultOption,'Max',buttonCounter,...
%     'Value',defaultOption,'Tag','selectButton',...
%     'Callback', {@optionAction,particleTraceFig,traceFigHandles,...
%     Particles,Spots,currentChannel,...
%     markerSize,lineWidth,lineWidthFitted,color1,color3,transparencyFitted,...
%     figTitleFontSize, axisFontSize});
% initialSlopeHistogram.Units = 'normalized';
% currentButtonYPoint = maxButtonYPoint- buttonSpacing*(buttonCounter-1)...
%     - currentSection*sectionSpacing;
% currentButtonXPoint = firstXPoint;
% initialSlopeHistogram.Position = [currentButtonXPoint currentButtonYPoint...
%     fullButtonWidth buttonHeight];
% 
% 
% 
% 
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

%% other functions (most will be placed in a new script)
% redrawing the trace of the particle that was choosen
function chooseParticle(source,event,Particles,Spots,currentChannel,...
    markerSize,lineWidth,lineWidthFitted,color1,color3,transparencyFitted,...
    figTitleFontSize,axisFontSize)
currentParticle = round(source.Value);

drawCurrentTrace(currentParticle,Particles,Spots,currentChannel,...
    markerSize,lineWidth,color1,color3,figTitleFontSize,axisFontSize)

end

% draw the current particle's trace
function drawCurrentTrace(currentParticle,Particles,Spots,currentChannel,...
    markerSize,lineWidth,color1,color3,figTitleFontSize,axisFontSize)

[Frame,AmpIntegral,AmpIntegral3,~,~,~,~,~,~,~,~,~]=...
    GetParticleTrace(currentParticle,...
    Particles{currentChannel},Spots{currentChannel});

% plotting the traces
cla

plot(Frame,AmpIntegral,'Marker','.','MarkerSize',markerSize,...
    'LineWidth',lineWidth,'Color',color1,'DisplayName','one slice')
hold on
plot(Frame,AmpIntegral3,'Marker','.','MarkerSize',markerSize,...
    'LineWidth',lineWidth,'Color',color3,'DisplayName','three slice')

legend();
title(['Particle ' num2str(currentParticle)],'FontSize',figTitleFontSize)
xlabel('Frame','FontSize',axisFontSize)
ylabel('Integrated Intensity','FontSize',axisFontSize)


xPadding = 0.3*length(Frame);
minX = Frame(1)-xPadding;
maxX = Frame(end)+xPadding;
xlim([minX maxX])

% minY = 0;
% maxY = calcualteMaxIntensityOfAll(Particles,Spots);
% ylim([minY maxY])
end

function initialSlopeOptionAction(source,event,particleTraceFig,traceFigureHandle,...
    Particles,Spots,currentChannel,...
    markerSize,lineWidth,lineWidthFitted,color1,color3,transparencyFitted...
    ,figTitleFontSize,axisFontSize)

figure(particleTraceFig);

%find the current particle
currentParticle = round(traceFigureHandle.slider.Value);

%find the values of the other check boxes

drawCurrentTrace(currentParticle,Particles,Spots,currentChannel,...
    markerSize,lineWidth,color1,color3,figTitleFontSize,axisFontSize)

if source.Value
    rescaling = 0;
    plotFittedLine(currentParticle,Particles,Spots,currentChannel,...
                lineWidthFitted,color1,color3,transparencyFitted,rescaling);
end

end

function toggleResizingOptionAction(source,event,particleTraceFig,traceFigureHandle,...
    Particles,Spots,currentChannel,...
    markerSize,lineWidth,lineWidthFitted,color1,color3,transparencyFitted,...
    figTitleFontSize,axisFontSize)

figure(particleTraceFig);

%find the current particle
currentParticle = round(traceFigureHandle.slider.Value);

%find the values of the other check boxes

drawCurrentTrace(currentParticle,Particles,Spots,currentChannel,...
    markerSize,lineWidth,color1,color3,figTitleFontSize,axisFontSize)

if source.Value
    rescaling = 1;
    plotFittedLine(currentParticle,Particles,Spots,currentChannel,...
                lineWidthFitted,color1,color3,transparencyFitted,rescaling);
else
    rescaling = 0;
    plotFittedLine(currentParticle,Particles,Spots,currentChannel,...
                lineWidthFitted,color1,color3,transparencyFitted,rescaling);
end

end

% function initialTimeHistogramOptionAction(source,event,particleTraceFig,traceFigHandles,...
%     Particles,Spots,currentChannel,...
%     markerSize,lineWidth,lineWidthFitted,color1,color3,transparencyFitted,...
%     figTitleFontSize, axisFontSize});
%    
% 
% 
% end



function optionAction(source,event,particleTraceFig,traceFigureHandle,...
    Particles,Spots,currentChannel,...
    markerSize,lineWidth,lineWidthFitted,color1,color3,transparencyFitted,...
    figTitleFontSize,axisFontSize)

figure(particleTraceFig);

%find the current particle
currentParticle = round(traceFigureHandle.slider.Value);

%find the values of the other check boxes

disp('This button works! Congrats!!!')

end


% calculate the initial slope of the current particle
function [coefficientsOfFittedLine,indexOfMax] =...
    initialSlopeCalc(frames,intensities)

[~, indexOfMax]= max(intensities);
regionOfRiseX = frames(1:indexOfMax);
regionOfRiseY = intensities(1:indexOfMax);

if length(regionOfRiseX) == 1
    coefficientsOfFittedLine = NaN;
else
    coefficientsOfFittedLine = polyfit(regionOfRiseX,regionOfRiseY,1);
end

end

% plot the fitted line
function plotFittedLine(currentParticle,Particles,Spots,currentChannel,...
    lineWidthFitted,color1,color3,transparencyFitted,rescaling)

[Frame,AmpIntegral,AmpIntegral3,~,~,~,~,~,~,~,~,~]=...
    GetParticleTrace(currentParticle,...
    Particles{currentChannel},Spots{currentChannel});

% Initial Slopes
[fittedEQ,indexOfMax] = ...
    initialSlopeCalc(Frame,AmpIntegral); % one slice
[fittedEQ3,indexOfMax3] = ...
    initialSlopeCalc(Frame,AmpIntegral3); % three slice

% plotting on the same figure
hold on

% plotting fitted lines
if ~isnan(fittedEQ)
    xSlope = Frame(1:indexOfMax);
    extrapolatedInitialTime = roots(fittedEQ);
    yFitted = polyval(fittedEQ,xSlope);
    
    xSlope3 = Frame(1:indexOfMax3);
    extrapolatedInitialTime3 = roots(fittedEQ3);
    yFitted3 = polyval(fittedEQ3,xSlope3);
    
    temp = plot([extrapolatedInitialTime xSlope],[0 yFitted],...
        'LineWidth',lineWidthFitted,...
        'Color',color1,'DisplayName','Fitted');
    temp.Color(4) = transparencyFitted;
    temp = plot([extrapolatedInitialTime3 xSlope3],[0 yFitted3],...
        'Color',color3,'DisplayName','Fitted3');
    temp.Color(4) = transparencyFitted;
    
    furtherPoint = min(extrapolatedInitialTime,extrapolatedInitialTime3);
    if rescaling
        maxXLimit = Frame(1)-furtherPoint + Frame(end);
        xlim([furtherPoint maxXLimit])
    else
        xPadding = 0.3*length(Frame);
        xlim([furtherPoint Frame(end)+xPadding])
        
    end
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






