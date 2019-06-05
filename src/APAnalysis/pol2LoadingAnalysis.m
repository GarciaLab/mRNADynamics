function pol2LoadingAnalysis(dataSet)

allColors = [0.9412    0.2353    0.2353;...
    0.9412    0.6118    0.2353;...
    0.9412    0.8941    0.2353;...
    0.5176    0.9412    0.2353;...
    0.2353    0.9412    0.3765;...
    0.2353    0.8941    0.9412;...
    0.2353    0.6588    0.9412;...
    0.2353    0.3294    0.9412;...
    0.3765    0.2353    0.9412;...
    0.5647    0.2353    0.9412;...
    0.7529    0.2353    0.9412;...
    0.9412    0.2353    0.8941;...
    1.0000    0.4392    0.8118;...
    0.9412    0.2353    0.6118;...
    1.0000    0.5804    0.5804;...
    1.0000    0.8314    0.5804;...
    0.8863    1.0000    0.5804;...
    0.5804    1.0000    0.6902;...
    0.5804    1.0000    0.9451;...
    0.5804    0.8863    1.0000];

data = LoadMS2Sets(dataSet);

% Checking to see which data sets has fittedLineEquations
% the data set will not be plotted if it does not have this variables
% and associated fields!
dataSetsToInclude = [];
for currentDataSet = 1:length(data)
    if isfield(data(currentDataSet),'fittedLineEquations')
        dataSetsToInclude = [dataSetsToInclude currentDataSet];
    end
end

%% Initializing and positioning figures
% figure for scatter of loading rate vs AP
scatterLoading = figure();
set(scatterLoading,'Position',[4.5000 43.5000 679 601])

% figure for box and whisker of loading rate vs AP
boxLoading = figure();
set(boxLoading,'Position',[693 43.5000 585.5000 215])

% figure for avg loading rate vs AP
% avgLoading = figure();
% set(avgLoading,'Position',[692 341 583 303.5000])

% figure for time on vs AP
% scatterTimeOn = figure();

% figure for avg time on vs AP
% avgTimeOn = figure();
% set(avgTimeOn,'Position',[691.5000 43.5000 583 303.5000])


%% compiling data and plotting on continous AP axis (all nc, scatter)
plotsLabeledLoading = plot(NaN,NaN);
namesLoading = {'placeHolder'};
allAPPositions = [];
allInitialSlopes = [];
allSlopeError = []; % the error is taken to be the norm of the residuals
allTimeOn = [];
allCorrespondingNC = [];

for currentDataSet = dataSetsToInclude
    numberOfParticles = length(data(currentDataSet).Particles);
    labelPlot = 0; %condition
    clear apPositions
    currentInitialSlopes = NaN(1,numberOfParticles);
    currentTimeOn = NaN(1,numberOfParticles);
    errorEstimations = NaN(1,numberOfParticles);
    currentNuclearCycle = NaN(1,numberOfParticles); % stores the corresponding nuclear cycle of the particle
    
    anaphaseBoundaries = [data(currentDataSet).nc9,data(currentDataSet).nc10,...
        data(currentDataSet).nc11,data(currentDataSet).nc12,...
        data(currentDataSet).nc13,data(currentDataSet).nc14];
    
    for currentParticle = 1:numberOfParticles
        apPositions(currentParticle) = mean(data(currentDataSet).Particles(currentParticle).APpos);
        firstFrame = data(currentDataSet).Particles(currentParticle).Frame(1);
        [~,sortedIndex] = sort([anaphaseBoundaries firstFrame]);
        tempIndex = find(sortedIndex == length(anaphaseBoundaries)+1);
        currentNuclearCycle(currentParticle) = sortedIndex(tempIndex-1) + 8;
        
        try
            tempSlope = ...
                data(currentDataSet).fittedLineEquations(currentParticle).Coefficients(1,1);
            tempTimeOn = roots(data(currentDataSet).fittedLineEquations(currentParticle).Coefficients(1,:));
            %             disp([num2str(currentDataSet) ', ' num2str(currentParticle) ': ' ...
            %                 num2str(tempTimeOn)])
            if tempSlope >= 0
                currentInitialSlopes(currentParticle) = tempSlope;
                errorEstimations(currentParticle) =...
                    data(currentDataSet).fittedLineEquations(currentParticle).ErrorEstimation(1).normr/...
                    data(currentDataSet).fittedLineEquations(currentParticle).numberOfParticlesUsedForFit(1);
                currentTimeOn(currentParticle) = tempTimeOn;
            end
        end
    end
    
    %Storing these values
    allAPPositions(end+1:end+numberOfParticles) = apPositions;
    allInitialSlopes(end+1:end+numberOfParticles) = currentInitialSlopes;
    allSlopeError(end+1:end+numberOfParticles) = errorEstimations;
    allTimeOn(end+1:end+numberOfParticles) = currentTimeOn;
    allCorrespondingNC(end+1:end+numberOfParticles) = currentNuclearCycle;
    
    %Plotting loading rate vs AP -----------------------------------------
    %Plotting the error bars
    figure(scatterLoading)
    hold on
    %     apPositions = 0.3.*ones(1,length(currentInitialSlopes));
    errorbar(apPositions, currentInitialSlopes,errorEstimations,'.',...
        'Color','black');%allColors(currentDataSet,:));
    %Plotting the values
    plotsLabeledLoading(end+1) = plot(apPositions,currentInitialSlopes,'.','MarkerSize',15,...
        'Color',allColors(currentDataSet,:));
    
    if sum(isnan(currentInitialSlopes)) ~= length(data(currentDataSet).Particles)
        hold on
        meanInitialSlope = nanmean(currentInitialSlopes);
        
        apBinOfMovie = data(currentDataSet).APbinID(max(data(currentDataSet).TotalEllipsesAP,[],2)>0);
        startX = min(apBinOfMovie);
        endX = max(apBinOfMovie);
        lengthX = endX-startX;
        paddedBoundary = 0.05;
        rectanglePosition = [ startX-paddedBoundary meanInitialSlope*0.99 ...
            lengthX+2*paddedBoundary meanInitialSlope*0.01];
        rectangle('Position',rectanglePosition,'Curvature',0.2,...
            'FaceColor',allColors(currentDataSet,:),...
            'EdgeColor',allColors(currentDataSet,:));
        
        %         plotsLabeled(end+1) = plot(nanmean(apPositions), nanmean(currentInitialSlopes),...
        %             '.','MarkerSize',30,'Color',allColors(currentDataSet,:));%h.Color);
        
        %         nameTemp = regexp(data(currentDataSet).Prefix,'\d*','Match');
        %         namesLoading{end+1} = nameTemp{end};
        %         disp(nameTemp{end})
    end
    hold off
    
    % Plotting time on vs AP ----------------------------------------------
    %     figure(scatterTimeOn)
    %     hold on
    
end
% labeling and titles
legend(plotsLabeledLoading(2:end),namesLoading(2:end), 'Interpreter', 'none');
xlim([0 1])
xlabel('Embryo Length (%)')
ylabel('Initial Rate (a.u./min)')
title('rate across AP for 1A3v7,all nc')

%% rate vs ap scatter with different nuclear cycles
ncOfInterest = [12 13 14];
for currentNC = ncOfInterest
    rateHistFig = figure();
    rateHistAxes = axes(rateHistFig);
    plotsLabeledLoading1 = plot(NaN,NaN);
    namesLoading = {'placeHolder'};
    %finding the indexes of the spots in the current nuclear cycle
    currentNCOnly = allCorrespondingNC == currentNC;
    
    %Plotting loading rate vs AP -----------------------------------------
    %Plotting the error bars
    
    errorbar(rateHistAxes,allAPPositions(currentNCOnly),...
        allInitialSlopes(currentNCOnly),allSlopeError(currentNCOnly),...
        '.','Color','black');
    hold('on')
    %Plotting the values
    plotsLabeledLoading1(end+1) =...
        plot(rateHistAxes,allAPPositions(currentNCOnly),...
        allInitialSlopes(currentNCOnly),'.','MarkerSize',15,...
        'Color','red');
    
    %         if ~isempty(allInitialSlopes(currentNCOnly))
    %             hold on
    %             meanInitialSlope = nanmean(allInitialSlopes(currentNCOnly));
    %
    %             apBinOfMovie = data(currentDataSet).APbinID(max(data(currentDataSet).TotalEllipsesAP,[],2)>0);
    %             startX = min(apBinOfMovie);
    %             endX = max(apBinOfMovie);
    %             lengthX = endX-startX;
    %             paddedBoundary = 0.05;
    %             rectanglePosition = [ startX-paddedBoundary meanInitialSlope*0.99 ...
    %                 lengthX+2*paddedBoundary meanInitialSlope*0.01];
    %             rectangle('Position',rectanglePosition,'Curvature',0.2,...
    %                 'FaceColor',allColors(currentDataSet,:),...
    %                 'EdgeColor',allColors(currentDataSet,:));
    %
    %             %         plotsLabeled(end+1) = plot(nanmean(apPositions), nanmean(currentInitialSlopes),...
    %             %             '.','MarkerSize',30,'Color',allColors(currentDataSet,:));%h.Color);
    %
    %             %         nameTemp = regexp(data(currentDataSet).Prefix,'\d*','Match');
    %             %         namesLoading{end+1} = nameTemp{end};
    %             %         disp(nameTemp{end})
    %         end
    hold off
    
    % Plotting time on vs AP ----------------------------------------------
    %     figure(scatterTimeOn)
    %     hold on
    
    legend(plotsLabeledLoading1(2:end),namesLoading(2:end), 'Interpreter', 'none');
    xlim([0 1])
    xlabel('Embryo Length (%)')
    ylabel('Initial Rate (a.u./min)')
    rateScatterTitle = ['rate across AP scatter for 1A3v7, nuclear cycle ',...
        num2str(currentNC)];
    title(gca, rateScatterTitle);
    set(gcf, 'NumberTitle', 'off','Name', rateScatterTitle);
end


%% plotting with binned ap
ap = data(1).APbinID;
apBinWidth = ap(2)-ap(1);
% apMid = ap+apBinWidth/2;
% apMidString = cellstr(string(ap+apBinWidth/2));
numAPBins = length(ap);

ncOfInterest = [12 13 14];
for currentNC = ncOfInterest
    meanInitialRateAP = NaN(1,numAPBins);
    seInitialRateAP = NaN(1,numAPBins);
    meanTimeOnAP = NaN(1,numAPBins);
    apBinGrouping = NaN(1,length(allInitialSlopes)); % stores corresponding bin string for allInitialSlopes
    
    % finding the indexes of the ones in the current nuclear cycle
    currentNCOnly = allCorrespondingNC == currentNC;
    subSetSlopeData = allInitialSlopes(currentNCOnly);
    subSetTimeOnData = allTimeOn(currentNCOnly);
    
    for currentAPBinIndex = 2:numAPBins % index of upper bound of ap
        apLowerBound = ap(currentAPBinIndex-1);
        apUpperBound = ap(currentAPBinIndex);
        apPositionsIncluded = (allAPPositions(currentNCOnly)>apLowerBound) & (allAPPositions(currentNCOnly)<apUpperBound);
        apBinGrouping(apPositionsIncluded) = round(mean([apLowerBound,apUpperBound]),2,'significant');
        
        meanInitialRateAP(currentAPBinIndex-1) = nanmean(subSetSlopeData(apPositionsIncluded));
        denomLoading = sqrt(sum(~isnan((subSetSlopeData(apPositionsIncluded)))));
        
        meanTimeOnAP(currentAPBinIndex-1) = nanmean(subSetTimeOnData(apPositionsIncluded));
        denomTimeOn = sqrt(sum(~isnan((subSetTimeOnData(apPositionsIncluded)))));
        
        if denomLoading==0
            denomLoading = 1;
        end
        seInitialRateAP(currentAPBinIndex-1) = nanstd(subSetSlopeData(apPositionsIncluded))/denomLoading; %au/min
        
        if denomTimeOn==0
            denomTimeON = 1;
        end
        seTimeOnAP(currentAPBinIndex-1) = nanstd(subSetTimeOnData(apPositionsIncluded))/denomTimeON; %min
        
        % need labeling? 
        % Do we need the bottom one?
    end
    
    %% rate vs. ap curve figure
    rateFig = figure();
    rateAxes = axes(rateFig);
    % % bar(ap+apBinWidth/2,meanInitialRateAP);
    idx = ~any(isnan(seInitialRateAP),1);
    x = ap+apBinWidth/2;
    er = errorbar(rateAxes,x(idx),meanInitialRateAP(idx), seInitialRateAP(idx));
    xlim(rateAxes,[.275,.65])
    ylim(rateAxes,[0, 1000])
    xlabel(rateAxes,'fraction anterior-posterior')
    ylabel(rateAxes,'pol II loading rate(a.u./min)')
    standardizeFigure(rateAxes, [], 'fontSize', 14)
    set(er, 'LineStyle', '-')
    rateTitleLine = ['rate across AP, nc', num2str(currentNC)];
    title(rateAxes, rateTitleLine)
    set(rateFig, 'NumberTitle', 'off','Name', rateTitleLine);

    
    %% time on vs AP figure
    timeOnFig = figure();
    timeOnAxes = axes(timeOnFig);
    idXTimeOn = ~any(isnan(seTimeOnAP),1);
    x = ap+apBinWidth/2;
    erTimeON = errorbar(timeOnAxes,x(idXTimeOn),meanTimeOnAP(idXTimeOn), seTimeOnAP(idXTimeOn));
    xlim(timeOnAxes,[.275 .65])
    ylim(timeOnAxes,[0 20])
    xlabel(timeOnAxes,'fraction anterior-posterior')
    ylabel(timeOnAxes,'time on (min)')
    standardizeFigure(timeOnAxes, [], 'fontSize', 14)
    set(erTimeON, 'LineStyle', '-')
    timeTitleLine = ['time on across AP, nc', num2str(currentNC)];
    title(timeOnAxes,timeTitleLine);
    set(timeOnFig, 'NumberTitle', 'off','Name', timeTitleLine);

end

%% rate vs ap BOX PLOT
currentNC = 12;
figure(boxLoading)
boxplot(allInitialSlopes,apBinGrouping,'PlotStyle','Compact')
boxTitle = ['rate across AP for 1A3v7, nuclear cycle ',...
    num2str(currentNC)];
title(gca, boxTitle);
set(gcf, 'NumberTitle', 'off','Name', boxTitle);

%% time on histogram figure (nc12,13,14)
ncOfInterest = [12 13 14];
for currentNC = ncOfInterest
    
    timeOnHistFig = figure();
    timeOnHistAxes = axes(timeOnHistFig);
    
    currentNCOnly = allCorrespondingNC == currentNC;
    cycleTime12 = 10; %mins
    nBins = cycleTime12*60 / 60; %one minute bins
    
    timeOnHist = histogram(timeOnHistAxes,allTimeOn(currentNCOnly), 'normalization', 'pdf', 'normalization', 'pdf', 'NumBins', nBins);
    
    xlabel(timeOnHistAxes,'time on (min)')
    ylabel(timeOnHistAxes,'frequency')
    timeOnHistTitle = {'time on frequency distribution for 1A3v7'; ['nuclear cycle '...
        num2str(currentNC) ', all AP bins']};
    title(timeOnHistAxes, timeOnHistTitle);
    xlim(timeOnHistAxes, [0, cycleTime12]) %mins
    standardizeFigure(timeOnHistAxes, [])
    %set(timeOnHistFig, 'NumberTitle', 'off','Name', timeOnHistTitle);
    
end

%% pol2 loading rate histogram (nc12,13,14)
ncOfInterest = [12 13 14];
for currentNC = ncOfInterest
    rateHistFig = figure();
    currentNCOnly = allCorrespondingNC == currentNC;
    rateHistAxes = axes(rateHistFig);
    rateHist = histogram(rateHistAxes,...
        allInitialSlopes(currentNCOnly), 'normalization', 'pdf', 'BinWidth', 50);
    hold on
    height = max(ylim);
    plot(rateHistAxes,...
        [1 1].*nanmean(allInitialSlopes(currentNCOnly)),[0 1].*height,...
        'DisplayName','Mean')
    plot(rateHistAxes,...
        [1 1].*nanmedian(allInitialSlopes(currentNCOnly)),[0 1].*height,...
        'DisplayName','Median')
    xlabel(rateHistAxes,'pol2 loading rate (a.u.)')
    ylabel(rateHistAxes,'frequency')
    rateHistTitle = {'loading rates: 1A3v7'; ['nuclear cycle '...
        num2str(currentNC) ', all AP bins']};
    title(rateHistAxes,rateHistTitle);
    standardizeFigure(rateHistAxes, [])
    legend('show')
    set(rateHistFig, 'NumberTitle', 'off','Name', rateHistTitle);

end

end

