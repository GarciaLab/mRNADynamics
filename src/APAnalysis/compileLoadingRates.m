function [allInitialSlopes,allTimeOn,allSlopeError,allCorrespondingNC] = ...
    compileLoadingRates(data,dataSetsToInclude)
% to be used in pol2Loading Anaylsis

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
        
        nameTemp = regexp(data(currentDataSet).Prefix,'\d*','Match');
%         namesLoading{end+1} = nameTemp{end};
%         disp(nameTemp{end})
    end
    hold off 
    
    % Plotting time on vs AP ----------------------------------------------
%     figure(scatterTimeOn)
%     hold on 
    
end
end