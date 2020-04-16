function dorsalResults = createDorsalResults(DataType, varargin)

displayTiles = false;

minNuclei =1; %minimum total nuclei for a bin to be plottable
minEmbryos = 1; %minimum number of nuclei per bin

for i = 1:length(varargin)
    if strcmpi(varargin{i}, 'displayTiles')
        displayTiles = true;
        tileFig = figure();
        holdFig = figure();
    end
end

[~, resultsFolder, prefixes] = getDorsalPrefixes(DataType);

load([resultsFolder,filesep,DataType,filesep,'dlfluobins.mat'], 'dlfluobins');

ch = 1;
nEmbryos = length(prefixes);
nBins = length(dlfluobins);

binFilter =  {ones(1, nBins), ones(1, nBins), ones(1, nBins)};
embryosPerBin = {zeros(1, nBins), zeros(1, nBins), zeros(1, nBins)};
allmrnasnomean = cell(3, nBins);

npartFluoEmbryo = {};
nschnitzFluoEmbryo = {};
fracFluoEmbryo = {};
dorsalResults = {};

compiledProjects = {};
for e = 1:nEmbryos
    
    load([resultsFolder,filesep,prefixes{e},filesep,'compiledProject.mat'], 'compiledProject');
    
    compiledProjects{e} = compiledProject;
    
    for nc = 12:14
        for bin = 1:nBins
            
            dorsalResults{nc-11}.allmrnasnomean{bin} = {};
            
            tempNucleiOfInterest = [];
            spotAccumulatedFluo = [];
            spotDurations = [];
            spotMaxFluos = [];
            spotTurnOnTimes = [];
            
            nucleiOfInterest= find( [compiledProject.cycle] == nc & [compiledProject.dorsalFluoBin] == bin );
            
            
            particlesOfInterest = 0;
            for n = 1:length(nucleiOfInterest)
                
                cn = compiledProject(nucleiOfInterest(n));
                
                tempNucleiOfInterest = [tempNucleiOfInterest, nucleiOfInterest(n)];
                
                if ~isempty(cn.particleFrames)
                    particlesOfInterest = particlesOfInterest + 1;
                    spotFrames = cn.particleFrames;
                    spotFluo = cn.particleFluo3D; %can also use .Fluo for 2d fluo
                    spotDurations = [spotDurations, cn.particleDuration];
                    spotTurnOnTimes = [spotTurnOnTimes, cn.particleTimeOn];
                    spotMaxFluos = [spotMaxFluos, cn.particleFluo95];
                    spotAccumulatedFluo = [spotAccumulatedFluo, cn.particleAccumulatedFluo];
                end
                
                %                 figure(holdFig)
                %                 yyaxis left
                %                 plot(cn.nuclearFrames, cn.FluoTimeTrace, '-g');
                %                 ylim([0, max(dlfluobins)]);
                %                 yyaxis right
                %                 plot(spotFrames, spotFluo, '.-r'); %spot intensity
                %                 %                         hold on
                %                 ylim([0, 1000]);
                %                 waitforbuttonpress;
                
            end
            
            
            nucleiOfInterest = tempNucleiOfInterest;
            
            dorsalResults{nc-11}.allmrnasnomean{bin} = [dorsalResults{nc-11}.allmrnasnomean{bin}, spotAccumulatedFluo];
            
            %             npartFluoEmbryo{nc-11}(bin, e) = particlesOfInterest;
            %             nschnitzFluoEmbryo{nc-11}(bin, e) = length(nucleiOfInterest);
            
            nschnitzFluoEmbryo{nc-11}(bin, e) = length(find( [compiledProject.cycle] == nc...
                & [compiledProject.dorsalFluoBin] == bin ));
            npartFluoEmbryo{nc-11}(bin, e) = length(find([compiledProject.cycle] ==...
                nc & [compiledProject.dorsalFluoBin] == bin  & ~cellfun(@isempty,...
                {compiledProject.particleFrames})));
            dorsalResults{nc-11}.allmrnasEmbryo(bin, e) = nanmean(spotAccumulatedFluo);
            dorsalResults{nc-11}.alldurationsEmbryo(bin, e) = nanmean(spotDurations);
           dorsalResults{nc-11}. allTurnOnsEmbryo(bin, e) = nanmean(spotTurnOnTimes);
            dorsalResults{nc-11}.allMaxFluoEmbryo(bin, e) = nanmean(spotMaxFluos);
            
            %             if nschnitzFluoEmbryo{nc-11}(bin, e) >= 1
            %                 embryosPerBin{nc-11}(bin) = embryosPerBin{nc-11}(bin) + 1;
            %             end
            % %
            %             if e == 2 & bin == 3 & nc == 12 %for debugging
            %                 1
            %             end
            %
        end
        
        
        dorsalResults{nc-11}.fracFluoEmbryo(:, e) = npartFluoEmbryo{nc-11}(:,e)./nschnitzFluoEmbryo{nc-11}(:,e);
        
    end
    
end

%average across embryos
nSchnitzBinTotal = {[], []};
minNucleiFilter = {[], []};
minEmbryoFilter = {[], []};
nSchnitzBinTotalWithZeros = {[], []};
binFilter = {[], []};
nEmbryos = {[], []};
embryoDim = 2;
binDim = 1;

for nc = 1:2
    
    nSchnitzBinTotal{nc} = nansum((nschnitzFluoEmbryo{nc}), embryoDim);
    nSchnitzBinTotalWithZeros{nc} = nansum((nschnitzFluoEmbryo{nc}), embryoDim);
    nSchnitzBinTotal{nc}(~nSchnitzBinTotal{nc}) = nan; %so i can divide without getting Infs
    
    minNucleiFilter{nc} =  nSchnitzBinTotal{nc} >= minNuclei;
    minEmbryoFilter{nc} = sum(logical(nschnitzFluoEmbryo{nc}), embryoDim) >= minEmbryos;
    nEmbryos{nc} = minEmbryoFilter{nc}.*sum(logical(nschnitzFluoEmbryo{nc}), embryoDim);
    binFilter{nc} = minEmbryoFilter{nc} .* minNucleiFilter{nc};
    binFilter{nc}(~binFilter{nc}) = nan;
    
    
    
    nSamples = 100;
    
    filteredWeightedMean = @(x) ((nansum(x.*nschnitzFluoEmbryo{nc},embryoDim))./nSchnitzBinTotal{nc}).*binFilter{nc};
    filteredWeightedSE = @(y) nanstd(bootstrp(nSamples, @(x) filteredWeightedMean(x), y,'Weights',nSchnitzBinTotalWithZeros{nc})', 0, embryoDim);
    
    %     filteredWeightedMean = @(x) nanmean(x,2).*binFilter{nc}';
    %      filteredWeightedSE = @(y) nanstd(bootstrp(nSamples, @(x) filteredWeightedMean(x), y), 0, 1);
    %         filteredWeightedSE = @(x) (nanstd(x,0, 2)./sqrt(embryosPerBin{nc}'))  .* binFilter{nc}';
    
    dorsalResults{nc}.meanFracFluoEmbryo = filteredWeightedMean(dorsalResults{nc}.fracFluoEmbryo);
    dorsalResults{nc}.seFracFluoEmbryo = filteredWeightedSE(dorsalResults{nc}.fracFluoEmbryo);
    
    dorsalResults{nc}.meanallmrnasEmbryo = filteredWeightedMean(dorsalResults{nc}.allmrnasEmbryo);
    dorsalResults{nc}.seallmrnasEmbryo = filteredWeightedSE(dorsalResults{nc}.allmrnasEmbryo);
    
    dorsalResults{nc}.meanalldurationsEmbryo = filteredWeightedMean(dorsalResults{nc}.alldurationsEmbryo);
    dorsalResults{nc}.sealldurationsEmbryo = filteredWeightedSE(dorsalResults{nc}.alldurationsEmbryo);
    
    dorsalResults{nc}.meanTurnOnsEmbryo = filteredWeightedMean(dorsalResults{nc}.allTurnOnsEmbryo);
    dorsalResults{nc}.seTurnOnsEmbryo = filteredWeightedSE(dorsalResults{nc}.allTurnOnsEmbryo);
    
    dorsalResults{nc}.meanAllMaxFluoEmbryo = filteredWeightedMean(dorsalResults{nc}.allMaxFluoEmbryo);
   dorsalResults{nc}. seAllMaxFluoEmbryo = filteredWeightedSE(dorsalResults{nc}.allMaxFluoEmbryo);
   
    dorsalResults{nc}.dorsalFluoBins  = dlfluobins;
    dorsalResults{nc}.DataType = DataType;
    
end
% 
% allmrnasnc12 = dorsalResults{nc}.allmrnasnomean{:};
% lens = [];
% for b = 1:nBins
%     lens(b) = length(allmrnasnc12{b});
% end
% allmrnasnc12mat = zeros(nBins, max(lens));
% for bin = 1:nBins
%     if ~isempty(allmrnasnc12{bin})
%         allmrnasnc12mat(bin, :) = padarray(allmrnasnc12{bin}',max(lens)-length(allmrnasnc12{bin}),NaN, 'post');
%     end
% end

save([resultsFolder,filesep,DataType,filesep,'dorsalResults.mat'], 'dorsalResults', '-v6')

end
