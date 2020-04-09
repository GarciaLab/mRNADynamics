function plotFracByDlFluo2(DataType, varargin)

fractionFlag = false;
mrnaFlag = false;
durationFlag =false;
maxFlag =false;
turnOnFlag = false;


displayTiles = false;

minNuclei = 3; %minimum total nuclei for a bin to be plottable
minEmbryos = 2; %minimum number of nuclei per bin

for i = 1:length(varargin)
    if strcmpi(varargin{i}, 'fraction')
        fractionFlag = true;
    elseif strcmpi(varargin{i}, 'mrna')
        mrnaFlag = true;
    elseif strcmpi(varargin{i}, 'duration')
        durationFlag = true;
    elseif strcmpi(varargin{i}, 'maxfluo')
        maxFlag = true;
    elseif strcmpi(varargin{i}, 'turnOn') || strcmpi(varargin{i}, 'timeOn')
        turnOnFlag = true;
    elseif strcmpi(varargin{i}, 'displayTiles')
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

compiledProjects = {};
for e = 1:nEmbryos
    
    load([resultsFolder,filesep,prefixes{e},filesep,'compiledProject.mat'], 'compiledProject');
    
    compiledProjects{e} = compiledProject;
    
    for nc = 12:14
        for bin = 1:nBins
            
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
            
            allmrnasnomean{nc-11,bin} = [allmrnasnomean{nc-11,bin}, spotAccumulatedFluo];
            
            %             npartFluoEmbryo{nc-11}(bin, e) = particlesOfInterest;
            %             nschnitzFluoEmbryo{nc-11}(bin, e) = length(nucleiOfInterest);
            nschnitzFluoEmbryo{nc-11}(bin, e) = length(find( [compiledProject.cycle] == nc & [compiledProject.dorsalFluoBin] == bin ));
            npartFluoEmbryo{nc-11}(bin, e) = length(find([compiledProject.cycle] == nc & [compiledProject.dorsalFluoBin] == bin  & ~cellfun(@isempty, {compiledProject.particleFrames})));
            allmrnasEmbryo{nc-11}(bin, e) = nanmean(spotAccumulatedFluo);
            alldurationsEmbryo{nc-11}(bin, e) = nanmean(spotDurations);
            allTurnOnsEmbryo{nc-11}(bin, e) = nanmean(spotTurnOnTimes);
            allMaxFluoEmbryo{nc-11}(bin, e) = nanmean(spotMaxFluos);
            
            %             if nschnitzFluoEmbryo{nc-11}(bin, e) >= 1
            %                 embryosPerBin{nc-11}(bin) = embryosPerBin{nc-11}(bin) + 1;
            %             end
            % %
            %             if e == 2 & bin == 3 & nc == 12 %for debugging
            %                 1
            %             end
            %
        end
        
        
        fracFluoEmbryo{nc-11}(:, e) = npartFluoEmbryo{nc-11}(:,e)./nschnitzFluoEmbryo{nc-11}(:,e);
        
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
    
    meanFracFluoEmbryo{nc} = filteredWeightedMean(fracFluoEmbryo{nc});
    seFracFluoEmbryo{nc} = filteredWeightedSE(fracFluoEmbryo{nc});
    
    meanallmrnasEmbryo{nc} = filteredWeightedMean(allmrnasEmbryo{nc});
    seallmrnasEmbryo{nc} = filteredWeightedSE(allmrnasEmbryo{nc});
    
    meanalldurationsEmbryo{nc} = filteredWeightedMean(alldurationsEmbryo{nc});
    sealldurationsEmbryo{nc} = filteredWeightedSE(alldurationsEmbryo{nc});
    
    meanTurnOnsEmbryo{nc} = filteredWeightedMean(allTurnOnsEmbryo{nc});
    seTurnOnsEmbryo{nc} = filteredWeightedSE(allTurnOnsEmbryo{nc});
    
    meanAllMaxFluoEmbryo{nc} = filteredWeightedMean(allMaxFluoEmbryo{nc});
    seAllMaxFluoEmbryo{nc} = filteredWeightedSE(allMaxFluoEmbryo{nc});
    
end

allmrnasnc12 = allmrnasnomean(1, :);
lens = [];
for b = 1:nBins
    lens(b) = length(allmrnasnc12{b});
end
allmrnasnc12mat = zeros(nBins, max(lens));
for bin = 1:nBins
    if ~isempty(allmrnasnc12{bin})
        allmrnasnc12mat(bin, :) = padarray(allmrnasnc12{bin}',max(lens)-length(allmrnasnc12{bin}),NaN, 'post');
    end
end


dvLim = [0, 3000];
% dvLim = ([0, max(dlfluobins)*1.1]);
%%
%plotting
if fractionFlag
    figure()
    axs = {};
    for cycle = 1:1
        
        %     axs{cycle} = subplot(1, 2, cycle);
        x = dlfluobins(:);
        y = meanFracFluoEmbryo{cycle}(:);
        z = seFracFluoEmbryo{cycle}(:);
        idx = ~any(isnan(y),2);
        x = x(idx); y = y(idx); z = z(idx);
%         base = 3.6;
%         logbase = @(x, b) log(x)/log(b);
        nObs = nEmbryos{nc}(idx);
        nObs = nObs(:);
%         weights = 1./ (nObs.*(z.^2));
        weights = ones(length(z), 1);
        [fit, modelW] = fitDorsalActivity(x, y, z, nObs);
        errorbar(x, y , z, '-o');
        hold on
        plot(x, modelW(fit,x./weights));
%         errorbar(logbase(dlfluobins(idx), base), yData(idx) , zData(idx), '-o');
        xlabel('dorsal concentration (au)');
        ylabel('fraction active nuclei');
%         set(gca, 'XScale', 'linear')
            set(gca, 'XScale', 'log')
        ylim([0, 1]);
        %         xlim([0, max(dlfluobins)*1.1]);
        title([DataType, ' nc',num2str(cycle+11)]);
        standardizeFigure(gca, []);
        xlim(dvLim);
        legend('data', 'fit')
%         plotDorsalActivity(dlfluobins, meanFracFluoEmbryo{cycle}, seFracFluoEmbryo{cycle}, nEmbryos{nc}, 'fraction competent', 1)
%         ax = gca;
%         ax.Children.MarkerSize = 3;
%         ax.Children.LineStyle = 'none';
%         ax.Children.LineWidth = 3;
%         ylim([0, 1])
%         
%         xticks = 1:length(dlfluobins);
%         set(gca, 'XTick', xticks);
%         
%         
%         
%         for j = 1:length(xticks)
%           xtl{j} = [num2str(base) '^' num2str(xticks(j))];
%         end
%         set(gca, 'XTickLabel', xtl)
%         mind = min(logbase(dlfluobins(idx), base));
%         maxd = max(logbase(dlfluobins(idx), base));
%         xlim([mind, maxd]);
        
                    
    end
end




allPoints = false;
%%
if mrnaFlag
    %accumulated mrna stuff
    figure()
    axsmrna = {};
    
    for cycle = 1:1
        
        %     axsmrna{cycle} = subplot(1, 2, cycle);
        if ~allPoints
            errorbar(dlfluobins, meanallmrnasEmbryo{cycle},seallmrnasEmbryo{cycle}, '-o');
        else
            %             plot(dlfluobins, allmrnas{cycle}, 'ko'); %this line plots the
            %             mean within an embryo
            plot(dlfluobins, allmrnasnc12mat, '.k');
            hold on
            for i = 1:length(dlfluobins)
                plot(dlfluobins(i), nanmean(allmrnasnc12mat(i, :)), 'ob');
                hold on
            end
        end
        xlabel('dorsal concentration (au)');
        ylabel('mean acccumulated mRNA from active nuclei(au)');
        ylim([0, max(meanallmrnasEmbryo{cycle}*1.1)]);
        xlim(dvLim);
        title([DataType, ' nc',num2str(cycle+11)]);
        standardizeFigure(gca, []);
        set(gca, 'XScale', 'log')
    end
end

if durationFlag
    
    figure()
    cycle = 1;
    errorbar(dlfluobins, meanalldurationsEmbryo{cycle},sealldurationsEmbryo{cycle}, '-o');
    xlabel('dorsal concentration (au)');
    ylabel('duration active nuclei (frames)');
    xlim(dvLim);
    title([DataType, ' nc',num2str(cycle+11), ' duration']);
    standardizeFigure(gca, []);
    set(gca, 'XScale', 'log')
end

if turnOnFlag
    
    figure()
    cycle = 1;
    errorbar(dlfluobins, meanTurnOnsEmbryo{cycle},seTurnOnsEmbryo{cycle}, '-o');
    xlabel('dorsal concentration (au)');
    ylabel('turn on times active nuclei (min)');
    xlim(dvLim);
    title([DataType, ' nc',num2str(cycle+11), ' turn on time']);
    standardizeFigure(gca, []);
    set(gca, 'XScale', 'log')
end

if maxFlag
    
    figure()
    cycle = 1;
    errorbar(dlfluobins, meanAllMaxFluoEmbryo{cycle},seAllMaxFluoEmbryo{cycle}, '-o');
    xlabel('dorsal concentration (au)');
    ylabel('95% intensity from spots in active nuclei (au)');
    xlim(dvLim);
    title([DataType, ' nc',num2str(cycle+11), ' 95% of brightest spots']);
    standardizeFigure(gca, []);
    set(gca, 'XScale', 'log')

    
end





end
