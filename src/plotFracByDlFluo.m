function [npartFluo, nschnitzFluo, axs, axsmrna] = plotFracByDlFluo(DataType, varargin)

fractionFlag = false;
mrnaFlag = false;
durationFlag = false;
maxFlag = false;
displayTiles = false;
turnOnFlag = false;

minNuclei = 3; %minimum nuclei for a bin to be plottable

for i = 1:length(varargin)
    if strcmpi(varargin{i}, 'fraction')
        fractionFlag = true;
    elseif strcmpi(varargin{i}, 'mrna')
        mrnaFlag = true;
    elseif strcmpi(varargin{i}, 'duration')
        durationFlag = true;
    elseif strcmpi(varargin{i}, 'maxfluo')
        maxFlag = true;
    elseif strcmpi(varargin{i}, 'turnOn')
        turnOnFlag = true;
    elseif strcmpi(varargin{i}, 'displayTiles')
        displayTiles = true;
        tileFig = figure();
        holdFig = figure();
    end
end

if ischar(DataType)
    [allData, Prefixes, resultsFolder] = LoadMS2Sets(DataType, 'noCompiledNuclei');
else
    allData = DataType;
    DataType = inputname(1);
end

load([resultsFolder,filesep,DataType,filesep,'dlfluobins.mat'], 'dlfluobins');
load([resultsFolder,filesep,Prefixes{1},filesep,'FrameInfo.mat'], 'FrameInfo')

ch = 1;
nEmbryos = length(allData);
nBins = length(dlfluobins);

npartFluo = {zeros(1, nBins), zeros(1, nBins), zeros(1, nBins)};
nschnitzFluo = {zeros(1, nBins), zeros(1, nBins), zeros(1, nBins)};
binFilter =  {ones(1, nBins), ones(1, nBins), ones(1, nBins)};
embryosPerBin = {zeros(1, nBins), zeros(1, nBins), zeros(1, nBins)};
allmrnasnomean = cell(3, nBins);

npartFluoEmbryo = {};
nschnitzFluoEmbryo = {};
fracFluoEmbryo = {};

for e = 1:nEmbryos
    
    schnitzcells = allData(e).Particles.schnitzcells;
    CompiledParticles = allData(e).Particles.CompiledParticles;
    load([resultsFolder,filesep,Prefixes{e},filesep,'FrameInfo.mat'], 'FrameInfo')

    ncFrames = [zeros(1,8), allData(e).Particles.nc9, allData(e).Particles.nc10, allData(e).Particles.nc11, allData(e).Particles.nc12, allData(e).Particles.nc13, allData(e).Particles.nc14]; 
    ncFrames(ncFrames==0) = 1;
    time = [FrameInfo.Time]/60; %frame times in minutes 
    ncTimes = time(ncFrames);
    
    for nc = 12:14
        for bin = 1:nBins
            
            tempParticlesFluo = [];
            accumulatedFluo = [];
            durations = [];
            maxFluos = [];
            turnOnTimes = [];
            
            if ~isempty(CompiledParticles{ch})
                particlesFluo = find([CompiledParticles{ch}.cycle] == nc & [CompiledParticles{ch}.dlfluobin] == bin);
            end
            schnitzesFluo = find([schnitzcells.cycle] == nc & [schnitzcells.dlfluobin] == bin...
                & [schnitzcells.Approved]);
            
            for p = 1:length(particlesFluo)
                if schnitzcells(CompiledParticles{ch}(particlesFluo(p)).schnitz).Approved
                    tempParticlesFluo = [tempParticlesFluo, particlesFluo(p)];
                    
                    fluoFrames = CompiledParticles{ch}(particlesFluo(p)).Frame;
                    fluo = CompiledParticles{ch}(particlesFluo(p)).Fluo3DRaw'; %can also use .Fluo for 2d fluo
                    
                    durations = [durations, time(max(fluoFrames)) - time(min(fluoFrames))];
                    
                    turnOnTimes = [turnOnTimes, time(min(fluoFrames)) - ncTimes(nc)];
                    
                    maxFluos = vertcat(maxFluos, fluo(fluo>=prctile(fluo,95)));
                    
                    if length(fluoFrames) > 1
                        accumulatedFluo = [accumulatedFluo, trapz(fluoFrames, fluo, 1)];
                    else
                        accumulatedFluo = [accumulatedFluo, fluo];
                    end
                    
                    if displayTiles && nc==12
%                         figure(tileFig)
%                         nexttile;
%                         %dorsal
%                         yyaxis left
%                         plot(schnitzcells(CompiledParticles{ch}(particlesFluo(p)).schnitz).frames, schnitzcells(CompiledParticles{ch}(particlesFluo(p)).schnitz).FluoTimeTrace, '-g');
%                         xticks([]);
%                         yticks([]);
%                         ylim([0, max(dlfluobins)]);
%                         yyaxis right
%                         plot(fluoFrames, fluo, '-r'); %spot intensity
%                         ylim([0, 1000]);
%                         %                         plot(midCycle, schnitzcells(s).FluoFeature, 'ob');
%                         %                         hold off
%                         xticks([]);
%                         yticks([]);
%                         %                         yticks([min(schnitzcells(s).FluoTimeTrace), max(schnitzcells(s).FluoTimeTrace)]);
                       
                        
                        
                        figure(holdFig)
                          yyaxis left
                        plot(schnitzcells(CompiledParticles{ch}(particlesFluo(p)).schnitz).frames, schnitzcells(CompiledParticles{ch}(particlesFluo(p)).schnitz).FluoTimeTrace, '-g');
                        ylim([0, max(dlfluobins)]);
                        yyaxis right
                        plot(fluoFrames, fluo, '.-r'); %spot intensity
%                         hold on
                        ylim([0, 1000]);
                        waitforbuttonpress;
                    end
                    
                end
            end
            
            particlesFluo = tempParticlesFluo;
            
            allmrnasnomean{nc-11,bin} = [allmrnasnomean{nc-11,bin}, accumulatedFluo];
            
            npartFluoEmbryo{nc-11}(bin, e) = length(particlesFluo);
            nschnitzFluoEmbryo{nc-11}(bin, e) = length(schnitzesFluo);
            allmrnasEmbryo{nc-11}(bin, e) = nanmean(accumulatedFluo);
            alldurationsEmbryo{nc-11}(bin, e) = nanmean(durations);
            allTurnOnsEmbryo{nc-11}(bin, e) = nanmean(turnOnTimes);
            allMaxFluoEmbryo{nc-11}(bin, e) = nanmean(maxFluos);
            
            if nschnitzFluoEmbryo{nc-11}(bin, e) >= minNuclei
                embryosPerBin{nc-11}(bin) = embryosPerBin{nc-11}(bin) + 1;
            end
            
        end
        
        
        fracFluoEmbryo{nc-11}(:, e) = npartFluoEmbryo{nc-11}(:,e)./nschnitzFluoEmbryo{nc-11}(:,e);
        
    end
    
end


for nc = 1:2
    
    binFilter{nc} =  double(embryosPerBin{nc} > 0);
    binFilter{nc}(~binFilter{nc}) = nan;
    
    filteredMean = @(x) nanmean(x,2).*binFilter{nc}';
    filteredSE = @(x) (nanstd(x,0, 2)./sqrt(embryosPerBin{nc}')).*binFilter{nc}';
    
    meanFracFluoEmbryo{nc} = filteredMean(fracFluoEmbryo{nc});
    seFracFluoEmbryo{nc} = filteredSE(fracFluoEmbryo{nc});
    
    meanallmrnasEmbryo{nc} = filteredMean(allmrnasEmbryo{nc});
    seallmrnasEmbryo{nc} = filteredSE(allmrnasEmbryo{nc});
    
    meanalldurationsEmbryo{nc} = filteredMean(alldurationsEmbryo{nc});
    sealldurationsEmbryo{nc} = filteredSE(alldurationsEmbryo{nc});
    
    meanTurnOnsEmbryo{nc} = filteredMean(allTurnOnsEmbryo{nc});
    seTurnOnsEmbryo{nc} = filteredSE(allTurnOnsEmbryo{nc});
    
    meanAllMaxFluoEmbryo{nc} = filteredMean(allMaxFluoEmbryo{nc});
    seAllMaxFluoEmbryo{nc} = filteredSE(allMaxFluoEmbryo{nc});
    
end

axs = {};
axsmrna = {};
%%
%plotting
if fractionFlag
    figure()
    axs = {};
    for cycle = 1:1
        
        %     axs{cycle} = subplot(1, 2, cycle);
        errorbar(dlfluobins, meanFracFluoEmbryo{cycle},seFracFluoEmbryo{cycle}, '-o');
        xlabel('dorsal concentration (au)');
        ylabel('fraction active nuclei');
        ylim([0, 1]);
        xlim([0, max(dlfluobins)*1.1]);
        title([DataType, ' nc',num2str(cycle+11)]);
        standardizeFigure(gca, []);
        
    end
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
        xlim([0, max(dlfluobins)*1.1]);
        title([DataType, ' nc',num2str(cycle+11)]);
        standardizeFigure(gca, []);
        
    end
end

if durationFlag
    
    figure()
    cycle = 1;
    errorbar(dlfluobins, meanalldurationsEmbryo{cycle},sealldurationsEmbryo{cycle}, '-o');
    xlabel('dorsal concentration (au)');
    ylabel('duration active nuclei (frames)');
    xlim([0, max(dlfluobins)*1.1]);
    title([DataType, ' nc',num2str(cycle+11), ' duration']);
    standardizeFigure(gca, []);
    
end

if turnOnFlag
    
    figure()
    cycle = 1;
    errorbar(dlfluobins, meanTurnOnsEmbryo{cycle},seTurnOnsEmbryo{cycle}, '-o');
    xlabel('dorsal concentration (au)');
    ylabel('turn on times active nuclei (min)');
    xlim([0, max(dlfluobins)*1.1]);
    title([DataType, ' nc',num2str(cycle+11), ' turn on time']);
    standardizeFigure(gca, []);
    
end

if maxFlag
    
    figure()
    cycle = 1;
    errorbar(dlfluobins, meanAllMaxFluoEmbryo{cycle},seAllMaxFluoEmbryo{cycle}, '-o');
    xlabel('dorsal concentration (au)');
    ylabel('95% intensity from spots in active nuclei (au)');
    xlim([0, max(dlfluobins)*1.1]);
    title([DataType, ' nc',num2str(cycle+11), ' 95% of brightest spots']);
    standardizeFigure(gca, []);
    
end





end
