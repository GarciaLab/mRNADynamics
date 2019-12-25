function [npartFluo, nschnitzFluo, axs, axsmrna] = plotFracByDlFluo(DataType, varargin)

fractionFlag = false;
mrnaFlag = false;
durationFlag = false;
maxFlag = false;

for i = 1:length(varargin)
    if strcmpi(varargin{i}, 'fraction')
        fractionFlag = true;
    elseif strcmpi(varargin{i}, 'mrna')
        mrnaFlag = true;
    elseif strcmpi(varargin{i}, 'duration')
        durationFlag = true;
    elseif strcmpi(varargin{i}, 'maxfluo')
        maxFlag = true;
    end
end

if ischar(DataType)
    [allData, Prefixes, resultsFolder] = LoadMS2Sets(DataType, 'noCompiledNuclei');
else
    allData = DataType;
    DataType = inputname(1);
end

load([resultsFolder,filesep,DataType,filesep,'dlfluobins.mat'], 'dlfluobins');

ch = 1;
nbins = length(dlfluobins);
npart = {zeros(1, nbins), zeros(1, nbins), zeros(1, nbins)};
nschnitz = {zeros(1, nbins), zeros(1, nbins), zeros(1, nbins)};

npartFluo = {zeros(1, nbins), zeros(1, nbins), zeros(1, nbins)};
nschnitzFluo = {zeros(1, nbins), zeros(1, nbins), zeros(1, nbins)};
allmrnas = {zeros(1, nbins), zeros(1, nbins), zeros(1, nbins)};
% allmrnasnomean = {zeros(1, nbins), zeros(1, nbins), zeros(1, nbins)};
allmrnasnomean = {[], [], []};

npartFluoEmbryo = {};
nschnitzFluoEmbryo = {};
fracFluoEmbryo = {};

for e = 1:length(allData)
    
    schnitzcells = allData(e).Particles.schnitzcells;
    CompiledParticles = allData(e).Particles.CompiledParticles;
    
    for nc = 12:14
        for bin = 1:nbins

            if ~isempty(CompiledParticles{ch})
                particlesFluo = find([CompiledParticles{ch}.cycle] == nc & [CompiledParticles{ch}.dlfluobin] == bin);
            end
            schnitzesFluo = find([schnitzcells.cycle] == nc & [schnitzcells.dlfluobin] == bin...
                & [schnitzcells.Approved]);

            tempParticlesFluo = [];
            particlesFluoCopy = particlesFluo;
            accumulatedmRNA = [];
            durations = [];
            maxFluos = [];
            
            for p = 1:length(particlesFluo)
                if schnitzcells(CompiledParticles{ch}(particlesFluo(p)).schnitz).Approved
                    tempParticlesFluo = [tempParticlesFluo, particlesFluo(p)];
                    
                    fluoFrames = CompiledParticles{ch}(particlesFluo(p)).Frame;
                    fluo = CompiledParticles{ch}(particlesFluo(p)).Fluo3DRaw'; %can also use .Fluo for 2d fluo
                    
                    durations = [durations, length(fluoFrames)];
                  
                    maxFluos = vertcat(maxFluos, fluo(fluo>=prctile(fluo,0)));
                    
                    if length(fluoFrames) > 1
                        accumulatedmRNA = [accumulatedmRNA, trapz(fluoFrames, fluo, 1)];
                    else
                        accumulatedmRNA = [accumulatedmRNA, fluo];
                    end
                    
                    
                end
            end
            
            particlesFluo = tempParticlesFluo;
            
            npartFluo{nc-11}(bin) = npartFluo{nc-11}(bin) + length(particlesFluo);
            nschnitzFluo{nc-11}(bin) = nschnitzFluo{nc-11}(bin) + length(schnitzesFluo);
            allmrnas{nc-11}(bin) = allmrnas{nc-11}(bin) + mean(accumulatedmRNA);
          
          
             
%              allmrnasnomean{nc-11,bin} = [allmrnasnomean{nc-11,bin}, mrnas];
            
            npartFluoEmbryo{nc-11}(bin, e) = length(particlesFluo);
            nschnitzFluoEmbryo{nc-11}(bin, e) = length(schnitzesFluo);
            allmrnasEmbryo{nc-11}(bin, e) = mean(accumulatedmRNA);
            alldurationsEmbryo{nc-11}(bin, e) = mean(durations);
            allMaxFluoEmbryo{nc-11}(bin, e) = mean(maxFluos);
            
        end
        
        fracFluo{nc-11} = npartFluo{nc-11}./nschnitzFluo{nc-11};
        
        fracFluoEmbryo{nc-11}(:, e) = npartFluoEmbryo{nc-11}(:,e)./nschnitzFluoEmbryo{nc-11}(:,e);
        
    end
    
end


for nc = 1:2
    
    meanFracFluoEmbryo{nc} = nanmean(fracFluoEmbryo{nc}, 2);
    nEmbryos = size(fracFluoEmbryo{nc}, 2);
    seFracFluoEmbryo{nc} = nanstd(fracFluoEmbryo{nc},0, 2)./sqrt(nEmbryos);
    
    meanallmrnasEmbryo{nc} = nanmean(allmrnasEmbryo{nc}, 2);
    seallmrnasEmbryo{nc} = nanstd(allmrnasEmbryo{nc},0, 2)./sqrt(nEmbryos);
    
    meanalldurationsEmbryo{nc} = nanmean(alldurationsEmbryo{nc}, 2);
     sealldurationsEmbryo{nc} = nanstd(alldurationsEmbryo{nc},0, 2)./sqrt(nEmbryos);
     
     meanAllMaxFluoEmbryo{nc} = nanmean(allMaxFluoEmbryo{nc}, 2);
     seAllMaxFluoEmbryo{nc} = nanstd(allMaxFluoEmbryo{nc},0, 2)./sqrt(nEmbryos);

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
% 
% allmrnasnc12 = allmrnasnomean(1, :);
% lens = [];
% for b = 1:nbins
%     lens(b) = length(allmrnasnc12{b});
% end
% allmrnasnc12mat = zeros(nbins, max(lens));
% for bin = 1:nbins
%     if ~isempty(allmrnasnc12{bin})
%         allmrnasnc12mat(bin, :) = padarray(allmrnasnc12{bin}',max(lens)-length(allmrnasnc12{bin}),NaN, 'post');
%     end
% end

%stacked bar

allPoints = true;
%%
figure();
if mrnaFlag
    %accumulated mrna stuff
    % figure()
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

        cycle = 1;
        errorbar(dlfluobins, meanalldurationsEmbryo{cycle},sealldurationsEmbryo{cycle}, '-o');
        xlabel('dorsal concentration (au)');
        ylabel('duration active nuclei (frames)');
        xlim([0, max(dlfluobins)*1.1]);
        title([DataType, ' nc',num2str(cycle+11), ' duration']);
        standardizeFigure(gca, []);

end

if maxFlag 

        cycle = 1;
        errorbar(dlfluobins, meanAllMaxFluoEmbryo{cycle},seAllMaxFluoEmbryo{cycle}, '-o');
        xlabel('dorsal concentration (au)');
        ylabel('95% intensity from spots in active nuclei (au)');
        xlim([0, max(dlfluobins)*1.1]);
        title([DataType, ' nc',num2str(cycle+11), ' 95% of brightest spots']);
        standardizeFigure(gca, []);

end





end
