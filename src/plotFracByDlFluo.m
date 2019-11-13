function [npartFluo, nschnitzFluo, axs, axsmrna] = plotFracByDlFluo(DataType)



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


npartFluoEmbryo = {};
nschnitzFluoEmbryo = {};
fracFluoEmbryo = {};

for e = 1:length(allData)
    
    schnitzcells = allData(e).Particles.schnitzcells;
    CompiledParticles = allData(e).Particles.CompiledParticles;
    
    for nc = 12:14
        for bin = 1:nbins
            
            %             particles = find([CompiledParticles{ch}.cycle] == nc & [CompiledParticles{ch}.dvbin] == bin);
            %             schnitzes = find([schnitzcells.cycle] == nc & [schnitzcells.dvbin] == bin...
            %                 & [schnitzcells.Approved]);
            if ~isempty(CompiledParticles{ch})
                particlesFluo = find([CompiledParticles{ch}.cycle] == nc & [CompiledParticles{ch}.dlfluobin] == bin);
            end
            schnitzesFluo = find([schnitzcells.cycle] == nc & [schnitzcells.dlfluobin] == bin...
                & [schnitzcells.Approved]);
            %             for p = 1:length(particles)
            %                 if ~schnitzcells(CompiledParticles{ch}(particles(p)).schnitz).Approved
            %                     particles(p) = [];
            %                 end
            %             end
            
            tempParticlesFluo = [];
            particlesFluoCopy = particlesFluo;
            mrnas = [];
            for p = 1:length(particlesFluo)
                if schnitzcells(CompiledParticles{ch}(particlesFluo(p)).schnitz).Approved
                    tempParticlesFluo = [tempParticlesFluo, particlesFluo(p)];
                    
                    mrnaframes = CompiledParticles{ch}(particlesFluo(p)).Frame;
                    mrna = CompiledParticles{ch}(particlesFluo(p)).Fluo';
                    if length(mrnaframes) > 1
                        mrnas = [mrnas, trapz(mrnaframes, mrna, 1)];
                    else
                        mrnas = [mrnas, mrna];
                    end
                end
            end
            
            particlesFluo = tempParticlesFluo;
            
            
            %             npart{nc-11}(bin) = npart{nc-11}(bin) + length(particles);
            %             nschnitz{nc-11}(bin) = nschnitz{nc-11}(bin) + length(schnitzes);
            npartFluo{nc-11}(bin) = npartFluo{nc-11}(bin) + length(particlesFluo);
            nschnitzFluo{nc-11}(bin) = nschnitzFluo{nc-11}(bin) + length(schnitzesFluo);
            allmrnas{nc-11}(bin) = allmrnas{nc-11}(bin) + mean(mrnas);
            
            npartFluoEmbryo{nc-11}(bin, e) = length(particlesFluo);
            nschnitzFluoEmbryo{nc-11}(bin, e) = length(schnitzesFluo);
            allmrnasEmbryo{nc-11}(bin, e) = mean(mrnas);
            
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
    
end


%%
%plotting

figure()
axs = {};
for cycle = 1:2
    
    axs{cycle} = subplot(1, 2, cycle);
    errorbar(dlfluobins, meanFracFluoEmbryo{cycle},seFracFluoEmbryo{cycle}, '-o');
    xlabel('dorsal concentration (au)');
    ylabel('fraction active nuclei');
    ylim([0, 1]);
    xlim([0, max(dlfluobins)*1.1]);
    title(['nc',num2str(cycle+11)]);
    standardizeFigure(gca, []);
    
end


%stacked bar
figure();
for cycle = 1:2
    subplot(1, 2, cycle)
    bardata = cat(1,npartFluo{cycle}, nschnitzFluo{cycle});
    bar(dlfluobins, bardata', 'stacked');
    leg = legend('particles', 'ellipses');
    xlabel('dorsal concentration (au)');
    ylabel('number');
    title(['nc',num2str(cycle+11)]);
    standardizeFigure(gca, leg);
end

%%
%accumulated mrna stuff
figure()
axsmrna = {};

for cycle = 1:2
    
    axsmrna{cycle} = subplot(1, 2, cycle);
    errorbar(dlfluobins, meanallmrnasEmbryo{cycle},seallmrnasEmbryo{cycle}, '-o');
    xlabel('dorsal concentration (au)');
    ylabel('mean acccumulated mRNA (au)');
    ylim([0, max(meanallmrnasEmbryo{cycle}*1.1)]);
    xlim([0, max(dlfluobins)*1.1]);
    title(['nc',num2str(cycle+11)]);
    standardizeFigure(gca, []);
    
end






end
