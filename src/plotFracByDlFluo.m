function [npart, nschnitz, npartFluo, nschnitzFluo] = plotFracByDlFluo(DataType)

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

for e = 1:length(allData);
    
    schnitzcells = allData(e).Particles.schnitzcells;
    CompiledParticles = allData(e).Particles.CompiledParticles;
        
    for nc = 12:14
        for bin = 1:nbins
            
%             particles = find([CompiledParticles{ch}.cycle] == nc & [CompiledParticles{ch}.dvbin] == bin);
%             schnitzes = find([schnitzcells.cycle] == nc & [schnitzcells.dvbin] == bin...
%                 & [schnitzcells.Approved]);
            particlesFluo = find([CompiledParticles{ch}.cycle] == nc & [CompiledParticles{ch}.dlfluobin] == bin);
            schnitzesFluo = find([schnitzcells.cycle] == nc & [schnitzcells.dlfluobin] == bin...
                & [schnitzcells.Approved]);
%             for p = 1:length(particles)
%                 keyboard
%                 if ~schnitzcells(CompiledParticles{ch}(p).schnitz).Approved
%                     keyboard
%                     particles(p) = [];
%                 end
%             end
            tempParticlesFluo = [];
            for p = 1:length(particlesFluo)
                if schnitzcells(CompiledParticles{ch}(p).schnitz).Approved
                    tempParticlesFluo = [tempParticlesFluo, particlesFluo(p)];
                end
            end
            
%             npart{nc-11}(bin) = npart{nc-11}(bin) + length(particles);
%             nschnitz{nc-11}(bin) = nschnitz{nc-11}(bin) + length(schnitzes);
            npartFluo{nc-11}(bin) = npartFluo{nc-11}(bin) + length(particlesFluo);
            nschnitzFluo{nc-11}(bin) = nschnitzFluo{nc-11}(bin) + length(schnitzesFluo);
        end
        
        fracFluo{nc-11} = npartFluo{nc-11}./nschnitzFluo{nc-11};
    end
end


figure()
for cycle = 1:3
    subplot(1, 3, cycle)
    plot(dlfluobins, fracFluo{cycle}, '-o');
    xlabel('dorsal concentration (au)');
    ylabel('fraction active nuclei');
    ylim([0, 1]);
    xlim([0, max(dlfluobins)*1.1]);
    title(['nc',num2str(cycle+11)]);
    standardizeFigure(gca, []);
end


%stacked bar
figure()
for cycle = 1:3
    subplot(1, 3, cycle)
    bardata = cat(1,npartFluo{cycle}, nschnitzFluo{cycle});
    bar(dlfluobins, bardata', 'stacked');
    legend('particles', 'ellipses')
    xlabel('dorsal concentration (au)');
    ylabel('number');
    title(['nc',num2str(cycle+11)]);
    standardizeFigure(gca, []);
end


end
