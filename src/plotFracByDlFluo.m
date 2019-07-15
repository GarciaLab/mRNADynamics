function [npart, nschnitz, npartFluo, nschnitzFluo] = plotFracByDlFluo(DataType)

if ischar(DataType)
    [allData, Prefixes, resultsFolder] = LoadMS2Sets(DataType);
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
            particles = find([CompiledParticles{ch}.cycle] == nc & [CompiledParticles{ch}.dvbin] == bin);
            schnitzes = find([schnitzcells.cycle] == nc & [schnitzcells.dvbin] == bin);
            particlesFluo = find([CompiledParticles{ch}.cycle] == nc & [CompiledParticles{ch}.dlfluobin] == bin);
            schnitzesFluo = find([schnitzcells.cycle] == nc & [schnitzcells.dlfluobin] == bin);
            
            npart{nc-11}(bin) = npart{nc-11}(bin) + length(particles);
            nschnitz{nc-11}(bin) = nschnitz{nc-11}(bin) + length(schnitzes);
            npartFluo{nc-11}(bin) = npartFluo{nc-11}(bin) + length(particlesFluo);
            nschnitzFluo{nc-11}(bin) = nschnitzFluo{nc-11}(bin) + length(schnitzesFluo);
        end
        
        fracFluo{nc-11} = npartFluo{nc-11}./nschnitzFluo{nc-11};
    end
end


figure()
for cycle = 1:3
    subplot(1, 3, cycle)
    plot(dlfluobins, fracFluo{cycle}, 'o');
    set(gca, 'YScale', 'log');
    xlabel('dorsal concentration (au)');
    ylabel('fraction active nuclei');
    title(['nc',num2str(cycle+11)]);
    standardizeFigure(gca, []);
end


end
