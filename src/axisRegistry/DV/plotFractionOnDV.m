function plotFractionOnDV(combinedEmbryo)
%Armando 2019 Edition

close all;

if ischar(combinedEmbryo)
    load(combinedEmbryo);
end

nBins = size(NParticlesDV, 2);
nFrames = length(ElapsedTime);
npDV = zeros(1, nBins);
neDV = zeros(1, nBins);
bins = 1:nBins;
ncs = [nc12, nc13, nc14, nFrames];

NEllipsesDV(NEllipsesDV==inf) = NaN;
NEllipsesDV(NEllipsesDV==0) = NaN;



figure()
for nc = 1:3
    
    for bin = bins
        npDV(bin) = nansum(NParticlesDV(ncs(nc):ncs(nc+1),bin));
        neDV(bin) = nanmax(NEllipsesDV(ncs(nc):ncs(nc+1),bin));
    end
    
    neDV(neDV==inf) = NaN;
    neDV(neDV==0) = NaN;
    
    dvAxis = bins/nBins;
    
    fractionOn = npDV./neDV;
    subplot(1,3, nc)
    plot(dvAxis, fractionOn);
    xlabel('fraction DV');
    ylabel('fraction active nuclei');
    xlim([0, 1]);
    title(['nc',num2str(nc+11)]);
    standardizeFigure(gca, []);
    
end

end