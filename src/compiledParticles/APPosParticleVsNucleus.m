function CompiledParticles = APPosParticleVsNucleus(NChannels, ...
    CompiledParticles, schnitzcells, EllipsePos, DropboxFolder, Prefix, SkipAll)
%APPOSPARTICLEVSNUCLEUS Summary of this function goes here
%   Detailed explanation goes here

for ChN=1:NChannels
    %How different are the AP positions of the nuclei to the particles as a
    %function of time? Let's save the information about nuclear position in the
    %CompiledParticles structure.
    for i=1:length(CompiledParticles{ChN})
        for j=1:length(CompiledParticles{ChN}(i).Frame)
            CurrentNucleus=CompiledParticles{ChN}(i).Nucleus;
            CurrentFrame=CompiledParticles{ChN}(i).Frame(j);
            
            CurrentEllipse=schnitzcells(CurrentNucleus).cellno(...
                find((schnitzcells(CurrentNucleus).frames)==CurrentFrame));
            
            if ~isempty(CurrentEllipse)
                CompiledParticles{ChN}(i).NuclearAP(j)=EllipsePos{CurrentFrame}(CurrentEllipse);
            else
                CompiledParticles{ChN}(i).NuclearAP(j)=nan;
            end
        end
    end
    
    if ~SkipAll
        figure(14)
        clf
        hold all
        for i=1:length(CompiledParticles{ChN})
            plot(CompiledParticles{ChN}(i).APpos-CompiledParticles{ChN}(i).NuclearAP)
        end
        hold off
        box on
        xlabel('Frame')
        ylabel('AP difference between spot and nucleus (x/L)')
        saveas(gca,[DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'APNucVsParticle_ch',...
            iIndex(ChN,2),'.tif'])
    end
end

end

