function [NEllipsesAP, OnRatioAP, MeanVectorAllAP, SEVectorAllAP,...
    EllipsesFilteredPos, FilteredParticlesPos, ParticleCountAP, ParticleCountProbAP] =...
    ...
    computeFractionFirstHalf(Particles, ChN, CompiledParticles, Ellipses, schnitzcells,...
    APbinID, pixelSize, FrameInfo, EllipsePos, DropboxFolder, Prefix, ElapsedTime, edgeWidth, SkipAll)

NEllipsesAP = [];
OnRatioAP = [];
MeanVectorAllAP = [];
SEVectorAllAP = [];
EllipsesFilteredPos = [];
FilteredParticlesPos = [];
ParticleCountAP = [];
ParticleCountProbAP = [];
EdgeWidth= edgeWidth/pixelSize;

if ~isempty(Particles{ChN})
    
    %I'll use the Ellipses structure to count nuclei. This is because
    %schnitzcells sometimes misses things at the edges.
    
    
    %First, add the corresponding Ellipse number to each frame of the
    %particles. Also save the information in a cell array. This should make
    %searching easier.
    ParticleNuclei = {};
    ParticleFrames = {};
    for i=1:length(Particles{ChN})
        %3/29/19 JL: Also make sure to check that the nucleus label
        %isn't 0.
        if ~isempty(Particles{ChN}(i).Nucleus) && Particles{ChN}(i).Nucleus ~=0
            ParticleNuclei{i}=...
                schnitzcells(Particles{ChN}(i).Nucleus).cellno(ismember(schnitzcells(Particles{ChN}(i).Nucleus).frames,...
                Particles{ChN}(i).Frame));
            ParticleFrames{i}=...
                schnitzcells(Particles{ChN}(i).Nucleus).frames(ismember(schnitzcells(Particles{ChN}(i).Nucleus).frames,...
                Particles{ChN}(i).Frame));
        else
            ParticleNuclei{i}=[];
            ParticleFrames{i}=[];
        end
    end
    
    %Do the analogous for CompiledParticles. We'll use this to estimate
    %fluorescence per ALL nuclei
    CompiledParticleNuclei = {};
    CompiledParticleFrames = {};
    for i=1:length(CompiledParticles{ChN})
        if ~isempty(Particles{ChN}(i).Nucleus)
            CompiledParticleNuclei{i}=...
                schnitzcells(CompiledParticles{ChN}(i).Nucleus).cellno(ismember(schnitzcells(CompiledParticles{ChN}(i).Nucleus).frames,...
                CompiledParticles{ChN}(i).Frame));
            CompiledParticleFrames{i}=...
                schnitzcells(CompiledParticles{ChN}(i).Nucleus).frames(ismember(schnitzcells(CompiledParticles{ChN}(i).Nucleus).frames,...
                CompiledParticles{ChN}(i).Frame));
        else
            CompiledParticleNuclei{i}=[];
            CompiledParticleFrames{i}=[];
        end
    end
    
    
    
    %For each frame find the number of ellipses that are outside of an area
    %delimited from the edge of the image.
    %The information in Ellipses is
    %(x, y, a, b, theta, maxcontourvalue, time, particle_id)
    
    %Initialize matrices where we will store the information of number of
    %particles vs. AP vs. time
    NEllipsesAP=zeros(length(Ellipses),length(APbinID));
    NParticlesEllipsesAP{ChN}=zeros(length(Ellipses),length(APbinID));
    %Fluorescence per all of nuclei
    MeanVectorAllAP{ChN}=nan(length(Ellipses),length(APbinID));
    SEVectorAllAP{ChN}=nan(length(Ellipses),length(APbinID));
    
    
    
    
    for i=1:length(Ellipses)
        CurrentEllipses=Ellipses{i};
        
        Radius=max(CurrentEllipses(:,3:4)')';
        
        EllipseFilter=(CurrentEllipses(:,1)-Radius-EdgeWidth>0)&...
            (CurrentEllipses(:,2)-Radius-EdgeWidth>0)&...
            (CurrentEllipses(:,1)+Radius+EdgeWidth<FrameInfo(1).PixelsPerLine)&...
            (CurrentEllipses(:,2)+Radius+EdgeWidth<FrameInfo(1).LinesPerFrame);
        
        
        %Figure out which particles are in this frame and have an approved
        %flag of 1 or 2. Note that we haven't yet checked if their
        %corresponding Ellipses have been approved.
        CurrentParticlesIndex=cellfun(@(x) find(x==i),ParticleFrames,...
            'UniformOutput',false);
        CurrentParticlesFilter=~cellfun(@isempty,CurrentParticlesIndex);
        ParticlesToCheck=find(CurrentParticlesFilter);
        
        
        %Find which of the particles in this frame are related to filtered
        %ellipses and save their corresonding AP information.
        
        FilteredParticlesPos=[];
        for j=1:length(ParticlesToCheck)
            if EllipseFilter(ParticleNuclei{ParticlesToCheck(j)}(CurrentParticlesIndex{ParticlesToCheck(j)}))&...
                    (Particles{ChN}(ParticlesToCheck(j)).Approved==1 | Particles{ChN}(ParticlesToCheck(j)).Approved==2)
                FilteredParticlesPos=[FilteredParticlesPos,...
                    EllipsePos{i}(ParticleNuclei{ParticlesToCheck(j)}(CurrentParticlesIndex{ParticlesToCheck(j)}))];
            end
        end
        
        %Count the number of filtered ellipses per AP bin
        EllipsesFilteredPos{i}=EllipsePos{i}(EllipseFilter);
        for j=1:length(EllipsesFilteredPos{i})
            NEllipsesAP(i,max(find(APbinID<=EllipsesFilteredPos{i}(j))))=...
                NEllipsesAP(i,max(find(APbinID<=EllipsesFilteredPos{i}(j))))+1;
        end
        
        
        
        %Count the number of filtered particles per AP bin.
        for j=1:length(FilteredParticlesPos)
            NParticlesEllipsesAP{ChN}(i,max(find(APbinID<=FilteredParticlesPos(j))))=...
                NParticlesEllipsesAP{ChN}(i,max(find(APbinID<=FilteredParticlesPos(j))))+1;
        end
        
        
        
        
        EllipsesFiltered{ChN}{i}=Ellipses{i}(EllipseFilter,:);
        
        NEllipsesFiltered{ChN}(i)=sum(EllipseFilter);
        
        
        %Calculate the fluorescence per ellipse. Here we'll draw the fluorescence from
        %CompiledParticles just to make sure that everything has been
        %quantified correctly.
        CurrentCompiledParticlesIndex=cellfun(@(x) find(x==i),CompiledParticleFrames,...
            'UniformOutput',false);
        CurrentCompiledParticlesFilter=~cellfun(@isempty,CurrentCompiledParticlesIndex);
        CompiledParticlesToCheck=find(CurrentCompiledParticlesFilter);
        
        %Find the fluorescence of each set of particles and the AP
        %positions of their Ellipses
        FluorescenceCompiledParticles=[];
        ErrorFluorescenceCompiledParticles=[];
        FilteredCompiledParticlesPos=[];
        for j=1:length(CompiledParticlesToCheck)
            if EllipseFilter(CompiledParticleNuclei{CompiledParticlesToCheck(j)}(CurrentCompiledParticlesIndex{CompiledParticlesToCheck(j)}))&...
                    (CompiledParticles{ChN}(CompiledParticlesToCheck(j)).Approved==1 | CompiledParticles{ChN}(CompiledParticlesToCheck(j)).Approved==2)
                FilteredCompiledParticlesPos=[FilteredCompiledParticlesPos,...
                    EllipsePos{i}(CompiledParticleNuclei{CompiledParticlesToCheck(j)}(CurrentCompiledParticlesIndex{CompiledParticlesToCheck(j)}))];
                
                FluorescenceCompiledParticles=[FluorescenceCompiledParticles,...
                    CompiledParticles{ChN}(CompiledParticlesToCheck(j)).Fluo(...
                    CurrentCompiledParticlesIndex{CompiledParticlesToCheck(j)})];
                
                ErrorFluorescenceCompiledParticles=[ErrorFluorescenceCompiledParticles,...
                    CompiledParticles{ChN}(CompiledParticlesToCheck(j)).FluoError];
            end
        end
        try
            %Sum the fluorescence values and divide by the number of ellipses
            for j=1:length(APbinID)
                if ~isnan(APbinArea(j))
                    APFilterTemp=(APbinID(j)<=FilteredCompiledParticlesPos)&...
                        (FilteredCompiledParticlesPos<APbinID(j+1));
                    
                    MeanVectorAllAP{ChN}(i,j)=sum(FluorescenceCompiledParticles(APFilterTemp))/NEllipsesAP(i,j);
                    SEVectorAllAP{ChN}(i,j)=sqrt(sum((ErrorFluorescenceCompiledParticles(APFilterTemp)/NEllipsesAP(i,j)).^2));
                end
            end
        catch
            %AR 10/22- not really sure why this fails.
        end
    end
    
    
    %Minimum number of nuclei to actually do the calculation
    MinNuclei=3;
    MinAPIndexProb=min(find(sum(NEllipsesAP>=MinNuclei)));
    MaxAPIndexProb=max(find(sum(NEllipsesAP>=MinNuclei)));
    
    MinNucleiFilter=NEllipsesAP>=MinNuclei;
    
    %Calculate the ratio
    OnRatioAP{ChN}=NParticlesEllipsesAP{ChN}./NEllipsesAP;
    %Filter out the elements that correspond to a number of nuclei below our
    %limit of MinNuclei
    OnRatioAP{ChN}(~MinNucleiFilter)=nan;
    OnRatioAP{ChN}=reshape(OnRatioAP{ChN},size(MinNucleiFilter));
    
    if ~SkipAll
        if MaxAPIndexProb>MinAPIndexProb
            colormap(jet(128));
            cmap=colormap;
            
            Color=cmap(round((APbinID(MinAPIndexProb:MaxAPIndexProb)-...
                APbinID(MinAPIndexProb))/...
                (APbinID(MaxAPIndexProb)-APbinID(MinAPIndexProb))*127)+1,:);
            figure(15)
            clf
            PlotHandle=[];
            hold on
            %             for j=MinAPIndexProb:MaxAPIndexProb
            %                 PlotHandle=[PlotHandle,...
            %                     plot(ElapsedTime,OnRatioAP{ChN}(:,j),'color',Color(j-MinAPIndexProb+1,:))];
            %             end
            hold off
            xlabel('Time (min)')
            ylabel('Fraction of on nuclei')
            h = colorbar;
            caxis([APbinID(MinAPIndexProb),APbinID(MaxAPIndexProb)])
            ylabel(h,'AP Position (x/L)')
            StandardFigure(PlotHandle,gca)
            try
                xlim([0,ElapsedTime(end)])
            catch
                warning('Y2K bug strikes again? ElapsedTime field of FrameInfo is broken.');
            end
            ylim([0,1.01])
            saveas(gca,[DropboxFolder,filesep,Prefix,filesep,'Probabilities',filesep,'ProbVsTimeVsAP.tif'])
        end
    end
    
    
end
end