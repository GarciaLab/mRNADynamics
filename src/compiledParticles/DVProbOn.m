function [NEllipsesDV, MeanVectorAllDV, SEVectorAllDV, OnRatioDV, ParticleCountDV, ...
    ParticleCountProbDV, TotalEllipsesDV, EllipsesOnDV, EllipsesFilteredPos, ...
    FilteredParticlesPos] = DVProbOn(NChannels, ...
    Particles, schnitzcells, CompiledParticles, Ellipses, FrameInfo, ...
    DropboxFolder, Prefix, ElapsedTime, DVbinID)
%DVPROBON Summary of this function goes here
%   Detailed explanation goes here

EdgeWidth=10; %AR 1/12/18 why 10?
    
for ChN=1:NChannels

    if ~isempty(Particles{ChN})

        %I'll use the Ellipses structure to count nuclei. This is because
        %schnitzcells sometimes misses things at the edges.


        %First, add the corresponding Ellipse number to each frame of the
        %particles. Also save the information in a cell array. This should make
        %searching easier.
        clear ParticleNuclei
        clear ParticleFrames
        for i=1:length(Particles{ChN})
            if ~isempty(Particles{ChN}(i).Nucleus)
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
        clear CompiledParticleNuclei
        clear CompiledParticleFrames
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
        %particles vs. DV vs. time
        NEllipsesDV=zeros(length(Ellipses),length(DVbinID));
        NParticlesEllipsesDV{ChN}=zeros(length(Ellipses),length(DVbinID));
        %Fluorescence per all of nuclei
        MeanVectorAllDV{ChN}=nan(length(Ellipses),length(DVbinID));
        SEVectorAllDV{ChN}=nan(length(Ellipses),length(DVbinID));


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
                        EllipsePos_DV{i}(ParticleNuclei{ParticlesToCheck(j)}(CurrentParticlesIndex{ParticlesToCheck(j)}))];
                end
            end

            %Count the number of filtered ellipses per DV bin
            EllipsesFilteredPos{i}=EllipsePos_DV{i}(EllipseFilter);
            for j=1:length(EllipsesFilteredPos{i})
                NEllipsesDV(i,max(find(DVbinID<=EllipsesFilteredPos{i}(j))))=...
                    NEllipsesDV(i,max(find(DVbinID<=EllipsesFilteredPos{i}(j))))+1;
            end


            %Count the number of filtered particles per DV bin.
            for j=1:length(FilteredParticlesPos)
                NParticlesEllipsesDV{ChN}(i,max(find(DVbinID<=FilteredParticlesPos(j))))=...
                    NParticlesEllipsesDV{ChN}(i,max(find(DVbinID<=FilteredParticlesPos(j))))+1;
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

            %Find the fluorescence of each set of particles and the DV
            %positions of their Ellipses
            FluorescenceCompiledParticles=[];
            ErrorFluorescenceCompiledParticles=[];
            FilteredCompiledParticlesPos=[];
            for j=1:length(CompiledParticlesToCheck)
                if EllipseFilter(CompiledParticleNuclei{CompiledParticlesToCheck(j)}(CurrentCompiledParticlesIndex{CompiledParticlesToCheck(j)}))&...
                        (CompiledParticles{ChN}(CompiledParticlesToCheck(j)).Approved==1 | CompiledParticles{ChN}(CompiledParticlesToCheck(j)).Approved==2)
                    FilteredCompiledParticlesPos=[FilteredCompiledParticlesPos,...
                        EllipsePos_DV{i}(CompiledParticleNuclei{CompiledParticlesToCheck(j)}(CurrentCompiledParticlesIndex{CompiledParticlesToCheck(j)}))];

                    FluorescenceCompiledParticles=[FluorescenceCompiledParticles,...
                        CompiledParticles{ChN}(CompiledParticlesToCheck(j)).Fluo(...
                            CurrentCompiledParticlesIndex{CompiledParticlesToCheck(j)})];

                    ErrorFluorescenceCompiledParticles=[ErrorFluorescenceCompiledParticles,...
                        CompiledParticles{ChN}(CompiledParticlesToCheck(j)).FluoError];
                end
            end



            %Sum the fluorescence values and divide by the number of ellipses
            for j=1:length(DVbinID)-1
                if ~isnan(DVbinArea(j))
                    DVFilterTemp=(DVbinID(j)<=FilteredCompiledParticlesPos)&...
                        (FilteredCompiledParticlesPos<DVbinID(j+1));

                    MeanVectorAllDV{ChN}(i,j)=sum(FluorescenceCompiledParticles(DVFilterTemp))/NEllipsesDV(i,j);
                    SEVectorAllDV{ChN}(i,j)=sqrt(sum((ErrorFluorescenceCompiledParticles(DVFilterTemp)/NEllipsesDV(i,j)).^2));
                end
            end
        end


        %Minimum number of nuclei to actually do the calculation
        MinNuclei=3;
        MinDVIndexProb=min(find(sum(NEllipsesDV>=MinNuclei)));
        MaxDVIndexProb=max(find(sum(NEllipsesDV>=MinNuclei)));

        MinNucleiFilter=NEllipsesDV>=MinNuclei;

        %Calculate the ratio
        OnRatioDV{ChN}=NParticlesEllipsesDV{ChN}./NEllipsesDV;
        %Filter out the elements that correspond to a number of nuclei below our
        %limit of MinNuclei
        OnRatioDV{ChN}(~MinNucleiFilter)=nan;
        OnRatioDV{ChN}=reshape(OnRatioDV{ChN},size(MinNucleiFilter));


        if MaxDVIndexProb>MinDVIndexProb
            colormap(jet(128));
            cmap=colormap;

            Color=cmap(round((DVbinID(MinDVIndexProb:MaxDVIndexProb)-...
                DVbinID(MinDVIndexProb))/...
                (DVbinID(MaxDVIndexProb)-DVbinID(MinDVIndexProb))*127)+1,:);
            figure(15)
            clf
            PlotHandle=[];
            hold on
%             for j=MinDVIndexProb:MaxDVIndexProb
%                 PlotHandle=[PlotHandle,...
%                     plot(ElapsedTime,OnRatioDV{ChN}(:,j),'color',Color(j-MinDVIndexProb+1,:))];
%             end
            hold off
            xlabel('Time (min)')
            ylabel('Fraction of on nuclei')
            h = colorbar;
            caxis([DVbinID(MinDVIndexProb),DVbinID(MaxDVIndexProb)])
            ylabel(h,'DV Position (x/L)')
            StandardFigure(PlotHandle,gca)
            xlim([0,ElapsedTime(end)])
            ylim([0,1.01])
            saveas(gca,[DropboxFolder,filesep,Prefix,filesep,'Probabilities',filesep,'ProbVsTimeVsDV.tif'])
        end



        %Now I want to compute the probability of nuclei being on in at least one
        %frame over the whole nc. This is a little bit tricky because I haven't
        %checked the tracking of the nuclei without particles. As a result, those
        %lineages are not complete and we could possibly overcount off nuclei if we
        %just looked at the schnitzcells we have right now. I'll fix this later,
        %but for now I'm going calculate the number of on nuclei per DV bin, where
        %I'll actually check how much area corresponds to that DV bin.





        %What's the probability of a nucleus being on for at least one frame in an nc?
        %This will be an array with rows for nc13 and nc14 and columns
        %corresponding to each DV bin.
        %This only works if we trust the tracking within one nc
        ParticleCountDV{ChN}=zeros(3,length(DVbinID));

        try
            for i=1:length(Particles{ChN})
                %See if the particle has either the flag 1 or 2
                if (Particles{ChN}(i).Approved==1)||(Particles{ChN}(i).Approved==2)

                    %Determine the nc so we can add to the right position of ParticleCountDV
                    if (FrameInfo(min(Particles{ChN}(i).Frame(Particles{ChN}(i).FrameApproved))).nc)>=12
                        CurrentNC=FrameInfo(min(Particles{ChN}(i).Frame(Particles{ChN}(i).FrameApproved))).nc;

                        %Now determine the DV bin this particle is on. We'll use
                        %the nucleus positioning. Also, in order to count it we
                        %need to make sure that it fulfills the same criteria we
                        %used to count NEllipsesDV
                        EllipseDVPosTemp=[];
                        EllipsePosXTemp=[];
                        EllipsePosYTemp=[];
                        RadiusTemp=[];
                        xPosTemp=[];
                        yPosTemp=[];
                        for j=1:length(schnitzcells(Particles{ChN}(i).Nucleus).frames)
                            CurrentEllipse=Ellipses{schnitzcells(Particles{ChN}(i).Nucleus).frames(j)}(schnitzcells(Particles{ChN}(i).Nucleus).cellno(j),:);
                            RadiusTemp=[RadiusTemp,max(CurrentEllipse(3:4))];

                            %Determine the DV position
                            EllipseDVPosTemp=[EllipseDVPosTemp,...
                                EllipsePos_DV{schnitzcells(Particles{ChN}(i).Nucleus).frames(j)}(schnitzcells(Particles{ChN}(i).Nucleus).cellno(j))];

                            %Determine the x and y positions
                            xPosTemp=[xPosTemp,CurrentEllipse(1)];
                            yPosTemp=[yPosTemp,CurrentEllipse(2)];
                        end

                        %Make sure that this nucleus was inside the limits at all
                        %points

                        if sum((xPosTemp-RadiusTemp-EdgeWidth>0)&...
                            (yPosTemp-RadiusTemp-EdgeWidth>0)&...
                            (xPosTemp+RadiusTemp+EdgeWidth<FrameInfo(1).PixelsPerLine)&...
                            (yPosTemp+RadiusTemp+EdgeWidth<FrameInfo(1).LinesPerFrame))==...
                            length(xPosTemp)

                            MeanDV=mean(EllipseDVPosTemp);

                            ParticleCountDV{ChN}(CurrentNC-11,max(find(DVbinID<=MeanDV)))=...
                                ParticleCountDV{ChN}(CurrentNC-11,max(find(DVbinID<=MeanDV)))+1;
                        end
                    end
                end
            end
        end

        %Calculate the probability using the mean number of nuclei per DV bin
        %in each nc. Note that we look at a reduced range within the nc to
        %reduce variability in counting at mitosis.
        ParticleCountProbDV{ChN}(:,1)=ParticleCountDV{ChN}(1,:)./mean(NEllipsesDV(nc12+5:nc13-5,:));
        ParticleCountProbDV{ChN}(:,2)=ParticleCountDV{ChN}(2,:)./mean(NEllipsesDV(nc13+5:nc14-5,:));
        ParticleCountProbDV{ChN}(:,3)=ParticleCountDV{ChN}(3,:)./...
            mean(NEllipsesDV(max(1,nc14-5):numFrames-5,:));
        % ES 2014-01-08: accounting for movies started fewer than 5 frames before
        % mitosis 13

        figure(16)   
        plot(DVbinID,ParticleCountProbDV{ChN}(:,1),'.-b')
        hold on
        plot(DVbinID,ParticleCountProbDV{ChN}(:,2),'.-k')
        plot(DVbinID,ParticleCountProbDV{ChN}(:,3),'.-r')
        hold off
        %ylim([0,max(ParticleCountDV(1,:))*2*1.1])
        xlabel('DV position (x)')
        ylabel('Active nuclei')
        title('Number active nuclei')
        legend('nc12', 'nc13', 'nc14')
        saveas(gca,[DropboxFolder,filesep,Prefix,filesep,'Probabilities',filesep,'ProbVsDV_ch',...
            iIndex(ChN,2),'.tif'])

    end


    %Use the alternative approach I used for the movies. We are going to
    %look at each nucleus towards the end of each nc and ask if they
    %correspond to an on or off particle in any frame previous to that one.


    %We'll go for 2.5 minutes before the next mitosis. I might relate this
    %later to the elongation time as a way to say that these particles
    %won't contribute to the total amount of mRNA produced anyway.
    FramesBack=ceil(2.5/mean(diff(ElapsedTime)));

    TotalEllipsesDV=zeros(length(DVbinID),3);
    EllipsesOnDV{ChN}=zeros(length(DVbinID),3);
    for nc=12:14

        %Figure out which frame we'll look at
        if nc==14
            FrameToUse=numFrames-FramesBack;
            %JAKE
            %FrameToUse = nc14+41;
        else
            FrameToUse=eval(['nc',num2str(nc+1)])-FramesBack;
        end

        if FrameToUse>0
            %Filter ellipses that are within the image
            CurrentEllipses=Ellipses{FrameToUse};

            Radius=max(CurrentEllipses(:,3:4)')';

            EllipseFilter=(CurrentEllipses(:,1)-Radius-EdgeWidth>0)&...
                (CurrentEllipses(:,2)-Radius-EdgeWidth>0)&...
                (CurrentEllipses(:,1)+Radius+EdgeWidth<FrameInfo(1).PixelsPerLine)&...
                (CurrentEllipses(:,2)+Radius+EdgeWidth<FrameInfo(1).LinesPerFrame);

            %Check if the filtered ellipses had an associated particle
            EllipsesToCheck=find(EllipseFilter);

            for j=1:length(EllipsesToCheck)
                %Find which DV bind we're in
                CurrentDVbin=max(find(DVbinID<EllipsePos_DV{FrameToUse}(EllipsesToCheck(j))));
                %Count the total amount of ellipses in the right DV bin
                TotalEllipsesDV(CurrentDVbin,nc-11)=TotalEllipsesDV(CurrentDVbin,nc-11)+1;



                %Find the schnitz this corresponds to
                for k=1:length(schnitzcells)

                    IndexToUse=find((schnitzcells(k).frames)==FrameToUse);
                    if ~isempty(IndexToUse)

                        %Check this schnitz for consistency with cellno.
                        %Otherwise fix it. I obtained the fixing code from
                        %TrackmRNADynamicsV2.m
                        if length(schnitzcells(k).frames)~=length(schnitzcells(k).cellno)
                           %If there number of frames is different from the number of
                           %cellno then use the cenx and ceny to find the cellno in
                           %Ellipses an repopulate this schnitz
                           if (length(schnitzcells(k).frames)==length(schnitzcells(k).cenx))&...
                                   (length(schnitzcells(k).frames)==length(schnitzcells(k).ceny))
                               for m=1:length(schnitzcells(k).frames)
                                    %The information in Ellipses is
                                    %(x, y, a, b, theta, maxcontourvalue, time, particle_id)
                                    MaxDistance=2;  %Maximum pixel distance to identify an
                                                    %ellipse with a schnitz
                                    Distances=sqrt((Ellipses{schnitzcells(k).frames(m)}(:,1)-...
                                       schnitzcells(k).cenx(m)).^2+...
                                       (Ellipses{schnitzcells(k).frames(m)}(:,2)-...
                                       schnitzcells(k).ceny(m)).^2);
                                    [MinValue,MinIndex]=min(Distances);

                                   %Make sure no other schnitz is associated to this
                                   %ellipse
                                   EllipseFoundElsewhere=0;
                                   for n=[1:k-1,k+1:length(schnitzcells)]
                                       %Only consider it if the schnitzcell is also valid!
                                       if (length(schnitzcells(n).frames)==length(schnitzcells(n).cellno))
                                           if sum(schnitzcells(n).frames==schnitzcells(k).frames(m))
                                              IndexToCheck=find(schnitzcells(n).frames==schnitzcells(k).frames(m));

                                              %The schnitz I'm comparing to
                                              %might also be screwed up.
                                              %I'd have to compare its cenx
                                              %and ceny to be sure
                                              try
                                                  if schnitzcells(k).cellno(IndexToCheck)==MinIndex
                                                        error('duplicated schnitz?')
                                                        DistancesK=sqrt((Ellipses{schnitzcells(k).frames(m)}(:,1)-...
                                                           schnitzcells(n).cenx(IndexToCheck)).^2+...
                                                           (Ellipses{schnitzcells(k).frames(m)}(:,2)-...
                                                           schnitzcells(n).ceny(IndexToCheck)).^2);          

                                                        [MinValueK,MinIndexK]=min(DistancesK);
                                                        if MinValue<MinValueK
                                                            schnitzcells(n).cellno(IndexToCheck)=[];                                        
                                                        end
                                                    end
                                              end
                                           end
                                       end
                                    end

                                    if ~EllipseFoundElsewhere
                                       schnitzcells(k).cellno(m)=MinIndex;
                                    else
                                       MinValue
                                       1+1; error('What to do here?')
                                    end
                               end

                           else
                               error('Cannnot rescue schnitz')
                           end
                        end

                    end

                    if schnitzcells(k).cellno(IndexToUse)==EllipsesToCheck(j)
                        %Now see if there is an associated particle with it
                        for m=1:length(CompiledParticles{ChN})
                            if CompiledParticles{ChN}(m).Nucleus==k
                                [j,k,m];
                                EllipsesOnDV{ChN}(CurrentDVbin,nc-11)=EllipsesOnDV{ChN}(CurrentDVbin,nc-11)+1;
                            end                        
                        end
                    end
               end
            end
        end
    end
    %{
    %JAKE: Eliminate cells with too less ellipse
    for i=1:length(TotalEllipsesDV)
        for j = 1:3
            if TotalEllipsesDV(i,j)<=3
                TotalEllipsesDV(i,j) = nan;
            end
        end
    end
    %}

    %JAKE: Quick fix, sometimes ratio>1
    for i=1:length(TotalEllipsesDV)
        for j = 1:3
            if TotalEllipsesDV(i,j)<=EllipsesOnDV{ChN}(i,j)
                EllipsesOnDV{ChN}(i,j) = TotalEllipsesDV(i,j);
            end
        end
    end

    figure(50)
    %plot(DVbinID,EllipsesOnDV{ChN}(:,1)./TotalEllipsesDV(:,1),'.-b') % fraction on nc 12
    hold on
    %plot(DVbinID,EllipsesOnDV{ChN}(:,2)./TotalEllipsesDV(:,2),'.-k') % fraction on nc 13
    plot(DVbinID,EllipsesOnDV{ChN}(:,3)./TotalEllipsesDV(:,3),'.-r') % fraction on nc 14
    hold off
    title('Fraction active nuclei')
    xlabel('DV (x)')
    ylabel('Fraction')
    axis([-700 0 0 1.05])
    %legend('nc12', 'nc13', 'nc14')
    %legend('nc13', 'nc14')
    legend('nc14')
    saveas(gca,[DropboxFolder,filesep,Prefix,filesep,'Probabilities',filesep,'FracVsDV_ch',...
            iIndex(ChN,2),'.tif'])
end
end

