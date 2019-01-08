function [NEllipsesAP, MeanVectorAllAP, SEVectorAllAP, EllipsesFilteredPos, ...
    FilteredParticlesPos, OnRatioAP, ParticleCountAP, ParticleCountProbAP, ...
    EllipsesOnAP, rateOnAP, rateOnAPCell, timeOnOnAP, timeOnOnAPCell, TotalEllipsesAP]...
    = computeAPFractionOn(NChannels, Particles, schnitzcells, ...
    CompiledParticles, Ellipses, APbinID, FrameInfo, ElapsedTime, DropboxFolder, ...
    Prefix, EllipsePos, nc12, nc13, nc14, numFrames, doSingleFits, SkipAll, APbinArea, pixelSize)

%computeAPFractionOn Calculate the fraction of transcribing nuclei using three
%different methods. 
%
%HGNOTE: I'm going to measure the probability of a nucleus having detectable
%expression as a function of time and AP. In order to do this I'll use
%Particles that have both the Approved flag set to 1 and 2. However, I'll
%also check that the nuclei are not too close to the edges.
%
%HGNOTE: I need a way to go back and check the nuclei that weren't on. Maybe
%I should move this to Check particles
%
%TO DO: place each method of calculation inside a subfunction that will go
%in this file for easier reading and manipulation. Also combine the DV
%version with this version. 


    EdgeWidth=2.12/pixelSize; %in microns. 2.12 is simply the number that results in 10 pixel width
    %when using a spatial resolution of 212nm. 

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
                xlim([0,ElapsedTime(end)])
                ylim([0,1.01])
                saveas(gca,[DropboxFolder,filesep,Prefix,filesep,'Probabilities',filesep,'ProbVsTimeVsAP.tif'])
            end




            %Now I want to compute the probability of nuclei being on in at least one
            %frame over the whole nc. This is a little bit tricky because I haven't
            %checked the tracking of the nuclei without particles. As a result, those
            %lineages are not complete and we could possibly overcount off nuclei if we
            %just looked at the schnitzcells we have right now. I'll fix this later,
            %but for now I'm going calculate the number of on nuclei per AP bin, where
            %I'll actually check how much area corresponds to that AP bin.





            %What's the probability of a nucleus being on for at least one frame in an nc?
            %This will be an array with rows for nc13 and nc14 and columns
            %corresponding to each AP bin.
            %This only works if we trust the tracking within one nc
            ParticleCountAP{ChN}=zeros(3,length(APbinID));

            try
                for i=1:length(Particles{ChN})
                    %See if the particle has either the flag 1 or 2
                    if (Particles{ChN}(i).Approved==1)||(Particles{ChN}(i).Approved==2)

                        %Determine the nc so we can add to the right position of ParticleCountAP
                        if (FrameInfo(min(Particles{ChN}(i).Frame(Particles{ChN}(i).FrameApproved))).nc)>=12
                            CurrentNC=FrameInfo(min(Particles{ChN}(i).Frame(Particles{ChN}(i).FrameApproved))).nc;

                            %Now determine the AP bin this particle is on. We'll use
                            %the nucleus positioning. Also, in order to count it we
                            %need to make sure that it fulfills the same criteria we
                            %used to count NEllipsesAP
                            EllipseAPPosTemp=[];
                            EllipsePosXTemp=[];
                            EllipsePosYTemp=[];
                            RadiusTemp=[];
                            xPosTemp=[];
                            yPosTemp=[];
                            for j=1:length(schnitzcells(Particles{ChN}(i).Nucleus).frames)
                                CurrentEllipse=Ellipses{schnitzcells(Particles{ChN}(i).Nucleus).frames(j)}(schnitzcells(Particles{ChN}(i).Nucleus).cellno(j),:);
                                RadiusTemp=[RadiusTemp,max(CurrentEllipse(3:4))];

                                %Determine the AP position
                                EllipseAPPosTemp=[EllipseAPPosTemp,...
                                    EllipsePos{schnitzcells(Particles{ChN}(i).Nucleus).frames(j)}(schnitzcells(Particles{ChN}(i).Nucleus).cellno(j))];

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

                                MeanAP=mean(EllipseAPPosTemp);

                                ParticleCountAP{ChN}(CurrentNC-11,max(find(APbinID<=MeanAP)))=...
                                    ParticleCountAP{ChN}(CurrentNC-11,max(find(APbinID<=MeanAP)))+1;
                            end
                        end
                    end
                end
            end

            %Calculate the probability using the mean number of nuclei per AP bin
            %in each nc. Note that we look at a reduced range within the nc to
            %reduce variability in counting at mitosis.
            ParticleCountProbAP{ChN}(:,1)=ParticleCountAP{ChN}(1,:)./mean(NEllipsesAP(nc12+5:nc13-5,:));
            ParticleCountProbAP{ChN}(:,2)=ParticleCountAP{ChN}(2,:)./mean(NEllipsesAP(nc13+5:nc14-5,:));
            ParticleCountProbAP{ChN}(:,3)=ParticleCountAP{ChN}(3,:)./...
                mean(NEllipsesAP(max(1,nc14-5):numFrames-5,:));
            % ES 2014-01-08: accounting for movies started fewer than 5 frames before
            % mitosis 13

            figure(16)
            plot(APbinID,ParticleCountProbAP{ChN}(:,1),'.-b')
            hold on
            plot(APbinID,ParticleCountProbAP{ChN}(:,2),'.-k')
            plot(APbinID,ParticleCountProbAP{ChN}(:,3),'.-r')
            hold off
            %ylim([0,max(ParticleCountAP(1,:))*2*1.1])
            xlabel('AP position (x/L)')
            ylabel('Active nuclei')
            title('Number active nuclei')
            legend('nc12', 'nc13', 'nc14')
            saveas(gca,[DropboxFolder,filesep,Prefix,filesep,'Probabilities',filesep,'ProbVsAP_ch',...
                iIndex(ChN,2),'.tif'])

        end


        %Use the alternative approach I used for the movies. We are going to
        %look at each nucleus towards the end of each nc and ask if they
        %correspond to an on or off particle in any frame previous to that one.


        %We'll go for 2.5 minutes before the next mitosis. I might relate this
        %later to the elongation time as a way to say that these particles
        %won't contribute to the total amount of mRNA produced anyway.
        FramesBack=ceil(2.5/mean(diff(ElapsedTime)));

        TotalEllipsesAP=zeros(length(APbinID),3);
        EllipsesOnAP{ChN}=zeros(length(APbinID),3);
        rateOnAP{ChN}=zeros(length(APbinID),3);
        rateOnAPCell{ChN}=cell(length(APbinID),3);     
        timeOnOnAP{ChN}=zeros(length(APbinID),3);
        timeOnOnAPCell{ChN}=cell(length(APbinID),3);

        for nc=12:14

            %Figure out which frame we'll look at
            if nc==14
                FrameToUse=numFrames-FramesBack;
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
                    %Find which AP bin we're in
                    CurrentAPbin=max(find(APbinID<EllipsePos{FrameToUse}(EllipsesToCheck(j))));
                    %Count the total amount of ellipses in the right AP bin
                    TotalEllipsesAP(CurrentAPbin,nc-11)=TotalEllipsesAP(CurrentAPbin,nc-11)+1;



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
                                        MaxDistance=.424/pixelSize; %in microns. Maximum distance to identify an
                                        %ellipse with a schnitz. .424
                                        %corresponds to 2 pixels when using
                                        %212nm resolution.
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
                                            MinValue;
                                            error('What to do here?')
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
                                    EllipsesOnAP{ChN}(CurrentAPbin,nc-11)=EllipsesOnAP{ChN}(CurrentAPbin,nc-11)+1;
                                    if doSingleFits
                                        rateOnAP{ChN}(CurrentAPbin,nc-11) = nansum([rateOnAP{ChN}(CurrentAPbin,nc-11),CompiledParticles{ChN}(m).singleTraceLoadingRate]);
                                        rateOnAPCell{ChN}{CurrentAPbin,nc-11} = [rateOnAPCell{ChN}{CurrentAPbin,nc-11},CompiledParticles{ChN}(m).singleTraceLoadingRate];
                                        timeOnOnAP{ChN}(CurrentAPbin,nc-11) = nansum([timeOnOnAP{ChN}(CurrentAPbin,nc-11),CompiledParticles{ChN}(m).singleTraceTimeOn]);
                                        timeOnOnAPCell{ChN}{CurrentAPbin,nc-11} = [timeOnOnAPCell{ChN}{CurrentAPbin,nc-11},CompiledParticles{ChN}(m).singleTraceTimeOn];
                                    end
                                end
                            end

                        end
                    end
                end
            end
        end

        if doSingleFits
            rateOnAP{ChN} = rateOnAP{ChN} ./ EllipsesOnAP{ChN};
            timeOnOnAP{ChN} = timeOnOnAP{ChN} ./ EllipsesOnAP{ChN};
        end

        if ~SkipAll
            fractionFig = figure();
            fractionAxes = axes(fractionFig);
            plot(fractionAxes,APbinID,EllipsesOnAP{ChN}(:,1)./TotalEllipsesAP(:,1),'.-b') % fraction on nc 12
            hold(fractionAxes,'on')
            plot(fractionAxes,APbinID,EllipsesOnAP{ChN}(:,2)./TotalEllipsesAP(:,2),'.-k') % fraction on nc 13
            plot(fractionAxes,APbinID,EllipsesOnAP{ChN}(:,3)./TotalEllipsesAP(:,3),'.-r') % fraction on nc 14
            hold(fractionAxes,'off')
            title(fractionAxes,'Fraction active nuclei')
            xlabel(fractionAxes,'AP (x/L)')
            ylabel(fractionAxes,'Fraction')
            legend(fractionAxes,'nc12', 'nc13', 'nc14')
        end
    end
end

