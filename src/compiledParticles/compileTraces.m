function [Particles, CompiledParticles, ncFilter, ncFilterID] =...
    compileTraces(NChannels, Particles, HistoneChannel, ...
    schnitzcells, minTime, ExperimentAxis, APbinID, APbinArea, CompiledParticles, ...
    Spots, SkipTraces, nc9, nc10, nc11, nc12, nc13, nc14, ncFilterID, ncFilter, ...
    ElapsedTime, intArea, Ellipses, EllipsePos, AssignedNuclei, PreProcPath, ...
    FilePrefix, Prefix, DropboxFolder, NDigits)
%COMPILETRACES Summary of this function goes here
%   Detailed explanation goes here
h=waitbar(0,'Compiling traces');
for ChN=1:NChannels
    k=1;
    for i=1:length(Particles{ChN})
        waitbar(i/length(Particles{ChN}),h)
        if (Particles{ChN}(i).Approved==1)
            
            for NCh=1:NChannels
                if ~isfield(Particles{NCh},'FrameApproved')
                    for i=1:length(Particles{NCh})
                        Particles{NCh}(i).FrameApproved=true(size(Particles{NCh}(i).Frame));
                    end
                end
            end
            
            %Which frames were approved manually?
            FrameFilter=Particles{ChN}(i).FrameApproved;
            %What is the first frame that was found, regardless of the column
            %condition?
            FirstFrame=Particles{ChN}(i).Frame(min(find(Particles{ChN}(i).FrameApproved)));
            
            %Check that for the remaining frames we got a good z-profile
            %                 for j=1:length(Particles{ChN}(i).Frame)
            %                     ZProfile=fad(ChN).channels(Particles{ChN}(i).Frame(j)).fits.shadowsDog{Particles{ChN}(i).Index(j)};
            %                     [Dummy,ZMax]=max(ZProfile);
            %                     if (ZMax==1)|(ZMax==length(ZProfile))
            %                         FrameFilter(j)=0;
            %                     end
            %                 end
            
            %Should I only keep traces of a certain length? We also just keep
            %the ones that have a real schnitz associated with them
            AnalyzeThisParticle=1;      %Flag to see if this particle should be analyzed.
            
            if HistoneChannel
                if ~((sum(FrameFilter)>0)&...
                        (~isempty(schnitzcells(Particles{ChN}(i).Nucleus).frames)))
                    AnalyzeThisParticle=0;
                end
            elseif ~(sum(FrameFilter)>0)
                AnalyzeThisParticle=0;
            elseif length(Particles{ChN}(i)) <  minTime
                AnalyzeThisParticle=0;
            end
            
            
            
            %See if this particle is in one of the approved AP bins
            try
                if strcmpi(ExperimentAxis,'AP') || strcmpi(ExperimentAxis,'DV')
                    CurrentAPbin=max(find(APbinID<mean(Particles{ChN}(i).APpos(FrameFilter))));
                    if isnan(APbinArea(CurrentAPbin))
                        AnalyzeThisParticle=0;
                    end
                end
            catch
                error(['You probably need to re-run AddParticlePosition again. If that ',...
                    'doesn''t fix things, talk to HG.'])
            end
            
            
            
            if AnalyzeThisParticle
                
                %Reference to the original Particles index
                CompiledParticles{ChN}(k).OriginalParticle=i;
                
                %Copy the filtered information
                CompiledParticles{ChN}(k).Frame=Particles{ChN}(i).Frame(FrameFilter);
                CompiledParticles{ChN}(k).Index=Particles{ChN}(i).Index(FrameFilter);
                CompiledParticles{ChN}(k).xPos=Particles{ChN}(i).xPos(FrameFilter);
                CompiledParticles{ChN}(k).yPos=Particles{ChN}(i).yPos(FrameFilter);
                %(MT, 2018-02-11) Hacky fix to get lattice to run - FIX LATER
                %CompiledParticles{ChN}(k).DVpos=Particles{ChN}(i).DVpos(FrameFilter);
                CompiledParticles{ChN}(k).FrameApproved = Particles{ChN}(i).FrameApproved;
                
                if strcmpi(ExperimentAxis,'AP')
                    CompiledParticles{ChN}(k).APpos=Particles{ChN}(i).APpos(FrameFilter);
                    
                    %Determine the particles average and median AP position
                    CompiledParticles{ChN}(k).MeanAP=mean(Particles{ChN}(i).APpos(FrameFilter));
                    CompiledParticles{ChN}(k).MedianAP=median(Particles{ChN}(i).APpos(FrameFilter));
                elseif strcmpi(ExperimentAxis,'DV')%&isfield(Particles,'APpos')
                    %AP information:
                    CompiledParticles{ChN}(k).APpos=Particles{ChN}(i).APpos(FrameFilter);
                    CompiledParticles{ChN}(k).MeanAP=mean(Particles{ChN}(i).APpos(FrameFilter));
                    CompiledParticles{ChN}(k).MedianAP=median(Particles{ChN}(i).APpos(FrameFilter));
                    %DV information:
                    CompiledParticles{ChN}(k).DVpos=Particles{ChN}(i).DVpos(FrameFilter);
                    CompiledParticles{ChN}(k).MeanDV=mean(Particles{ChN}(i).DVpos(FrameFilter));
                    CompiledParticles{ChN}(k).MedianDV=median(Particles{ChN}(i).DVpos(FrameFilter));
                    
                end
                
                %If we have the histone channel we will actually replace the AP
                %position by the position of the nucleus where the particle was
                %found. If there is no nucleus (like when a particle survives
                %past the nuclear division) we will still use the actual particle
                %position.
                if HistoneChannel&strcmpi(ExperimentAxis,'AP')
                    %Save the original particle position
                    CompiledParticles{ChN}(k).APposParticle=CompiledParticles{ChN}(k).APpos;
                    
                    FramesToCheck=schnitzcells(Particles{ChN}(i).Nucleus).frames(...
                        ismember(schnitzcells(Particles{ChN}(i).Nucleus).frames,Particles{ChN}(i).Frame(FrameFilter)));
                    EllipsesToCheck=schnitzcells(Particles{ChN}(i).Nucleus).cellno(...
                        ismember(schnitzcells(Particles{ChN}(i).Nucleus).frames,Particles{ChN}(i).Frame(FrameFilter)));
                    
                    for j=1:length(FramesToCheck)
                        IndexToChange=find(CompiledParticles{ChN}(k).Frame==FramesToCheck(j));
                        CompiledParticles{ChN}(k).APPos(IndexToChange)=EllipsePos{FramesToCheck(j)}(EllipsesToCheck(j));
                    end
                end
                
                %First frame it was detected at
                CompiledParticles{ChN}(k).FirstFrame=FirstFrame;
                CompiledParticles{ChN}(k).Approved=Particles{ChN}(i).Approved;
                
                %Copy the fit results if they are there
                if isfield(Particles,'Fit')
                    CompiledParticles{ChN}(k).Fit=Particles{ChN}(i).Fit;
                end
                
                
                %Extract information from Spots about fluorescence and background
                [Frame,AmpIntegral, AmpIntegral3, AmpIntegral5, AmpGaussian,...
                    Off, ErrorIntegral,ErrorGauss,optFit1, FitType, ErrorIntegral3, ErrorIntegral5,backGround3]...
                    = GetParticleTrace(k,CompiledParticles{ChN},Spots{ChN});
                CompiledParticles{ChN}(k).Fluo= AmpIntegral;
                CompiledParticles{ChN}(k).Fluo3= AmpIntegral3;
                CompiledParticles{ChN}(k).Fluo5= AmpIntegral5;
                CompiledParticles{ChN}(k).FluoGauss= AmpGaussian;
                CompiledParticles{ChN}(k).Off=Off;
                CompiledParticles{ChN}(k).FluoError=ErrorIntegral(1); % SEANCHANGED
                CompiledParticles{ChN}(k).optFit1=optFit1;
                CompiledParticles{ChN}(k).FitType=FitType;
                
                
                %Determine the nc where this particle was born
                try
                    CompiledParticles{ChN}(k).nc=FrameInfo(CompiledParticles{ChN}(k).Frame(1)).nc;
                catch
                end
                
                if HistoneChannel
                    CompiledParticles{ChN}(k).Nucleus=Particles{ChN}(i).Nucleus;
                    %We have two fields with the same information:
                    %"Nucleus" and "schnitz". In future versions we'll get rid of
                    %"Nucleus"
                    CompiledParticles{ChN}(k).schnitz=Particles{ChN}(i).Nucleus;
                    
                    %Save lineage information in terms of particles
                    if ~isempty(schnitzcells(Particles{ChN}(i).Nucleus).P)
                        if isempty(find(AssignedNuclei{ChN}==schnitzcells(Particles{ChN}(i).Nucleus).P))
                            CompiledParticles{ChN}(k).PParticle=0;
                        else
                            CompiledParticles{ChN}(k).PParticle=find(AssignedNuclei{ChN}==schnitzcells(Particles{ChN}(i).Nucleus).P);
                        end
                    else
                        CompiledParticles{ChN}(k).PParticle=[];
                    end
                    
                    if ~isempty(schnitzcells(Particles{ChN}(i).Nucleus).D)
                        if isempty(find(AssignedNuclei{ChN}==schnitzcells(Particles{ChN}(i).Nucleus).D))
                            CompiledParticles{ChN}(k).DParticle=0;
                        else
                            CompiledParticles{ChN}(k).DParticle=find(AssignedNuclei{ChN}==schnitzcells(Particles{ChN}(i).Nucleus).D);
                        end
                    else
                        CompiledParticles{ChN}(k).DParticle=[];
                    end
                    
                    if ~isempty(schnitzcells(Particles{ChN}(i).Nucleus).E)
                        if isempty(find(AssignedNuclei{ChN}==schnitzcells(Particles{ChN}(i).Nucleus).E))
                            CompiledParticles{ChN}(k).EParticle=0;
                        else
                            CompiledParticles{ChN}(k).EParticle=find(AssignedNuclei{ChN}==schnitzcells(Particles{ChN}(i).Nucleus).E);
                        end
                    else
                        CompiledParticles{ChN}(k).EParticle=[];
                    end
                    
                    %Save information about the nucleus birth and death
                    CompiledParticles{ChN}(k).NucStart=schnitzcells(Particles{ChN}(i).Nucleus).frames(1);
                    CompiledParticles{ChN}(k).NucEnd=schnitzcells(Particles{ChN}(i).Nucleus).frames(end);
                    
                    
                end
                
                
                
                %Plot and save this trace together with its offset value
                
                if ~SkipTraces
                    if ~isnan(nc9)|~isnan(nc10)|~isnan(nc11)|~isnan(nc12)|~isnan(nc13)|~isnan(nc14)
                        %ncFilterID just tells you the identity of the different
                        %filters stored in the cell ncFilter
                        ncFilterID=[];
                        if nc9~=0
                            ncFilterID=9;
                        end
                        if nc10~=0
                            ncFilterID=[ncFilterID,10];
                        end
                        if nc11~=0
                            ncFilterID=[ncFilterID,11];
                        end
                        if nc12~=0
                            ncFilterID=[ncFilterID,12];
                        end
                        if nc13~=0
                            ncFilterID=[ncFilterID,13];
                        end
                        if nc14~=0
                            ncFilterID=[ncFilterID,14];
                        end
                        %Add the first nc
                        ncFilterID=[min(ncFilterID)-1,ncFilterID];
                        
                        
                        %Create the filter
                        % for ChN=1:NChannels
                        
                        if isempty(CompiledParticles)==1
                            error(['No compiled particles found in channel ',num2str(ChN),'. Did you mean to run the code with ApproveAll?'])
                        end
                        
                        ncFilter=false(length(CompiledParticles{ChN})...
                            ,length(ncFilterID)); %AR 6/16/17: I think multi-channel data might require this to be a cell? Something for the future.
                        
                        for i=1:length(CompiledParticles{ChN})
                            %Sometimes CompiledParticles{1}(i).nc is empty. This is because of some
                            %problem with FrameInfo! In that case we'll pull the information out of
                            %the XLS file.
                            if ~isfield(CompiledParticles{ChN}(i), 'nc')
                                CompiledParticles{ChN}(i).nc = [];
                            end
                            if ~isempty(CompiledParticles{ChN}(i).nc)
                                ncFilter(i,find(CompiledParticles{ChN}(i).nc==ncFilterID))=true;
                            else
                                ncsFound=find(CompiledParticles{ChN}(i).Frame(1)>=[nc9,nc10,nc11,nc12,nc13,nc14]);
                                if ncsFound(end)==1
                                    CompiledParticles{ChN}(i).nc=9;
                                    ncFilter(i,ncFilterID==9)=true;
                                elseif ncsFound(end)==2
                                    CompiledParticles{ChN}(i).nc=10;
                                    ncFilter(i,ncFilterID==10)=true;
                                elseif ncsFound(end)==3
                                    CompiledParticles{ChN}(i).nc=11;
                                    ncFilter(i,ncFilterID==11)=true;
                                elseif ncsFound(end)==4
                                    CompiledParticles{ChN}(i).nc=12;
                                    ncFilter(i,ncFilterID==12)=true;
                                elseif ncsFound(end)==5
                                    CompiledParticles{ChN}(i).nc=13;
                                    ncFilter(i,ncFilterID==13)=true;
                                elseif ncsFound(end)==6
                                    CompiledParticles{ChN}(i).nc=14;
                                    ncFilter(i,ncFilterID==14)=true;
                                end
                            end
                        end
                        %end
                    end
                    
                    figure(2)
                    left_color = [213,108,85]/255;
                    right_color = [0, 0, 0]/255;
                    set(gcf,'defaultAxesColorOrder',[left_color; right_color]);
                    %Size of the snippet for each frame
                    SnippetSize=31; %AR 3/15/16: Why is this 31?
                    %Width of the screen
                    ScreenWidth=get( 0, 'ScreenSize' );
                    ScreenWidth=ScreenWidth(3);
                    
                    %Figure out the arrangement of snippets
                    NFrames=length(CompiledParticles{ChN}(k).Frame);
                    
                    NRows=ceil(NFrames/ScreenWidth*SnippetSize)*2;
                    NCols=max([2,ceil(NFrames/NRows)]);
                    
                    %Actual total number of rows
                    TotalRows=14;
                    
                    subplot(TotalRows,NCols,[1:((TotalRows-NRows)*NCols)])
                    
                    
                    
                    %Top left plot
                    FilterMatrix=zeros((TotalRows-NRows),NCols);
                    FilterMatrix(:,1:ceil(NCols/2))=1;
                    subplot(TotalRows,NCols,find(FilterMatrix'))
                    
                    %                     yyaxis left
                    errorbar(ElapsedTime(CompiledParticles{ChN}(k).Frame),...
                        CompiledParticles{ChN}(k).Fluo,ones(size(CompiledParticles{ChN}(k).Fluo))*...
                        CompiledParticles{ChN}(k).FluoError,...
                        '.-r');
                    ylabel('fluorescence (au)')
                    hold on
                    
                    %                     yyaxis right
                    plot(ElapsedTime(CompiledParticles{ChN}(k).Frame),...
                        CompiledParticles{ChN}(k).Off*intArea,'.-g');
                    if ~isempty(CompiledParticles{ChN}(k).optFit1)
                        
                        if strcmp(CompiledParticles{ChN}(k).FitType,'spline')
                            SplineValues=ppval(CompiledParticles{ChN}(k).optFit1,double(CompiledParticles{ChN}(k).Frame));
                        elseif strcmp(CompiledParticles{ChN}(k).FitType,'mean')
                            SplineValues=ones(size(CompiledParticles{ChN}(k).Frame))*CompiledParticles{ChN}(k).optFit1;
                        elseif strcmp(CompiledParticles{ChN}(k).FitType,'line')
                            SplineValues=polyval(CompiledParticles{ChN}(k).optFit1,CompiledParticles{ChN}(k).Frame);
                        end
                        
                        %                         yyaxis right
                        plot(ElapsedTime(CompiledParticles{ChN}(k).Frame),SplineValues*intArea,'-b')
                        try
                            title(['particle ',num2str(k),'(',num2str(i),'), nc',num2str(CompiledParticles{ChN}(k).nc),', Ch: ',num2str(ChN)])
                        catch
                        end
                    else
                        title(['particle ',num2str(k),'(',num2str(i),'), nc',num2str(CompiledParticles{1}(k).nc),', Ch: ',num2str(ChN),...
                            ' - WARNING: No offset fit'])
                    end
                    hold off
                    legend({'Particle','Offset','Offset fit'},'Location','Best')
                    xlabel('time (min)')
                    axis square
                    set(gca, 'Position', get(gca, 'OuterPosition') - ...
                        get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
                    drawnow
                    
                    
                    
                    
                    %Top right plot
                    if HistoneChannel
                        subplot(TotalRows,NCols,find(~FilterMatrix'))
                        
                        if length(CompiledParticles{ChN}(k).Frame)>1
                            colormap(jet(128));
                            cmap=colormap;
                            
                            ColorTime=[];
                            for j=1:length(CompiledParticles{ChN}(k).Frame)
                                ColorTime(j,:)= cmap(round((j-1)/...
                                    (length(CompiledParticles{ChN}(k).Frame)-1)*127+1),:);
                            end
                        else
                            ColorTime=[];
                            ColorTime(1,:)=[1,0,0];
                        end
                        
                        
                        
                        hold on
                        for j=1:length(CompiledParticles{ChN}(k).Frame)
                            PosSchnitz=find((schnitzcells(CompiledParticles{ChN}(k).Nucleus).frames)==...
                                CompiledParticles{ChN}(k).Frame(j));
                            PosEllipse=schnitzcells(CompiledParticles{ChN}(k).Nucleus).cellno(PosSchnitz);
                            CurrEllipse=Ellipses{CompiledParticles{ChN}(k).Frame(j)}(PosEllipse,:);
                            
                            if ~isempty(CurrEllipse)
                                EllipseHandle=ellipse(CurrEllipse(3),...
                                    CurrEllipse(4),...
                                    CurrEllipse(5),...
                                    0,0,[],[],gca);
                                set(EllipseHandle,'color',ColorTime(j,:))
                                plot(CompiledParticles{ChN}(k).xPos(j)-CurrEllipse(1),...
                                    CompiledParticles{ChN}(k).yPos(j)-CurrEllipse(2),'o','color',...
                                    ColorTime(j,:))
                            else
                                PosSchnitz=length(schnitzcells(CompiledParticles{ChN}(k).Nucleus).frames);
                                PosEllipse=schnitzcells(CompiledParticles{ChN}(k).Nucleus).cellno(PosSchnitz);
                                CurrEllipse=...
                                    Ellipses{schnitzcells(CompiledParticles{ChN}(k).Nucleus).frames(PosSchnitz)}(PosEllipse,:);
                                
                                
                                
                                plot(CompiledParticles{ChN}(k).xPos(j)-CurrEllipse(1),...
                                    CompiledParticles{ChN}(k).yPos(j)-CurrEllipse(2),'o','color',...
                                    ColorTime(j,:))
                            end
                            
                            
                            
                        end
                        hold off
                        box on
                        xlabel('x position (pixels)')
                        ylabel('y position (pixels)')
                        axis square
                        set(gca, 'Position', get(gca, 'OuterPosition') - ...
                            get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
                        set(gca, 'YDir', 'reverse')
                        BarHandle = colorbar;
                        set(BarHandle,'YTick',[])
                        
                        %%%%%AR 6/16/17: This still doesn't work in 2017. We need to find a
                        %%%%%replacement.
                        %                         try
                        %                             if ~isempty(cbfreeze(BarHandle))
                        %                                 BarHandle=cbfreeze(BarHandle);
                        %                             else
                        %                                 warning('Issue with cbfreeze.m. Skipping it. The color bar will not reflect time appropriately.')
                        %                             end
                        %                             ylabel(BarHandle,'Time')
                        %                         catch
                        %                             warning('Issue with cbfreeze.m. Skipping it. The color bar will not reflect time appropriately. This is an issue of Matlab 2014b.')
                        %                         end
                        %%%%%%%%%%%%%%%%
                        if strcmpi(ExperimentAxis,'AP')
                            title(['Mean AP: ',num2str(CompiledParticles{ChN}(k).MeanAP)])
                        end
                        drawnow
                    end 
                    
                    %Snippets
                    for j=1:NFrames
                        subplot(TotalRows,NCols,(TotalRows-NRows)*NCols+j)
                        spotFrame = CompiledParticles{ChN}(k).Frame(j);
                        [x,y,z]=SpotsXYZ(Spots{ChN}(spotFrame));
                        if ~isempty(x)
                            xTrace=x(CompiledParticles{ChN}(k).Index(j));
                            yTrace=y(CompiledParticles{ChN}(k).Index(j));
                            zTrace=z(CompiledParticles{ChN}(k).Index(j));
                            Image=imread([PreProcPath,filesep,FilePrefix(1:end-1),filesep,...
                                FilePrefix,iIndex(CompiledParticles{ChN}(k).Frame(j),NDigits),'_z',iIndex(zTrace,2),...
                                '_ch',iIndex(ChN,2),'.tif']);
                            [ImRows,ImCols]=size(Image);
                            ImageSnippet=zeros(SnippetSize,SnippetSize);
                            yRange=round(yTrace)+[-(SnippetSize-1)/2:(SnippetSize-1)/2];
                            yFilter=(yRange>0)&(yRange<=ImRows);
                            xRange=round(xTrace)+[-(SnippetSize-1)/2:(SnippetSize-1)/2];
                            xFilter=(xRange>0)&(xRange<=ImCols);
                            ImageSnippet(yFilter,xFilter)=Image(yRange(yFilter),...
                                xRange(xFilter));
                            imshow(ImageSnippet,[],'Border','Tight','InitialMagnification',200)
                            set(gca, 'Position', get(gca, 'OuterPosition') - ...
                                get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
                            if HistoneChannel
                                %Plot the corresponding nucleus
                                CurrentSchnitz=schnitzcells(CompiledParticles{ChN}(k).Nucleus);
                                if sum((CurrentSchnitz.frames)==CompiledParticles{ChN}(k).Frame(j))==1
                                    hold on
                                    EllipseNumber=CurrentSchnitz.cellno(...
                                        (CurrentSchnitz.frames)==CompiledParticles{ChN}(k).Frame(j));
                                    CurrEllipse=Ellipses{CompiledParticles{ChN}(k).Frame(j)}(EllipseNumber,:);
                                    EllipseHandle=ellipse(CurrEllipse(3),...
                                        CurrEllipse(4),...
                                        CurrEllipse(5),...
                                        CurrEllipse(1)-xTrace+(SnippetSize-1)/2,...
                                        CurrEllipse(2)-yTrace+(SnippetSize-1)/2,[],[],gca);
                                    %set(EllipseHandle,'color',ColorTime(j,:))
                                    set(EllipseHandle,'color','g')
                                    hold off
                                end
                            end
                        end
                    end
                    set(gcf,'Position',[1,41,1280,684])
                    
                    drawnow
                    if isfield(CompiledParticles{ChN}(k), 'nc')
                        saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'ParticleTraces',filesep,iIndex(k,NDigits),...
                            '(',num2str(i),')-nc',...
                            num2str(CompiledParticles{ChN}(k).nc),'_ch',iIndex(ChN,2),'.tif'])
                    end
                    close(2)
                end
                
                
                k=k+1;
                
            end
        end
    end
end
close(h)
end
