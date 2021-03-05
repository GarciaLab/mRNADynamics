function [Particles, CompiledParticles, ncFilter, ncFilterID] =...
    ...
    compileTraces(NChannels, Particles, HistoneChannel, ...
    ...
    schnitzcells, minTime, ExperimentAxis, APbinID, APbinArea, CompiledParticles, ...
    Spots, SkipTraces, ncFilterID, ncFilter, ...
    ElapsedTime, Ellipses, EllipsePos, PreProcPath, ...
    FilePrefix, Prefix, DropboxFolder, numFrames, manualSingleFits,...
    edgeWidth, pixelSize_nm, coatChannels, fullEmbryoExists, liveExperiment)

%COMPILETRACES Summary of this function goes here
%   Detailed explanation goes here

liveExperiment = LiveExperiment(Prefix);

FrameInfo = getFrameInfo(liveExperiment);

anaphaseFrames = liveExperiment.anaphaseFrames';
nc9 = anaphaseFrames(1); nc10 = anaphaseFrames(2); nc11 = anaphaseFrames(3);
nc12 = anaphaseFrames(4); nc13 = anaphaseFrames(5); nc14 = anaphaseFrames(6);

NDigits = liveExperiment.nDigits;

DropboxFolder = liveExperiment.userResultsFolder;
PreProcPath = liveExperiment.userPreFolder;
pixelSize_nm = liveExperiment.pixelSize_nm;
SnippetSize = 2 * (floor(1300 / (2 * pixelSize_nm))) + 1; % nm. note that this is forced to be odd

NChannels = length(coatChannels);

if ~SkipTraces
    movieMat = getMovieMat(liveExperiment); 
end

h = waitbar(0,'Compiling traces');
for ChN=1:NChannels
    k=1;
    for i=1:length(Particles{ChN})
        try waitbar(i/length(Particles{ChN}),h); end
        if (Particles{ChN}(i).Approved==1)
            
            for NCh=1:NChannels
                if ~isfield(Particles{NCh},'FrameApproved')
                    for m=1:length(Particles{NCh})
                        Particles{NCh}(m).FrameApproved=true(size(Particles{NCh}(m).Frame));
                    end
                end
            end
            
            %Which frames were approved manually?
            FrameFilter=Particles{ChN}(i).FrameApproved;
            %What is the first frame that was found, regardless of the column
            %condition?
            FirstFrame=Particles{ChN}(i).Frame(min(find(Particles{ChN}(i).FrameApproved)));
            
            
            %Should I only keep traces of a certain length? We also just keep
            %the ones that have a real schnitz associated with them
            AnalyzeThisParticle=1;      %Flag to see if this particle should be analyzed.
            
            if HistoneChannel
                %3/29/19 JL Bug workaround: if there are any nuclei with
                %label 0 the following line breaks, due to zero indexing.
                %For now, just skip.
                try
                    if Particles{ChN}(i).Nucleus == 0
                        AnalyzeThisParticle=0;
                    elseif ~((sum(FrameFilter)>0)&...
                            (~isempty(schnitzcells(Particles{ChN}(i).Nucleus).frames)))
                        AnalyzeThisParticle=0;
                    end
                catch
                    error('Particle possibly not associated with nucleus. May lead to inaccuracies in fraction calcluations.')
                    
                end
            elseif ~(sum(FrameFilter)>0)
                AnalyzeThisParticle=0;
            elseif length(Particles{ChN}(i)) <  minTime
                AnalyzeThisParticle=0;
            end
            
            
            
            %See if this particle is in one of the approved AP bins
            try
                if (strcmpi(ExperimentAxis,'AP') || strcmpi(ExperimentAxis,'DV')) && fullEmbryoExists
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
                if isfield(Particles{ChN}(k),'zPos')
                    CompiledParticles{ChN}(k).zPos=Particles{ChN}(i).zPos(FrameFilter);
                end
                %(MT, 2018-02-11) Hacky fix to get lattice to run - FIX LATER
                %CompiledParticles{ChN}(k).DVpos=Particles{ChN}(i).DVpos(FrameFilter);
                CompiledParticles{ChN}(k).FrameApproved = Particles{ChN}(i).FrameApproved;
                
                if strcmpi(ExperimentAxis,'AP') && fullEmbryoExists
                    CompiledParticles{ChN}(k).APpos=Particles{ChN}(i).APpos(FrameFilter);
                    
                    %Determine the particles average and median AP position
                    CompiledParticles{ChN}(k).MeanAP=mean(Particles{ChN}(i).APpos(FrameFilter));
                    CompiledParticles{ChN}(k).MedianAP=median(Particles{ChN}(i).APpos(FrameFilter));
                elseif strcmpi(ExperimentAxis,'DV') && fullEmbryoExists %&isfield(Particles,'APpos')
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
                if HistoneChannel&strcmpi(ExperimentAxis,'AP') && fullEmbryoExists
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
                
                % Save the manually fitted initial slope and T_ON if it
                % exists
                if manualSingleFits && Particles{ChN}(i).fitApproved==1
                    CompiledParticles{ChN}(k).fittedSlope = Particles{ChN}(i).fittedSlope;
                    CompiledParticles{ChN}(k).fittedTon = Particles{ChN}(i).fittedTON;
                else
                    CompiledParticles{ChN}(k).fittedSlope = [];
                    CompiledParticles{ChN}(k).fittedTon = [];
                end
                
%                 %Extract position info and general Gauss3D fit info
%                 try
%                     [~,gx_vec,gy_vec,gz_vec,g_fits_cell,f3_vec,f3_raw_vec]=...
%                         getGauss3DFitInfo(k,CompiledParticles{ChN},Spots{ChN});
%                     threeDFlag = ~all(isnan(gx_vec));
%                     
%                 catch
%                     
%                     warning('Didn''t have complete Gauss 3D info');
%                     threeDFlag = false;
%                     
%                 end
                
%                 threeDFlag = false; %setting this to false until gaussian fits are fixed
                
%                 if threeDFlag
%                     CompiledParticles{ChN}(k).xPosGauss3D = gx_vec;
%                     CompiledParticles{ChN}(k).yPosGauss3D = gy_vec;
%                     CompiledParticles{ChN}(k).zPosGauss3D = gz_vec;
%                     CompiledParticles{ChN}(k).Fluo3DGauss = f3_vec;
%                     CompiledParticles{ChN}(k).Fluo3DRaw = f3_raw_vec;
%                     CompiledParticles{ChN}(k).zPosGauss3D = gz_vec;
%                     CompiledParticles{ChN}(k).fitParamsGauss3D = g_fits_cell;
%                 end

                % add 3D info
                threeDFlag = isfield(Spots{ChN}(Particles{ChN}(i).Frame(1)).Fits(Particles{ChN}(i).Index(1)),'GaussIntensity3DRaw');                
                if threeDFlag
                    CompiledParticles{ChN}(k).Fluo3DGauss = NaN(size(Particles{ChN}(i).Frame));
                    CompiledParticles{ChN}(k).Fluo3DRaw = NaN(size(Particles{ChN}(i).Frame));
                    CompiledParticles{ChN}(k).yPos3D = NaN(size(Particles{ChN}(i).Frame));
                    CompiledParticles{ChN}(k).xPos3D = NaN(size(Particles{ChN}(i).Frame));
                    CompiledParticles{ChN}(k).zPos3D = NaN(size(Particles{ChN}(i).Frame));
                    fitInfoCell = cell(size(Particles{ChN}(i).Frame));
                    for f = 1:length(Particles{ChN}(i).Frame)
                        spot = SpotsCh(Particles{ChN}(i).Frame(f)).Fits(Particles{ChN}(i).Index(f));

                        CompiledParticles{ChN}(k).Fluo3DRaw(f) = spot.GaussIntensity3DRaw;
                        CompiledParticles{ChN}(k).Fluo3DGauss(f) = spot.Gauss3DIntensity;  
                        CompiledParticles{ChN}(k).yPos3D(f) = spot.GaussPos3D(1);
                        CompiledParticles{ChN}(k).xPos3D(f) = spot.GaussPos3D(2);
                        CompiledParticles{ChN}(k).zPos3D(f) = spot.GaussPos3D(3);
                        fitInfoCell{f} = spot.SpotFitInfo3D;
                    end
                                                CompiledParticles{ChN}(k).SpotFitInfo3D = fitInfoCell;
                end 
                %Extract information from Spots about fluorescence and background
                plotTraceSettings = PlotTraceSettings();
                
                [Frame, AmpGaussian, Off, ErrorGauss, optFit1, FitType,...
                    AmpDog, AmpDogMax, ampdog3, ampdog3Max]...
                    = GetParticleTrace(k,CompiledParticles{ChN},Spots{ChN}, plotTraceSettings, false);
                CompiledParticles{ChN}(k).Fluo = plotTraceSettings.AmpIntegral;
                CompiledParticles{ChN}(k).Fluo3 = plotTraceSettings.AmpIntegral3;
                CompiledParticles{ChN}(k).FluoGauss= AmpGaussian;
                CompiledParticles{ChN}(k).Off=Off;
                if ~isempty(plotTraceSettings.ErrorIntegral)
                    CompiledParticles{ChN}(k).FluoError = plotTraceSettings.ErrorIntegral(1); % SEANCHANGED
                else
                    CompiledParticles{ChN}(k).FluoError = NaN;
                end
                CompiledParticles{ChN}(k).optFit1=optFit1;
                CompiledParticles{ChN}(k).FitType=FitType;
                CompiledParticles{ChN}(k).FluoDog = AmpDog;
                CompiledParticles{ChN}(k).FluoDogMax = AmpDogMax;
                ampIntegralGauss3DAux = plotTraceSettings.AmpIntegralGauss3D;
                CompiledParticles{ChN}(k).FluoGauss3D = ampIntegralGauss3DAux';
                CompiledParticles{ChN}(k).FluoGauss3DError = plotTraceSettings.ErrorIntegralGauss3D;
                CompiledParticles{ChN}(k).ampdog3 = ampdog3;
                CompiledParticles{ChN}(k).ampdog3Max = ampdog3Max;
                
                
                
                
                %Determine the nc where this particle was born
                try
                    CompiledParticles{ChN}(k).nc=FrameInfo(CompiledParticles{ChN}(k).Frame(1)).nc;
                    CompiledParticles{ChN}(k).cycle=FrameInfo(CompiledParticles{ChN}(k).Frame(1)).nc;
                end
                
                if HistoneChannel
                    CompiledParticles{ChN}(k).Nucleus=Particles{ChN}(i).Nucleus;
                    
                    %We have two fields with the same information:
                    %"Nucleus" and "schnitz". In future versions we'll get rid of
                    %"Nucleus"
                    CompiledParticles{ChN}(k).schnitz=Particles{ChN}(i).Nucleus;
                    
                    assert(CompiledParticles{ChN}(k).schnitz <= length(schnitzcells));
                    
                    %Save lineage information in terms of particles
                    if isfield(schnitzcells, 'P')
                        
                        if ~isempty(schnitzcells(Particles{ChN}(i).Nucleus).P)
                            if isempty(find([Particles{ChN}.Nucleus]==schnitzcells(Particles{ChN}(i).Nucleus).P))
                                CompiledParticles{ChN}(k).PParticle=0;
                            else
                                CompiledParticles{ChN}(k).PParticle=find([Particles{ChN}.Nucleus]==schnitzcells(Particles{ChN}(i).Nucleus).P);
                            end
                        else
                            CompiledParticles{ChN}(k).PParticle=[];
                        end
                        
                        if ~isempty(schnitzcells(Particles{ChN}(i).Nucleus).D)
                            if isempty(find([Particles{ChN}.Nucleus]==schnitzcells(Particles{ChN}(i).Nucleus).D))
                                CompiledParticles{ChN}(k).DParticle=0;
                            else
                                CompiledParticles{ChN}(k).DParticle=find([Particles{ChN}.Nucleus]==schnitzcells(Particles{ChN}(i).Nucleus).D);
                            end
                        else
                            CompiledParticles{ChN}(k).DParticle=[];
                        end
                        
                        if ~isempty(schnitzcells(Particles{ChN}(i).Nucleus).E)
                            if isempty(find([Particles{ChN}.Nucleus]==schnitzcells(Particles{ChN}(i).Nucleus).E))
                                CompiledParticles{ChN}(k).EParticle=0;
                            else
                                CompiledParticles{ChN}(k).EParticle=find([Particles{ChN}.Nucleus]==schnitzcells(Particles{ChN}(i).Nucleus).E);
                            end
                        else
                            CompiledParticles{ChN}(k).EParticle=[];
                        end
                        
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
                    
                    
                    titlestr = ['ParticleTraces',filesep,iIndex(k,NDigits),...
                        '(',num2str(i),')-nc',...
                        num2str(CompiledParticles{ChN}(k).nc),'_ch',iIndex(ChN,2)];
                    
                    figure('Name',titlestr,'NumberTitle','off');
                    
                    left_color = [213,108,85]/255;
                    right_color = [0, 0, 0]/255;
                    set(gcf,'defaultAxesColorOrder',[left_color; right_color]);
                    %Size of the snippet for each frame
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
%                     errorbar(ElapsedTime(CompiledParticles{ChN}(k).Frame),...
%                         CompiledParticles{ChN}(k).Fluo,ones(size(CompiledParticles{ChN}(k).Fluo))*...
%                         CompiledParticles{ChN}(k).FluoError,...
%                         '.-r');
                    errorbar(ElapsedTime(CompiledParticles{ChN}(k).Frame),...
                        CompiledParticles{ChN}(k).Fluo3, 3*CompiledParticles{ChN}(k).FluoError *ones(size(CompiledParticles{ChN}(k).Fluo)),...
                        '.-r');
                    ylabel('fluorescence (au)')
                    %                     hold on
                    
                    %                     yyaxis right
                    %                     plot(ElapsedTime(CompiledParticles{ChN}(k).Frame),...
                    %                         CompiledParticles{ChN}(k).Off*intArea,'.-g');
                    %                     if ~isempty(CompiledParticles{ChN}(k).optFit1)
                    %
                    %                         if strcmp(CompiledParticles{ChN}(k).FitType,'spline')
                    %                             SplineValues=ppval(CompiledParticles{ChN}(k).optFit1,double(CompiledParticles{ChN}(k).Frame));
                    %                         elseif strcmp(CompiledParticles{ChN}(k).FitType,'mean')
                    %                             SplineValues=ones(size(CompiledParticles{ChN}(k).Frame))*CompiledParticles{ChN}(k).optFit1;
                    %                         elseif strcmp(CompiledParticles{ChN}(k).FitType,'line')
                    %                             SplineValues=polyval(CompiledParticles{ChN}(k).optFit1,CompiledParticles{ChN}(k).Frame);
                    %                         end
                    %
                    %                         %                         yyaxis right
                    %                         plot(ElapsedTime(CompiledParticles{ChN}(k).Frame),SplineValues*intArea,'-b')
                    %                         try
                    %                             title(['particle ',num2str(k),'(',num2str(i),'), nc',num2str(CompiledParticles{ChN}(k).nc),', Ch: ',num2str(ChN)])
                    %                         catch
                    %                         end
                    %                     else
                    %                         title(['particle ',num2str(k),'(',num2str(i),'), nc',num2str(CompiledParticles{1}(k).nc),', Ch: ',num2str(ChN),...
                    %                             ' - WARNING: No offset fit'])
                    %                     end
                    %                     hold off
                    %                     legend({'Particle','Offset','Offset fit'},'Location','Best')
                    xlabel('time (min)')
                    axis square
                    set(gca, 'Position', get(gca, 'OuterPosition') - ...
                        get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
                    standardizeFigure(gca, []);
                    drawnow
                    
                    
                    
                    %                     %%
                    %                     %Top right plot
                    %                     if HistoneChannel
                    %                         subplot(TotalRows,NCols,find(~FilterMatrix'))
                    %
                    %                         if length(CompiledParticles{ChN}(k).Frame)>1
                    %                             colormap(jet(128));
                    %                             cmap=colormap;
                    %
                    %                             ColorTime=[];
                    %                             for j=1:length(CompiledParticles{ChN}(k).Frame)
                    %                                 ColorTime(j,:)= cmap(round((j-1)/...
                    %                                     (length(CompiledParticles{ChN}(k).Frame)-1)*127+1),:);
                    %                             end
                    %                         else
                    %                             ColorTime=[];
                    %                             ColorTime(1,:)=[1,0,0];
                    %                         end
                    %
                    %
                    %
                    %                         hold on
                    %                         for j=1:length(CompiledParticles{ChN}(k).Frame)
                    %                             PosSchnitz=find((schnitzcells(CompiledParticles{ChN}(k).Nucleus).frames)==...
                    %                                 CompiledParticles{ChN}(k).Frame(j));
                    %                             PosEllipse=schnitzcells(CompiledParticles{ChN}(k).Nucleus).cellno(PosSchnitz);
                    %                             CurrEllipse=Ellipses{CompiledParticles{ChN}(k).Frame(j)}(PosEllipse,:);
                    %
                    %                             if ~isempty(CurrEllipse)
                    %                                 EllipseHandle=ellipse(CurrEllipse(3),...
                    %                                     CurrEllipse(4),...
                    %                                     CurrEllipse(5),...
                    %                                     0,0,[],[],gca);
                    %                                 set(EllipseHandle,'color',ColorTime(j,:))
                    %                                 plot(CompiledParticles{ChN}(k).xPos(j)-CurrEllipse(1),...
                    %                                     CompiledParticles{ChN}(k).yPos(j)-CurrEllipse(2),'o','color',...
                    %                                     ColorTime(j,:))
                    %                             else
                    %                                 PosSchnitz=length(schnitzcells(CompiledParticles{ChN}(k).Nucleus).frames);
                    %                                 PosEllipse=schnitzcells(CompiledParticles{ChN}(k).Nucleus).cellno(PosSchnitz);
                    %                                 CurrEllipse=...
                    %                                     Ellipses{schnitzcells(CompiledParticles{ChN}(k).Nucleus).frames(PosSchnitz)}(PosEllipse,:);
                    %
                    %
                    %
                    %                                 plot(CompiledParticles{ChN}(k).xPos(j)-CurrEllipse(1),...
                    %                                     CompiledParticles{ChN}(k).yPos(j)-CurrEllipse(2),'o','color',...
                    %                                     ColorTime(j,:))
                    %                             end
                    %
                    %
                    %
                    %                         end
                    %                         hold off
                    %                         box on
                    %                         xlabel('x position (pixels)')
                    %                         ylabel('y position (pixels)')
                    %                         axis square
                    %                         set(gca, 'Position', get(gca, 'OuterPosition') - ...
                    %                             get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
                    %                         set(gca, 'YDir', 'reverse')
                    %                         BarHandle = colorbar;
                    %                         set(BarHandle,'YTick',[])
                    %
                    %                         if strcmpi(ExperimentAxis,'AP')
                    %                             title(['Mean AP: ',num2str(CompiledParticles{ChN}(k).MeanAP)])
                    %                         end
                    %                         drawnow
                    %                     end
                    %%
                    %Snippets
                    for j=1:NFrames
                        subplot(TotalRows,NCols,(TotalRows-NRows)*NCols+j)
                        spotFrame = CompiledParticles{ChN}(k).Frame(j);
                        [x,y,z]=SpotsXYZ(Spots{ChN}(spotFrame));
                        %                             [x,y,z]=SpotsXYZ(Spots{ChN}(spotFrame), true);
                        if ~isempty(x)
                            if threeDFlag
                                xTrace = round(CompiledParticles{ChN}(k).xPosGauss3D(j));
                                yTrace = round(CompiledParticles{ChN}(k).yPosGauss3D(j));
                                zTrace = round(CompiledParticles{ChN}(k).zPosGauss3D(j));
                            else
                                yTrace=x(CompiledParticles{ChN}(k).Index(j));
                                xTrace=y(CompiledParticles{ChN}(k).Index(j));
                                zTrace=z(CompiledParticles{ChN}(k).Index(j));
                            end
%                             try
                                if isempty(movieMat)
                                    Image = getMovieSlice(liveExperiment, CompiledParticles{ChN}(k).Frame(j), coatChannels(ChN), zTrace );
                                else
                                    Image = movieMat(:, :, zTrace, CompiledParticles{ChN}(k).Frame(j), coatChannels(ChN));
                                end
%                           
                                
                                [ImRows,ImCols]=size(Image);
                                ImageSnippet=zeros(SnippetSize,SnippetSize);
                                yRange= round(yTrace)+ [-(SnippetSize-1)/2:(SnippetSize-1)/2];
                                yFilter=(yRange>0)&(yRange<=ImRows);
                                xRange=round(xTrace)+[-(SnippetSize-1)/2:(SnippetSize-1)/2];
                                xFilter=(xRange>0)&(xRange<=ImCols);
                                
                                ImageSnippet(yFilter,xFilter)=Image(yRange(yFilter),...
                                    xRange(xFilter));
                                
                                imshow(ImageSnippet,[],'Border','Tight','InitialMagnification',200)
                                
                                set(gca, 'Position', get(gca, 'OuterPosition') - ...
                                    get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
                                
                                showNucleus = false;
                                if HistoneChannel && showNucleus
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
%                             end %try
                        end
                    end
                    set(gcf,'Position',[1,41,1280,684])
                    
                    drawnow
               
                    if isfield(CompiledParticles{ChN}(k), 'nc')
                        saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,titlestr,'.tif'])
                    end
                    close(gcf)
                end
                
                
                k=k+1;
                
            end
        end
    end
end

try close(h); end

end

