function [NEllipsesAP, MeanVectorAllAP, SEVectorAllAP, EllipsesFilteredPos, ...
    FilteredParticlesPos, OnRatioAP, ParticleCountAP, ParticleCountProbAP, ...
    EllipsesOnAP, rateOnAP, rateOnAPCell, timeOnOnAP, timeOnOnAPCell, TotalEllipsesAP, rateOnAPManual, rateOnAPCellManual, timeOnOnAPManual, timeOnOnAPCellManual]...
    ...
    = computeAPFractionOn(NChannels, Particles, schnitzcells, ...
    ...
    CompiledParticles, Ellipses, APbinID, FrameInfo, ElapsedTime, DropboxFolder, ...
    Prefix, EllipsePos, nc12, nc13, nc14, numFrames, SkipFits, SkipAll,...
    APbinArea, pixelSize, manualSingleFits, edgeWidth)
%
% computeAPFractionOn(varargin)
%
% DESCRIPTION
% This function calculates the fraction of transcribing nuclei using three
% different methods.
%
% ARGUMENTS
% NChannels : number of channels
%
% Particles : a structure of the properties of the particles
%
% schnitzcells :
%
% CompiledParticles : a structure for storing properties about the
%                     compiled particles
%
% Ellipses :
%
% APbinID : an array containing the id of each ap bin
%
% FrameInfo : temporal information about the movie
%
% ElapsedTime : an array of elapsed time since the start of the movie
%
% DropboxFolder : dropbox folder for the movie
%
% Prefix : prefix of the movie
%
% EllipsePos : positions of the ellipses
%
% nc12 : frame number of the start of nuclear cycle 12
%
% nc13 : frame number of the start of nuclear cycle 13
%
% nc14 : frame number of the start of nuclear cycle 14
%
% numFrames : number of frames in the movie
%
% SkipFits : a flag to not compile single trace fits generated by
%            fitSingleTraces
%
% SkipAll : a flag to skip all that can be skipped
%
% APbinArea : the area of an ap bin
%
% pixelSize : the pixel size of the data set
%
% manualSingleFits : a flag to compile values from manually generated
%                    single trace fits
%
% OPTIONS
% None

% Author (contact): Hernan Garcia (hggarcia@berkeley.edu)
% Created:
% Last Updated: 3/8/19 (AR)
%
% Documented by: Emma Luu (emma_luu@berkeley.edu)

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

EdgeWidth= edgeWidth/pixelSize; %in microns. 2.12 is simply the number that results in 10 pixel width
%when using a spatial resolution of 212nm.


continueanyway = 0; %3/29/19 JL: Workaround to skip schnitzcell rescuing error.
%Skips the calculation of AP fraction ON if the schnitzcell can't be
%rescued.

if ~SkipAll
    fractionFig = figure();
    fractionAxes = axes(fractionFig);
end

EllipsesOnAP = cell(1, NChannels);

for ChN=1:NChannels
    
    EllipsesOnAP{ChN}=zeros(length(APbinID),3);
    TotalEllipsesAP=zeros(length(APbinID),3);
    
    [NEllipsesAP, OnRatioAP, MeanVectorAllAP, SEVectorAllAP,...
        EllipsesFilteredPos, FilteredParticlesPos, ParticleCountAP, ParticleCountProbAP] =...
        ...
        computeFractionFirstHalf(Particles, ChN, CompiledParticles, Ellipses,...
        schnitzcells, APbinID, pixelSize, FrameInfo, EllipsePos,...
        DropboxFolder, Prefix, ElapsedTime, edgeWidth, SkipAll);
    
    %Use the alternative approach I used for the movies. We are going to
    %look at each nucleus towards the end of each nc and ask if they
    %correspond to an on or off particle in any frame previous to that one.
    
    
    %We'll go for 2.5 minutes before the next mitosis. I might relate this
    %later to the elongation time as a way to say that these particles
    %won't contribute to the total amount of mRNA produced anyway.
    FramesBack=ceil(2.5/mean(diff(ElapsedTime)));
    if FramesBack < 0
        warning('Y2K bug strikes again? ElapsedTime field of FrameInfo is broken.');
        warning('Setting FramesBack to arbitrary number for ellipses and fraction calculations.');
        FramesBack = 15;
    end
    
    %single fits generated by fitSingleTraces.m
    rateOnAP{ChN}=zeros(length(APbinID),3);
    rateOnAPCell{ChN}=cell(length(APbinID),3);
    timeOnOnAP{ChN}=zeros(length(APbinID),3);
    timeOnOnAPCell{ChN}=cell(length(APbinID),3);
    
    %yjk single fits (manually fitted with CheckParticleTracking)
    rateOnAPManual{ChN}=zeros(length(APbinID),3);
    rateOnAPCellManual{ChN}=cell(length(APbinID),3);
    timeOnOnAPManual{ChN}=zeros(length(APbinID),3);
    timeOnOnAPCellManual{ChN}=cell(length(APbinID),3);
    
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
            
            if EdgeWidth == 0
                Radius = 0;
            end
            
            EllipseFilter=(CurrentEllipses(:,1)-Radius-EdgeWidth>0)&...
                (CurrentEllipses(:,2)-Radius-EdgeWidth>0)&...
                (CurrentEllipses(:,1)+Radius+EdgeWidth<FrameInfo(1).PixelsPerLine)&...
                (CurrentEllipses(:,2)+Radius+EdgeWidth<FrameInfo(1).LinesPerFrame);
            
            %Check if the filtered ellipses had an associated particle
            EllipsesToCheck=find(EllipseFilter);
            
            for j=1:length(EllipsesToCheck)
                %Find which AP bin we're in
                CurrentAPbin=find(APbinID<EllipsePos{FrameToUse}(EllipsesToCheck(j)), 1, 'last' );
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
                                if continueanyway == 0
                                    continueanywaycheck = input('Cannot rescue schnitz. Enter y to continue anyway.','s');
                                    if strcmp('y',continueanywaycheck)
                                        continueanyway = 1;
                                    else
                                        error('Cannnot rescue schnitz')
                                    end
                                else
                                end
                            end
                        end
                        
                    end
                    
                    if schnitzcells(k).cellno(IndexToUse)==EllipsesToCheck(j)
                        %Now see if there is an associated particle with it
                        for m=1:length(CompiledParticles{ChN})
                            if CompiledParticles{ChN}(m).Nucleus==k
                                EllipsesOnAP{ChN}(CurrentAPbin,nc-11)=EllipsesOnAP{ChN}(CurrentAPbin,nc-11)+1;
                                if ~SkipFits % if the fitting was done by fitSingleTraces, save the information
                                    rateOnAP{ChN}(CurrentAPbin,nc-11) = nansum([rateOnAP{ChN}(CurrentAPbin,nc-11),CompiledParticles{ChN}(m).singleTraceLoadingRate]);
                                    rateOnAPCell{ChN}{CurrentAPbin,nc-11} = [rateOnAPCell{ChN}{CurrentAPbin,nc-11},CompiledParticles{ChN}(m).singleTraceLoadingRate];
                                    timeOnOnAP{ChN}(CurrentAPbin,nc-11) = nansum([timeOnOnAP{ChN}(CurrentAPbin,nc-11),CompiledParticles{ChN}(m).singleTraceTimeOn]);
                                    timeOnOnAPCell{ChN}{CurrentAPbin,nc-11} = [timeOnOnAPCell{ChN}{CurrentAPbin,nc-11},CompiledParticles{ChN}(m).singleTraceTimeOn];
                                    
                                end
                                if manualSingleFits
                                    rateOnAPManual{ChN}(CurrentAPbin,nc-11) = nansum([rateOnAPManual{ChN}(CurrentAPbin,nc-11),CompiledParticles{ChN}(m).fittedSlope]);
                                    rateOnAPCellManual{ChN}{CurrentAPbin,nc-11} = [rateOnAPCellManual{ChN}{CurrentAPbin,nc-11},CompiledParticles{ChN}(m).fittedSlope];
                                    timeOnOnAPManual{ChN}(CurrentAPbin,nc-11) = nansum([timeOnOnAPManual{ChN}(CurrentAPbin,nc-11),CompiledParticles{ChN}(m).fittedTon]);
                                    timeOnOnAPCellManual{ChN}{CurrentAPbin,nc-11} = [timeOnOnAPCellManual{ChN}{CurrentAPbin,nc-11},CompiledParticles{ChN}(m).fittedTon];
                                end
                            end
                        end
                        
                    end
                end
            end
        end
    end
    
    if ~SkipFits
        rateOnAP{ChN} = rateOnAP{ChN} ./ EllipsesOnAP{ChN};
        timeOnOnAP{ChN} = timeOnOnAP{ChN} ./ EllipsesOnAP{ChN};
    end
    if manualSingleFits
        rateOnAPManual{ChN} = rateOnAPManual{ChN} ./ EllipsesOnAP{ChN};
        timeOnOnAPManual{ChN} = timeOnOnAPManual{ChN} ./ EllipsesOnAP{ChN};
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
