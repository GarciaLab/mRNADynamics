function [hisOverlayHandle, ellipseHandles] =...
    displayOverlays(overlayAxes, cptState, SpeedMode, ShowThreshold2,...
    Overlay, numFrames, UseSchnitz, ZoomRange, fish, multiAx,...
    HisOverlayFigAxes, hisOverlayHandle, ellipseHandles, ImageHisMat)

% JP: refactor this references to call new CPTState methods
Particles = cptState.Particles;
CurrentFrame = cptState.CurrentFrame;
CurrentChannel = cptState.CurrentChannel;
CurrentParticle = cptState.CurrentParticle;
schnitzcells = cptState.schnitzcells;
Ellipses = cptState.Ellipses;

EllipseHandle=[];
EllipseHandleYellow=[];
EllipseHandleBlue=[];
EllipseHandleWhite=[];
EllipseHandleGreen=[];

if ~isempty(multiAx)
    multiView = true;
else
    multiView = false;
end

%Get the coordinates of all the spots in this frame
[x,y,z]=SpotsXYZ(cptState.Spots{CurrentChannel}(CurrentFrame));
%Pull out the right particle if it exists in this frame
CurrentParticleIndex=...
    Particles{CurrentChannel}(CurrentParticle).Index(Particles{CurrentChannel}(CurrentParticle).Frame==...
    CurrentFrame);
xTrace=x(CurrentParticleIndex);
yTrace=y(CurrentParticleIndex);

numParticles = length(Particles{CurrentChannel});

ApprovedParticles=[Particles{CurrentChannel}.Approved];

%These are the positions of all the approved, disapproved and
%unflagged particles

%Approved particles
IndexApprovedParticles=[];
for i=1:numParticles
    if sum(Particles{CurrentChannel}(i).Frame==CurrentFrame)&&...
            sum(Particles{CurrentChannel}(i).Approved==1)
        IndexApprovedParticles=[IndexApprovedParticles,...
            Particles{CurrentChannel}(i).Index(Particles{CurrentChannel}(i).Frame==CurrentFrame)];
    end
end
xApproved=x(IndexApprovedParticles);
yApproved=y(IndexApprovedParticles);

%Disapproved particles
IndexDisapprovedParticles=[];
for i=1:numParticles
    if sum(Particles{CurrentChannel}(i).Frame==CurrentFrame)&&sum(Particles{CurrentChannel}(i).Approved==-1)
        IndexDisapprovedParticles=[IndexDisapprovedParticles,...
            Particles{CurrentChannel}(i).Index(Particles{CurrentChannel}(i).Frame==CurrentFrame)];
    end
end
xDisapproved=x(IndexDisapprovedParticles);
yDisapproved=y(IndexDisapprovedParticles);

%Non-flagged particles (these are particles that have not been
%processed)
IndexNonFlaggedParticles=[];
for i=1:numParticles
    if sum(Particles{CurrentChannel}(i).Frame==CurrentFrame)&&...
            ~(sum(Particles{CurrentChannel}(i).Approved==-1)||sum(Particles{CurrentChannel}(i).Approved==1))
        IndexNonFlaggedParticles=[IndexNonFlaggedParticles,...
            Particles{CurrentChannel}(i).Index(Particles{CurrentChannel}(i).Frame==CurrentFrame)];
    end
end
xNonFlagged=x(IndexNonFlaggedParticles);
yNonFlagged=y(IndexNonFlaggedParticles);

%Show all particles in regular mode
if ~SpeedMode
    plot(overlayAxes,xNonFlagged,yNonFlagged,'ow')
    plot(overlayAxes,xApproved,yApproved,'ob')
    plot(overlayAxes,xDisapproved,yDisapproved,'^r')
    %plot(overlayAxes,x, y, 'sw')
end
%Always show current particle. this indicates the x-y center of the spot
%within the brightest z-slice and may differ from the position
%shown in the snippet image, which is centered at the position with
%the current z-slice.
plot(overlayAxes,xTrace,yTrace,'og')
hold(overlayAxes,'off')

if isfield(cptState.FrameInfo, 'nc')
    set(Overlay,'Name',['Particle: ',num2str(CurrentParticle),'/',num2str(numParticles),...
        ', Frame: ',num2str(CurrentFrame),'/',num2str(numFrames),...
        ', Z: ',num2str(cptState.CurrentZ),'/',num2str(cptState.ZSlices),' nc: ', num2str(cptState.FrameInfo(CurrentFrame).nc),...
        ', Ch: ',num2str(CurrentChannel)])
end
if UseSchnitz
    %Show all the nuclei in regular mode
    if ~SpeedMode
        hold(overlayAxes,'on')
        
        EllipseHandle=notEllipse(Ellipses{CurrentFrame}(:,3),...
            Ellipses{CurrentFrame}(:,4),...
            Ellipses{CurrentFrame}(:,5),...
            Ellipses{CurrentFrame}(:,1)+1,...
            Ellipses{CurrentFrame}(:,2)+1,'r',10, overlayAxes);
        hold(overlayAxes,'off')


        %Show the ones that have been approved

        hold(overlayAxes,'on')
        schnitzCellNo=[];
        for i=1:numParticles
            if Particles{CurrentChannel}(i).Approved==1
                try
                    schnitzIndex=find((schnitzcells(Particles{CurrentChannel}(i).Nucleus).frames)==CurrentFrame);
                    schnitzCellNo=[schnitzCellNo,schnitzcells(Particles{CurrentChannel}(i).Nucleus).cellno(schnitzIndex)];
                catch
                    %can't identify the nucleus for this particle.
                end
            end
        end

        EllipseHandleBlue=notEllipse(Ellipses{CurrentFrame}(schnitzCellNo,3),...
            Ellipses{CurrentFrame}(schnitzCellNo,4),...
            Ellipses{CurrentFrame}(schnitzCellNo,5),...
            Ellipses{CurrentFrame}(schnitzCellNo,1)+1,...
            Ellipses{CurrentFrame}(schnitzCellNo,2)+1,'b',10, overlayAxes);
        hold(overlayAxes,'off')
    end

    %Show the corresponding nucleus
    if ~isempty(Particles{CurrentChannel}(CurrentParticle).Nucleus) && Particles{CurrentChannel}(CurrentParticle).Nucleus > 0
        SchnitzIndex=find(schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).frames==CurrentFrame);
        NucleusIndex=schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).cellno(SchnitzIndex);

        if ~isempty(NucleusIndex)
            hold(overlayAxes,'on')
            EllipseHandleGreen=ellipse(Ellipses{CurrentFrame}(NucleusIndex,3),...
                Ellipses{CurrentFrame}(NucleusIndex,4),...
                Ellipses{CurrentFrame}(NucleusIndex,5),...
                Ellipses{CurrentFrame}(NucleusIndex,1)+1,...
                Ellipses{CurrentFrame}(NucleusIndex,2)+1,[],10, overlayAxes);
            set(EllipseHandleGreen,'Color','g')
            hold(overlayAxes,'off')
        else
            %('Error: Particle without an associated nucleus?')
        end



        %Show the daughter nuclei if applicable
        if isfield(schnitzcells, 'E')
            DaughterE=schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).E;
            DaughterD=schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).D;
            Mother=schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).P;
        else
            DaughterE = 0;
            DaughterD = 0;
            Mother = 0;
        end

        if DaughterE~=0
            SchnitzIndex=find(schnitzcells(DaughterE).frames==CurrentFrame);
            NucleusIndex=schnitzcells(DaughterE).cellno(SchnitzIndex);

            if ~isempty(NucleusIndex)
                hold(overlayAxes,'on')
                EllipseHandleWhite=[EllipseHandleWhite,ellipse(Ellipses{CurrentFrame}(NucleusIndex,3),...
                    Ellipses{CurrentFrame}(NucleusIndex,4),...
                    Ellipses{CurrentFrame}(NucleusIndex,5),...
                    Ellipses{CurrentFrame}(NucleusIndex,1)+1,...
                    Ellipses{CurrentFrame}(NucleusIndex,2)+1, [],10,overlayAxes)];

                hold(overlayAxes,'off')
            else
                %('Error: Particle without an associated nucleus?')
            end
        end

        if DaughterD~=0
            SchnitzIndex=find(schnitzcells(DaughterD).frames==CurrentFrame);
            NucleusIndex=schnitzcells(DaughterD).cellno(SchnitzIndex);

            if ~isempty(NucleusIndex)
                hold(overlayAxes,'on')
                EllipseHandleWhite=[EllipseHandleWhite,ellipse(Ellipses{CurrentFrame}(NucleusIndex,3),...
                    Ellipses{CurrentFrame}(NucleusIndex,4),...
                    Ellipses{CurrentFrame}(NucleusIndex,5),...
                    Ellipses{CurrentFrame}(NucleusIndex,1)+1,...
                    Ellipses{CurrentFrame}(NucleusIndex,2)+1,[],10,overlayAxes)];
                hold(overlayAxes,'off')
            else
                %('Error: Particle without an associated nucleus?')
            end
        end

        if ~isempty(EllipseHandleWhite)
            set(EllipseHandleWhite,'Color','w')
        end

        %Show the mother nucleus if applicable
    

        if Mother~=0
            SchnitzIndex=find(schnitzcells(Mother).frames==CurrentFrame);
            NucleusIndex=schnitzcells(Mother).cellno(SchnitzIndex);

            if ~isempty(NucleusIndex)
                hold(overlayAxes,'on')
                EllipseHandleYellow=ellipse(Ellipses{CurrentFrame}(NucleusIndex,3),...
                    Ellipses{CurrentFrame}(NucleusIndex,4),...
                    Ellipses{CurrentFrame}(NucleusIndex,5),...
                    Ellipses{CurrentFrame}(NucleusIndex,1)+1,...
                    Ellipses{CurrentFrame}(NucleusIndex,2)+1,[],10,overlayAxes);
                set(EllipseHandleYellow,'Color','y')
                hold(overlayAxes,'off')
            else
                %('Error: Particle without an associated nucleus?')
            end
        end

    else
        if cptState.UseHistoneOverlay
            warning('This particle does not have an associated nucleus.');
        end
    end
end

if ApprovedParticles(CurrentParticle)==1
    set(Overlay,'Color','g')
elseif ApprovedParticles(CurrentParticle)==-1
    set(Overlay,'Color','r')
elseif ApprovedParticles(CurrentParticle)==2
    set(Overlay,'Color','y')
else
    set(Overlay,'Color','default')
end

%Show the particles that were under threshold 2.
if ShowThreshold2
    %Get the positions of all the spots in this frame
    [x2,y2]=SpotsXYZ(cptState.Spots{CurrentChannel}(CurrentFrame));
    %Filter those that were under threshold 2.
    CurrentSpotFilter=...
        ~logical(cptState.SpotFilter{CurrentChannel}(CurrentFrame,~isnan(cptState.SpotFilter{CurrentChannel}(CurrentFrame,:))));
    x2=x2(CurrentSpotFilter);
    y2=y2(CurrentSpotFilter);

    hold(overlayAxes,'on')
    plot(overlayAxes,x2,y2,'sr')
    hold(overlayAxes,'off')
end

if cptState.ZoomMode
    %Find the closest frame
    [~,MinIndex]=min((Particles{CurrentChannel}(CurrentParticle).Frame-CurrentFrame).^2);
    if length(MinIndex)>1
        MinIndex=MinIndex(1);
    end
    [cptState.xForZoom,cptState.yForZoom]=...
        SpotsXYZ(cptState.Spots{CurrentChannel}(Particles{CurrentChannel}(CurrentParticle).Frame(MinIndex)));

    cptState.xForZoom=cptState.xForZoom(Particles{CurrentChannel}(CurrentParticle).Index(MinIndex));
    cptState.yForZoom=cptState.yForZoom(Particles{CurrentChannel}(CurrentParticle).Index(MinIndex));

    try
        xlim(overlayAxes,[cptState.xForZoom-ZoomRange,cptState.xForZoom+ZoomRange])
        ylim(overlayAxes,[cptState.yForZoom-ZoomRange/2,cptState.yForZoom+ZoomRange/2])
        if multiView
           for z = 1:length(multiAx)
               for f = 1:length(multiAx)
                    xlim(multiAx{z, f},[cptState.xForZoom-ZoomRange,cptState.xForZoom+ZoomRange])
                    ylim(multiAx{z, f},[cptState.yForZoom-ZoomRange/2,cptState.yForZoom+ZoomRange/2]) 
               end
           end
        end
    catch
        %something's outside the limits of the image
    end
end

if cptState.GlobalZoomMode
    xlim(overlayAxes,[cptState.xForZoom-ZoomRange,cptState.xForZoom+ZoomRange])
    ylim(overlayAxes,[cptState.yForZoom-ZoomRange/2,cptState.yForZoom+ZoomRange/2])
    if multiView
           for z = 1:length(multiAx)
               for f = 1:length(multiAx)
                   xlim(multiAx{z, f},[cptState.xForZoom-ZoomRange,cptState.xForZoom+ZoomRange])
                    ylim(multiAx{z, f},[cptState.yForZoom-ZoomRange/2,cptState.yForZoom+ZoomRange/2]) 
               end
           end
    end
end

if cptState.UseHistoneOverlay

    if isempty(cptState.DisplayRange)
        HisOverlayImageMat=cat(3,mat2gray(ImageHisMat),mat2gray(cptState.ImageMat),zeros(size(cptState.ImageMat)));
    else
        HisOverlayImageMat=cat(3,mat2gray(ImageHisMat,double(cptState.DisplayRange)),mat2gray(cptState.ImageMat),zeros(size(cptState.ImageMat)));
    end
    
    hisOverlayHandle = imshow(HisOverlayImageMat,[],'Border','Tight','Parent',HisOverlayFigAxes);
    
    hold(HisOverlayFigAxes,'on')
    if ~SpeedMode
        plot(HisOverlayFigAxes,xNonFlagged,yNonFlagged,'ow')
        plot(HisOverlayFigAxes,xApproved,yApproved,'ob')
    end
    plot(HisOverlayFigAxes,xTrace,yTrace,'og')
    hold(HisOverlayFigAxes,'off')

    if ShowThreshold2
        hold(HisOverlayFigAxes,'on')
        plot(HisOverlayFigAxes,x2,y2,'sw')
        hold(HisOverlayFigAxes,'off')
    end


    if UseSchnitz

        copyobj(EllipseHandle,HisOverlayFigAxes)
        copyobj(EllipseHandleBlue,HisOverlayFigAxes)
        copyobj(EllipseHandleGreen,HisOverlayFigAxes)
        copyobj(EllipseHandleWhite,HisOverlayFigAxes)
        copyobj(EllipseHandleYellow,HisOverlayFigAxes)

    end

    if cptState.ZoomMode || cptState.GlobalZoomMode
        try
            xlim(HisOverlayFigAxes,[cptState.xForZoom-ZoomRange,cptState.xForZoom+ZoomRange])
            ylim(HisOverlayFigAxes,[cptState.yForZoom-ZoomRange/2,cptState.yForZoom+ZoomRange/2])
        catch
            %something's outside the limits of the image
        end
    end

end

ellipseHandles = {EllipseHandle,EllipseHandleYellow,EllipseHandleBlue,EllipseHandleWhite,EllipseHandleGreen};

end

