function [ImageHis, xForZoom, yForZoom] =...
    ...
    displayOverlays(...
    ...
    overlayAxes, Image, SpeedMode, FrameInfo, Particles, ...
    Spots, CurrentFrame, ShowThreshold2, ...
    Overlay, CurrentChannel, CurrentParticle, ZSlices, CurrentZ, numFrames, ...
    schnitzcells, UseSchnitz, DisplayRange, Ellipses, SpotFilter, ZoomMode, GlobalZoomMode, ...
    ZoomRange, xForZoom, yForZoom, UseHistoneOverlay, HisOverlayFigAxes, HisPath1, HisPath2)

%PLOTFRAME Summary of this function goes here
%   Detailed explanation goes here


EllipseHandle=[];
EllipseHandleYellow=[];
EllipseHandleBlue=[];
EllipseHandleWhite=[];
EllipseHandleGreen=[];

%Get the coordinates of all the spots in this frame
[x,y,z]=SpotsXYZ(Spots{CurrentChannel}(CurrentFrame));
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

if isfield(FrameInfo, 'nc')
    set(Overlay,'Name',['Particle: ',num2str(CurrentParticle),'/',num2str(numParticles),...
        ', Frame: ',num2str(CurrentFrame),'/',num2str(numFrames),...
        ', Z: ',num2str(CurrentZ),'/',num2str(ZSlices),' nc: ', num2str(FrameInfo(CurrentFrame).nc),...
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
            Ellipses{CurrentFrame}(:,2)+1,'r',50, overlayAxes);
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
            Ellipses{CurrentFrame}(schnitzCellNo,2)+1,'b',50, overlayAxes);
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
                Ellipses{CurrentFrame}(NucleusIndex,2)+1,[],[], overlayAxes);
            set(EllipseHandleGreen,'Color','g')
            hold(overlayAxes,'off')
        else
            %('Error: Particle without an associated nucleus?')
        end



        %Show the daughter nuclei if applicable
        DaughterE=schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).E;
        DaughterD=schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).D;


        if DaughterE~=0
            SchnitzIndex=find(schnitzcells(DaughterE).frames==CurrentFrame);
            NucleusIndex=schnitzcells(DaughterE).cellno(SchnitzIndex);

            if ~isempty(NucleusIndex)
                hold(overlayAxes,'on')
                EllipseHandleWhite=[EllipseHandleWhite,ellipse(Ellipses{CurrentFrame}(NucleusIndex,3),...
                    Ellipses{CurrentFrame}(NucleusIndex,4),...
                    Ellipses{CurrentFrame}(NucleusIndex,5),...
                    Ellipses{CurrentFrame}(NucleusIndex,1)+1,...
                    Ellipses{CurrentFrame}(NucleusIndex,2)+1, [],[],overlayAxes)];

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
                    Ellipses{CurrentFrame}(NucleusIndex,2)+1,[],[],overlayAxes)];
                hold(overlayAxes,'off')
            else
                %('Error: Particle without an associated nucleus?')
            end
        end

        if ~isempty(EllipseHandleWhite)
            set(EllipseHandleWhite,'Color','w')
        end

        %Show the mother nucleus if applicable
        Mother=schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).P;

        if Mother~=0
            SchnitzIndex=find(schnitzcells(Mother).frames==CurrentFrame);
            NucleusIndex=schnitzcells(Mother).cellno(SchnitzIndex);

            if ~isempty(NucleusIndex)
                hold(overlayAxes,'on')
                EllipseHandleYellow=ellipse(Ellipses{CurrentFrame}(NucleusIndex,3),...
                    Ellipses{CurrentFrame}(NucleusIndex,4),...
                    Ellipses{CurrentFrame}(NucleusIndex,5),...
                    Ellipses{CurrentFrame}(NucleusIndex,1)+1,...
                    Ellipses{CurrentFrame}(NucleusIndex,2)+1,[],[],overlayAxes);
                set(EllipseHandleYellow,'Color','y')
                hold(overlayAxes,'off')
            else
                %('Error: Particle without an associated nucleus?')
            end
        end

    else
        if UseHistoneOverlay
            warning('This particle does not have an associated nucleus.');
        end
    end
end

if ApprovedParticles(CurrentParticle)==1
    set(Overlay,'Color','g')
    %         set(TraceFig,'Color','g')
elseif ApprovedParticles(CurrentParticle)==-1
    set(Overlay,'Color','r')
    %         set(TraceFig,'Color','r')
elseif ApprovedParticles(CurrentParticle)==2
    set(Overlay,'Color','y')
    %         set(TraceFig,'Color','y')
else
    set(Overlay,'Color','default')
    %         set(TraceFig,'Color','default')
end

%Show the particles that were under threshold 2.
if ShowThreshold2
    %Get the positions of all the spots in this frame
    [x2,y2]=SpotsXYZ(Spots{CurrentChannel}(CurrentFrame));
    %Filter those that were under threshold 2.
    CurrentSpotFilter=...
        ~logical(SpotFilter{CurrentChannel}(CurrentFrame,~isnan(SpotFilter{CurrentChannel}(CurrentFrame,:))));
    x2=x2(CurrentSpotFilter);
    y2=y2(CurrentSpotFilter);

    hold(overlayAxes,'on')
    plot(overlayAxes,x2,y2,'sr')
    hold(overlayAxes,'off')
end

if ZoomMode
    %Find the closest frame
    [~,MinIndex]=min((Particles{CurrentChannel}(CurrentParticle).Frame-CurrentFrame).^2);
    if length(MinIndex)>1
        MinIndex=MinIndex(1);
    end
    [xForZoom,yForZoom]=...
        SpotsXYZ(Spots{CurrentChannel}(Particles{CurrentChannel}(CurrentParticle).Frame(MinIndex)));

    xForZoom=xForZoom(Particles{CurrentChannel}(CurrentParticle).Index(MinIndex));
    yForZoom=yForZoom(Particles{CurrentChannel}(CurrentParticle).Index(MinIndex));

    try
        xlim(overlayAxes,[xForZoom-ZoomRange,xForZoom+ZoomRange])
        ylim(overlayAxes,[yForZoom-ZoomRange/2,yForZoom+ZoomRange/2])
    catch
        %something's outside the limits of the image
    end
end

if GlobalZoomMode
    xlim(overlayAxes,[xForZoom-ZoomRange,xForZoom+ZoomRange])
    ylim(overlayAxes,[yForZoom-ZoomRange/2,yForZoom+ZoomRange/2])
end

if UseHistoneOverlay
    try
        ImageHis=imread(HisPath1);
    catch %Had to do this for KITP
        ImageHis=imread(HisPath2);
    end

    if isempty(DisplayRange)
        HisOverlayImage=cat(3,mat2gray(ImageHis),mat2gray(Image),zeros(size(Image)));
    else
        HisOverlayImage=cat(3,mat2gray(ImageHis,double(DisplayRange)),mat2gray(Image),zeros(size(Image)));
    end
    
    imshow(HisOverlayImage,[],'Border','Tight','Parent',HisOverlayFigAxes)


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

    %         set(HisOverlayFigAxes,'Name',['Particle: ',num2str(CurrentParticle),'/',num2str(numParticles),...
    %             ', Frame: ',num2str(CurrentFrame),'/',num2str(numFrames),...
    %             ', Z: ',num2str(CurrentZ),'/',num2str(ZSlices),' nc: ', num2str(FrameInfo(CurrentFrame).nc),...
    %             ' Ch: ',num2str(CurrentChannel)])

    if ZoomMode || GlobalZoomMode
        try
            xlim(HisOverlayFigAxes,[xForZoom-ZoomRange,xForZoom+ZoomRange])
            ylim(HisOverlayFigAxes,[yForZoom-ZoomRange/2,yForZoom+ZoomRange/2])
        catch
            %something's outside the limits of the image
        end
    end

end

end

