function displayOverlays(overlayAxes, cptState, SpeedMode, ShowThreshold2,...
    Overlay, numFrames, UseSchnitz, ZoomRange, fish, multiAx,...
    HisOverlayFigAxes, ImageHisMat)

EllipseHandle = [];
EllipseHandleYellow = [];
EllipseHandleBlue = [];
EllipseHandleWhite = [];
EllipseHandleGreen = [];

ExperimentType2spot = strcmpi(cptState.liveExperiment.experimentType, '2spot');

if ~isempty(multiAx)
    multiView = true;
else
    multiView = false;
end

% Get the coordinates of all the spots in this frame
[x,y,z] = SpotsXYZ(cptState.getCurrentFrameSpots());
% Pull out the right particle if it exists in this frame
CurrentParticleIndex = cptState.getCurrentParticleIndex();
xTrace = x(CurrentParticleIndex);
yTrace = y(CurrentParticleIndex);


ApprovedParticles = [cptState.getCurrentChannelParticles().Approved];

% Approved particles
[xApproved, yApproved] = cptState.getApprovedParticles(x, y);


% Disapproved particles
[xDisapproved, yDisapproved] = cptState.getDisapprovedParticles(x, y);

% Non-flagged particles (these are particles that have not been processed)
[xNonFlagged, yNonFlagged] = cptState.getNonFlaggedParticles(x, y);

% G. Martini 2/13/21
% Attempting to speed up plotting by only keeping tracking of spots and
% ellipses within the zoom image when in zoom mode
if cptState.ZoomMode
    if ~isempty(xTrace)
        ApprovedInZoomFrame = (abs(xApproved-xTrace) <= ZoomRange) & (abs(yApproved-yTrace) <= ZoomRange);
        DisapprovedInZoomFrame = (abs(xDisapproved-xTrace) <= ZoomRange) & (abs(yDisapproved-yTrace) <= ZoomRange);
        NonFlaggedInZoomFrame = (abs(xNonFlagged-xTrace) <= ZoomRange) & (abs(yNonFlagged-yTrace) <= ZoomRange);
        
    else
        [~,MinIndex]=min((cptState.getCurrentParticle().Frame - cptState.CurrentFrame).^2);
        if length(MinIndex)>1
            MinIndex=MinIndex(1);
        end
        currentChannelSpots = cptState.getCurrentChannelSpots();
        [xForZoom,yForZoom]=...
            SpotsXYZ(currentChannelSpots(cptState.getCurrentParticle().Frame(MinIndex)));
        xForZoom = xForZoom(cptState.getCurrentParticle().Index(MinIndex));
        yForZoom = yForZoom(cptState.getCurrentParticle().Index(MinIndex));
        ApprovedInZoomFrame = (abs(xApproved-xForZoom) <= ZoomRange) & (abs(yApproved-yForZoom) <= ZoomRange);
        DisapprovedInZoomFrame = (abs(xDisapproved-xForZoom) <= ZoomRange) & (abs(yDisapproved-yForZoom) <= ZoomRange);
        NonFlaggedInZoomFrame = (abs(xNonFlagged-xForZoom) <= ZoomRange) & (abs(yNonFlagged-yForZoom) <= ZoomRange);
    end
    xApproved = xApproved(ApprovedInZoomFrame);
    yApproved = yApproved(ApprovedInZoomFrame);
    
    xDisapproved = xDisapproved(DisapprovedInZoomFrame);
    yDisapproved = yDisapproved(DisapprovedInZoomFrame);
    
    xNonFlagged = xNonFlagged(NonFlaggedInZoomFrame);
    yNonFlagged = yNonFlagged(NonFlaggedInZoomFrame);
end

% Show all particles in regular mode
if ~SpeedMode
    plot(overlayAxes,xNonFlagged,yNonFlagged,'ow')
    plot(overlayAxes,xApproved,yApproved,'ob')
    plot(overlayAxes,xDisapproved,yDisapproved,'^r')
end
% Always show current particle. this indicates the x-y center of the spot
% within the brightest z-slice and may differ from the position
% shown in the snippet image, which is centered at the position with
% the current z-slice.
plot(overlayAxes,xTrace,yTrace,'og')
if ExperimentType2spot
    TwinParticleIndex = cptState.getTwinParticleIndex();
    xTwin = x(TwinParticleIndex);
    yTwin = y(TwinParticleIndex);
    plot(overlayAxes,xTwin,yTwin,'sc')
end
hold(overlayAxes,'off')

if isfield(cptState.FrameInfo, 'nc')
    if ~cptState.UseTwinTraces
    set(Overlay,'Name',['Particle: ',num2str(cptState.CurrentParticle),'/',num2str(cptState.numParticles()),...
        ', Frame: ',num2str(cptState.CurrentFrame),'/',num2str(numFrames),...
        ', Z: ',num2str(cptState.CurrentZ),'/',num2str(cptState.ZSlices),' nc: ', num2str(cptState.FrameInfo(cptState.CurrentFrame).nc),...
        ', Ch: ',num2str(cptState.CurrentChannel)])
    else
      numParticles = length([cptState.Particles{cptState.CurrentChannelIndex}(:).Nucleus]);
      counter = 0;
      for cp = 1:cptState.CurrentParticle
          if ~isempty(cptState.Particles{cptState.CurrentChannelIndex}(cp).Nucleus)
              counter = counter + 1;
          end
      end
      
      Unchecked = length(find(~isnan([cptState.Particles{cptState.CurrentChannelIndex}(:).Schnitz]) & ...
          ([cptState.Particles{cptState.CurrentChannelIndex}(:).Approved] == 0)));
          
              
      set(Overlay,'Name',['Particle: ',num2str(counter),'(', num2str(cptState.CurrentParticle),')/',num2str(numParticles),...
          '(', num2str(Unchecked), ')',...
        ', Frame: ',num2str(cptState.CurrentFrame),'/',num2str(numFrames),...
        ', Z: ',num2str(cptState.CurrentZ),'/',num2str(cptState.ZSlices),' nc: ', num2str(cptState.FrameInfo(cptState.CurrentFrame).nc),...
        ', Ch: ',num2str(cptState.CurrentChannel)])  
    end
end
if UseSchnitz
    % Show all the nuclei in regular mode
    if ~SpeedMode
        EllipseHandle = notEllipseCPT(cptState, 'r', 10, overlayAxes, xTrace, yTrace, ZoomRange);
        
        % Show the ones that have been approved
        schnitzCellNo=[];
        
        for i=1:cptState.numParticles()
            currentChannelParticles = cptState.getCurrentChannelParticles();
            if currentChannelParticles(i).Approved == 1
                try
                    schnitzIndex = find((cptState.schnitzcells(currentChannelParticles(i).Nucleus).frames) == cptState.CurrentFrame);
                    schnitzCellNo = [schnitzCellNo,cptState.schnitzcells(currentChannelParticles(i).Nucleus).cellno(schnitzIndex)];
                catch
                    % can't identify the nucleus for this particle.
                end
            end
        end
        
        EllipseHandleBlue = notEllipseCellCPT(cptState, schnitzCellNo, 'b', 10, overlayAxes, xTrace, yTrace, ZoomRange);
    end
    
    % Show the corresponding nucleus
    if ~isempty(cptState.getCurrentParticle().Nucleus) &&...
            cptState.getCurrentParticle().Nucleus > 0
        SchnitzIndex = cptState.schnitzcells(...
            cptState.getCurrentParticle().Nucleus).frames == cptState.CurrentFrame;
        EllipsesIndex = cptState.schnitzcells(...
            cptState.getCurrentParticle().Nucleus).cellno(SchnitzIndex);
        
        if ~isempty(EllipsesIndex)
            EllipseHandleGreen = ellipseCellCPT(cptState, EllipsesIndex, 'g', 10, overlayAxes);
        end
        
        % Show the daughter nuclei if applicable
        [DaughterE, DaughterD, Mother] = cptState.getMotherDaughters();
        
        try
            if DaughterE~=0
                SchnitzIndex = find(cptState.schnitzcells(DaughterE).frames == cptState.CurrentFrame);
                EllipsesIndex = cptState.schnitzcells(DaughterE).cellno(SchnitzIndex);
                
                if ~isempty(EllipsesIndex)
                    EllipseHandleWhite = [EllipseHandleWhite,ellipseCellCPT(cptState, EllipsesIndex, 'w', 10, overlayAxes)];
                end
            end
            
            if DaughterD~=0
                SchnitzIndex = find(cptState.schnitzcells(DaughterD).frames == cptState.CurrentFrame);
                EllipsesIndex = cptState.schnitzcells(DaughterD).cellno(SchnitzIndex);
                
                if ~isempty(EllipsesIndex)
                    EllipseHandleWhite = [EllipseHandleWhite,ellipseCellCPT(cptState, EllipsesIndex, 'w', 10, overlayAxes)];
                end
            end
        catch
            warning('inconsistency found in schnitzcells. proceeding anyway.');
        end
        
        try
            %Show the mother nucleus if applicable
            if Mother~=0
                SchnitzIndex = cptState.schnitzcells(Mother).frames == cptState.CurrentFrame;
                EllipsesIndex = cptState.schnitzcells(Mother).cellno(SchnitzIndex);
                
                if ~isempty(EllipsesIndex)
                    EllipseHandleYellow=ellipseCellCPT(cptState, EllipsesIndex, 'y', 10, overlayAxes);
                end
            end
        catch
            warning('inconsistency found in schnitzcells. proceeding anyway.');
        end
        
    else
        if cptState.UseHistoneOverlay
            warning('This particle does not have an associated nucleus.');
        end
    end
end

if ApprovedParticles(cptState.CurrentParticle) == 1
    set(Overlay,'Color','g')
elseif ApprovedParticles(cptState.CurrentParticle) == -1
    set(Overlay,'Color','r')
elseif ApprovedParticles(cptState.CurrentParticle) == 2
    set(Overlay,'Color','y')
else
    set(Overlay,'Color','default')
end

%Show the particles that were under threshold 2.
if ShowThreshold2
    %Get the positions of all the spots in this frame
    [x2,y2]=SpotsXYZ(cptState.Spots{cptState.CurrentChannelIndex}(cptState.CurrentFrame));
    %Filter those that were under threshold 2.
    CurrentSpotFilter=...
        ~logical(cptState.SpotFilter{cptState.CurrentChannelIndex}(cptState.CurrentFrame,~isnan(cptState.SpotFilter{cptState.CurrentChannelIndex}(cptState.CurrentFrame,:))));
    x2=x2(CurrentSpotFilter);
    y2=y2(CurrentSpotFilter);
    
    
    if cptState.ZoomMode
        if ~isempty(xTrace)
            Threshold2ParticlesInZoomFrame = (abs(x2-xTrace) <= ZoomRange) & (abs(y2-yTrace) <= ZoomRange);
        else
            Threshold2ParticlesInZoomFrame = (abs(x2-xForZoom) <= ZoomRange) & (abs(y2-yForZoom) <= ZoomRange);
        end
        x2 = x2(Threshold2ParticlesInZoomFrame);
        y2 = y2(Threshold2ParticlesInZoomFrame);
    end
    
    hold(overlayAxes,'on')
    plot(overlayAxes,x2,y2,'sr')
    hold(overlayAxes,'off')
end

if cptState.ZoomMode
    % Find the closest frame
    [~,MinIndex]=min((cptState.getCurrentParticle().Frame - cptState.CurrentFrame).^2);
    if length(MinIndex)>1
        MinIndex=MinIndex(1);
    end
    currentChannelSpots = cptState.getCurrentChannelSpots();
    [cptState.xForZoom,cptState.yForZoom]=...
        SpotsXYZ(currentChannelSpots(cptState.getCurrentParticle().Frame(MinIndex)));
    
    cptState.xForZoom = cptState.xForZoom(cptState.getCurrentParticle().Index(MinIndex));
    cptState.yForZoom = cptState.yForZoom(cptState.getCurrentParticle().Index(MinIndex));
    
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
    
    imshow(HisOverlayImageMat,[],'Border','Tight','Parent',HisOverlayFigAxes);
    
    hold(HisOverlayFigAxes,'on')
    if ~SpeedMode
        plot(HisOverlayFigAxes,xNonFlagged,yNonFlagged,'ow')
        plot(HisOverlayFigAxes,xApproved,yApproved,'ob')
    end
    plot(HisOverlayFigAxes,xTrace,yTrace,'og')
    if ExperimentType2spot
        plot(HisOverlayFigAxes,xTwin,yTwin,'oc')
    end
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

end
