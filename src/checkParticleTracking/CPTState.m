classdef CPTState < handle
    properties
        Spots
        Particles
        SpotFilter
        schnitzcells
        FrameInfo
        ImageMat

        Ellipses
        nucleiModified
        
        Frames
        CurrentFrame
        PreviousFrame
        
        ManualZFlag
        ZSlices
        CurrentZ
        CurrentZIndex
        
        CurrentParticleIndex
        CurrentParticle
        PreviousParticle
        lastParticle

        CurrentChannel
        CurrentChannelIndex
        PreviousChannel
        PreviousChannelIndex
        coatChannel
        
        FrameIndicesToFit
        Coefficients
        fitApproved
        lineFitted
        
        ZoomMode
        GlobalZoomMode
        xForZoom
        yForZoom

        DisplayRange
        DisplayRangeSpot
        UseHistoneOverlay
        ImageHis
        HideApprovedFlag

        nameSuffix

	    nWorkers
        plot3DGauss

        projectionMode
    end
    
    methods
        function this = CPTState(Spots, Particles, SpotFilter, schnitzcells, Ellipses,...
                FrameInfo, UseHistoneOverlay, nWorkers, plot3DGauss, projectionMode)
            
            this.Spots = Spots;
            this.Particles = Particles;
            this.SpotFilter = SpotFilter;
            this.schnitzcells = schnitzcells;
            this.FrameInfo = FrameInfo;
            this.ImageMat = [];
            
            this.Ellipses = Ellipses;
            this.nucleiModified = false;

            this.Frames = [];
            this.CurrentFrame = 0;
            this.PreviousFrame = this.CurrentFrame;
            
            this.ManualZFlag = 0;
            numberZSlices = this.FrameInfo(1).NumberSlices;
            this.ZSlices = numberZSlices + 2; %Note that the blank slices are included
            this.CurrentZ = round(this.ZSlices / 2);
        
            this.CurrentParticle = 1;
            this.PreviousParticle = 1;
            this.lastParticle = 0; %this gets flagged if there's a drop to one particle within the Particles structure.

            this.CurrentChannel = 1;
            this.CurrentChannelIndex = 1;
            this.PreviousChannel = this.CurrentChannel;
            this.PreviousChannelIndex = this.CurrentChannelIndex;

           
            this.FrameIndicesToFit = 0; % index of the current particle that were used for fitting
            this.Coefficients = []; % coefficients of the fitted line
            this.fitApproved = 0; %JP: I think functions should use this instead of calculating fitApproved on their own
            this.lineFitted = 0; % equals 1 if a line has been fitted

            this.ZoomMode = 0;
            this.GlobalZoomMode = 0;
            this.xForZoom = 0;
            this.yForZoom = 0;

            this.DisplayRange = [];
            this.DisplayRangeSpot = [];
            this.UseHistoneOverlay = UseHistoneOverlay;
            this.HideApprovedFlag = 0;

            this.nameSuffix = '';

            this.nWorkers = nWorkers;
            this.plot3DGauss = plot3DGauss;

            this.projectionMode = projectionMode;
        end

        function numParticles = numParticles(this)
            numParticles = length(this.Particles{this.CurrentChannelIndex});
        end

        function numValidFrames = numValidFrames(this)
            numValidFrames = length({this.Spots{1}.Fits});
        end

        function currentSpots = getCurrentChannelSpots(this)
            currentSpots = this.Spots{this.CurrentChannelIndex};
        end

        function currentFrameSpots = getCurrentFrameSpots(this)
            currentSpots = this.getCurrentChannelSpots();
            currentFrameSpots = currentSpots(this.CurrentFrame);
        end

        function currentParticles = getCurrentChannelParticles(this)
            currentParticles = this.Particles{this.CurrentChannelIndex};
        end

        function currentParticle = getCurrentParticle(this)
            currentParticles = this.getCurrentChannelParticles();
            currentParticle = currentParticles(this.CurrentParticle);
        end

        function currentParticleFit = getCurrentParticleFit(this)
            currentFrameSpots = this.getCurrentFrameSpots();
            currentParticleFit = currentFrameSpots.Fits(this.CurrentParticleIndex);
        end

        function currentXDoG = getCurrentXDoG(this)
            currentFit = this.getCurrentParticleFit();
            currentXDoG = double(currentFit.xDoG(this.CurrentZIndex));
        end

        function currentYDoG = getCurrentYDoG(this)
            currentFit = this.getCurrentParticleFit();
            currentYDoG = double(currentFit.yDoG(this.CurrentZIndex));
        end

        function currentXFit = getCurrentXFit(this)
            currentFit = this.getCurrentParticleFit();
            currentXFit = double(currentFit.xFit(this.CurrentZIndex));
        end

        function currentYFit = getCurrentYFit(this)
            currentFit = this.getCurrentParticleFit();
            currentYFit = double(currentFit.yFit(this.CurrentZIndex));
        end
    end
end

