classdef CPTState < handle
    properties
        Spots
        Particles
        SpotFilter
        schnitzcells
        FrameInfo
        
        Frames
        CurrentFrame
        PreviousFrame
        
        ManualZFlag
        ZSlices
        CurrentZ
        
        CurrentParticleIndex
        CurrentParticle
        PreviousParticle
        lastParticle

        CurrentChannel
        PreviousChannel
        coatChannel
        
        FrameIndicesToFit
        Coefficients
        fitApproved
        
        ZoomMode
        GlobalZoomMode
        xForZoom
        yForZoom

        DisplayRange
        UseHistoneOverlay
        ImageHis

        no_clicking

        nameSuffix

	nWorkers
        plot3DGauss
    end
    
    methods
        function this = CPTState(Spots, Particles, SpotFilter, schnitzcells, FrameInfo, UseHistoneOverlay, nWorkers, plot3DGauss)
            this.Spots = Spots;
            this.Particles = Particles;
            this.SpotFilter = SpotFilter;
            this.schnitzcells = schnitzcells;
            this.FrameInfo = FrameInfo;
            
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
            this.PreviousChannel = this.CurrentChannel;
           
            this.FrameIndicesToFit = 0; % index of the current particle that were used for fitting
            this.Coefficients = []; % coefficients of the fitted line
            this.fitApproved = 0; %JP: I don't think is this in use, will check.
            
            this.ZoomMode = 0;
            this.GlobalZoomMode = 0;
            this.xForZoom = 0;
            this.yForZoom = 0;

            this.DisplayRange = [];
            this.UseHistoneOverlay = UseHistoneOverlay;

            this.no_clicking = false;

            this.nameSuffix = '';

            this.nWorkers = nWorkers;
            this.plot3DGauss = plot3DGauss;
        end

        function numParticles = numParticles(this)
            numParticles = length(this.Particles{this.CurrentChannel});
        end
    end
end

