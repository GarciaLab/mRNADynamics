classdef CPTState < handle
    properties
        Spots
        Particles
        CurrentFrame
        PreviousFrame
        ManualZFlag
        CurrentZ
        ZSlices
        CurrentParticle
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
    end
    
    methods
        function this = CPTState(Spots, Particles, numberZSlices)
            this.Spots = Spots;
            this.Particles = Particles;
            this.CurrentFrame = 0;
            this.ManualZFlag = 0;
            this.PreviousFrame = this.CurrentFrame;
            this.ZSlices = numberZSlices + 2; %Note that the blank slices are included
            this.CurrentZ = round(this.ZSlices / 2);
            this.CurrentParticle = 1;
            this.CurrentChannel = 1;
            this.PreviousChannel = this.CurrentChannel;
            this.FrameIndicesToFit = 0; % index of the current particle that were used for fitting
            this.Coefficients = []; % coefficients of the fitted line
            this.fitApproved = 0; %JP: I don't think is this in use, will check.
            this.ZoomMode = 0;
            this.GlobalZoomMode = 0;
            this.xForZoom = 0;
            this.yForZoom = 0;
        end

        function numParticles = numParticles(this)
            numParticles = length(this.Particles{this.CurrentChannel});
        end
    end
end

