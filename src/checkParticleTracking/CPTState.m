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
    end
    
    methods
        function this = CPTState(Spots, Particles, CurrentFrame, ManualZFlag, numberZSlices, CurrentParticle, CurrentChannel)
            this.Spots = Spots;
            this.Particles = Particles;
            this.CurrentFrame = CurrentFrame;
            this.ManualZFlag = ManualZFlag;
            this.PreviousFrame = CurrentFrame;
            this.ZSlices = numberZSlices + 2; %Note that the blank slices are included
            this.CurrentZ = round(this.ZSlices / 2);
            this.CurrentParticle = CurrentParticle;
            this.CurrentChannel = CurrentChannel;
            this.PreviousChannel = CurrentChannel;
        end

        function numParticles = numParticles(this)
            numParticles = length(this.Particles{this.CurrentChannel});
        end
    end
end

