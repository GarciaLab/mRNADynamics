classdef CPTState < handle
    properties
        Spots
        CurrentFrame
        PreviousFrame
        ManualZFlag
        CurrentZ
        ZSlices
    end
    
    methods
        function this = CPTState(Spots, CurrentFrame, PreviousFrame, ManualZFlag, numberZSlices)
            this.Spots = Spots;
            this.CurrentFrame = CurrentFrame;
            this.ManualZFlag = ManualZFlag;
            this.PreviousFrame = PreviousFrame;
            this.ZSlices = numberZSlices + 2; %Note that the blank slices are included
            this.CurrentZ = round(this.ZSlices / 2);
        end
    end
end

