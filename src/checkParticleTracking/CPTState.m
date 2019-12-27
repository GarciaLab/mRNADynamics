classdef CPTState < handle
    properties
        CurrentFrame
        PreviousFrame
        ManualZFlag
        CurrentZ
    end
    
    methods
        function this = CPTState(CurrentFrame, PreviousFrame, ManualZFlag, CurrentZ)
            this.CurrentFrame = CurrentFrame;
            this.ManualZFlag = ManualZFlag;
            this.PreviousFrame = PreviousFrame;
            this.CurrentZ = CurrentZ;
        end
    end
end

