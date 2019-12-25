classdef CPTState < handle
    properties
        CurrentFrame
        PreviousFrame
        ManualZFlag
    end
    
    methods
        function this = CPTState(CurrentFrame, PreviousFrame, ManualZFlag)
            this.CurrentFrame = CurrentFrame;
            this.ManualZFlag = ManualZFlag;
            this.PreviousFrame = PreviousFrame;
        end
    end
end

