function [CurrentFrame, ManualZFlag] = changeToFlaggedFrame(FrameRange, cptState)
%CHANGEFRAME 
%   Change the frame of a particle in CheckParticleTracking
    
    % find flagged frames
    CH = cptState.CurrentChannelIndex;
    CP = cptState.CurrentParticle;
    FrameOptions = ~cptState.Particles{CH}(CP).FrameApproved & ...
      ismember(cptState.Particles{CH}(CP).Frame,FrameRange);
    
    if ~isempty(FrameOptions)
        CandidateFrames = cptState.Particles{CH}(CP).Frame(FrameOptions);
        [~,mi] = min(abs(cptState.CurrentFrame-CandidateFrames));
        CurrentFrame = CandidateFrames(mi);
    else 
        CurrentFrame = cptState.CurrentFrame;    
    end
    
    ManualZFlag = 0;
    
end

