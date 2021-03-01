function cptState=removeCurrentFrame(cptState)

CurrentParticle = cptState.CurrentParticle;
Particles = cptState.Particles{cptState.CurrentChannelIndex};
CurrentFrame = cptState.CurrentFrame;
if isempty(find(Particles(CurrentParticle).Frame  == CurrentFrame, 1))
    disp("Can't remove current frame from trace because there is no current frame data.")
else
    %Separate the particle trace at the specified position
    fieldNames = fields(Particles);
    fieldNames(cellfun(@(x) strcmpi(x, 'Approved'), fieldNames) | cellfun(@(x) strcmpi(x, 'Schnitz'), fieldNames) | cellfun(@(x) strcmpi(x, 'Nucleus'), fieldNames) ) = [];
    
    FrameFilter=Particles(CurrentParticle).Frame ~= CurrentFrame;
    
    
    %Create a gap for a new particle
    NewParticles(1:CurrentParticle)=Particles(1:CurrentParticle);
    NewParticles(CurrentParticle+2:length(Particles)+1)=Particles(CurrentParticle+1:end);
    
    for f = 1:length(fieldNames)
        
        %Delete the information in the current particle
        
        NewParticles(CurrentParticle).(fieldNames{f}) = NewParticles(CurrentParticle).(fieldNames{f})(FrameFilter);
        NewParticles(CurrentParticle).Approved=0;
        
        %Move the information to the new particle
        
        NewParticles(CurrentParticle+1).(fieldNames{f})=Particles(CurrentParticle).(fieldNames{f})(~FrameFilter);
        NewParticles(CurrentParticle+1).Approved=0;
        NewParticles(CurrentParticle+1).Nucleus=[];
        if isfield(NewParticles, 'Schnitz')
            NewParticles(CurrentParticle+1).Schnitz = NaN;
        end
        
        
    end
    
    cptState.CurrentParticle = CurrentParticle;
    cptState.Particles{cptState.CurrentChannelIndex} = NewParticles;
    cptState.CurrentFrame = CurrentFrame;
    cptState.PreviousParticle = 1;
end