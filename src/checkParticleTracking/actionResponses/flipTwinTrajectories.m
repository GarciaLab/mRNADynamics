function [Particles, CurrentParticle] = flipTwinTrajectories(cptState)
CurrentParticle = cptState.CurrentParticle;
CurrentFrame = cptState.CurrentFrame;
TwinParticle = cptState.TwinParticle;
Particles = cptState.Particles{cptState.CurrentChannelIndex};
Schnitz = Particles(cptState.CurrentParticle).Schnitz;
FrameFilter = (Particles(CurrentParticle).Frame < CurrentFrame);
if isempty(TwinParticle)
    if sum(FrameFilter) & sum(~FrameFilter)
        Particles=SeparateParticleTraces(CurrentParticle,CurrentFrame,Particles);
        Particles(CurrentParticle+1).Approved=0;
        Particles(CurrentParticle+1).Nucleus = Particles(CurrentParticle).Nucleus;
        Particles(CurrentParticle+1).Schnitz = Particles(CurrentParticle).Schnitz;
        TwinParticle = CurrentParticle + 1;
    else
        disp('Cannot create twin particle at current frame.')
    end
else
    FrameFilterTwin = (Particles(TwinParticle).Frame < CurrentFrame);
    if sum(FrameFilter) & sum(~FrameFilter) 
        Particles=SeparateParticleTraces(CurrentParticle,CurrentFrame,Particles);
        Particles(CurrentParticle+1).Approved=0;
        Particles(CurrentParticle+1).Nucleus = Particles(CurrentParticle).Nucleus;
        Particles(CurrentParticle+1).Schnitz = Particles(CurrentParticle).Schnitz;
        if TwinParticle > CurrentParticle
            TwinParticle = TwinParticle + 1;
        end
        if sum(FrameFilterTwin) & sum(~FrameFilterTwin)
             Particles=SeparateParticleTraces(TwinParticle,CurrentFrame,Particles);
             if CurrentParticle > TwinParticle
                CurrentParticle= CurrentParticle+ 1;
             end
            
              Particles=JoinParticleTraces(CurrentParticle,TwinParticle+1,Particles);
              if TwinParticle + 1 < CurrentParticle 
                  CurrentParticle = CurrentParticle -1;
              end
              Particles=JoinParticleTraces(TwinParticle,CurrentParticle+1,Particles);
              if CurrentParticle + 1 < TwinParticle
                  TwinParticle = TwinParticle -1;
              end
        elseif sum(FrameFilterTwin) & ~sum(~FrameFilterTwin)
            Particles=JoinParticleTraces(TwinParticle,CurrentParticle+1,Particles);
            if TwinParticle > CurrentParticle +1
                TwinParticle = TwinParticle -1;
            end
        elseif ~sum(FrameFilterTwin) & sum(~FrameFilterTwin)
            Particles=JoinParticleTraces(CurrentParticle,TwinParticle,Particles);
            if CurrentParticle > TwinParticle
                CurrentParticle = CurrentParticle -1;
            end  
            TwinParticle = CurrentParticle + 1;
        end
    elseif sum(FrameFilter) & ~sum(~FrameFilter) 
        if sum(FrameFilterTwin) & sum(~FrameFilterTwin)
             Particles=SeparateParticleTraces(TwinParticle,CurrentFrame,Particles);
             if CurrentParticle > TwinParticle
                CurrentParticle= CurrentParticle+ 1;
             end
             Particles=JoinParticleTraces(CurrentParticle,TwinParticle+1,Particles);
             if TwinParticle + 1 < CurrentParticle
                 CurrentParticle = CurrentParticle -1;
             end
        elseif ~sum(FrameFilterTwin) & sum(~FrameFilterTwin)
            Particles=JoinParticleTraces(CurrentParticle,TwinParticle,Particles);
            if TwinParticle < CurrentParticle
                CurrentParticle = CurrentParticle -1;
            end
            TwinParticle = [];
        else
            disp('Cannot switch twin particle at current frame.')
            
        end
    elseif ~sum(FrameFilter) & sum(~FrameFilter)
        if sum(FrameFilterTwin) & sum(~FrameFilterTwin)
             Particles=SeparateParticleTraces(TwinParticle,CurrentFrame,Particles);
             if CurrentParticle > TwinParticle
                CurrentParticle= CurrentParticle+ 1;
             end
             Particles=JoinParticleTraces(CurrentParticle,TwinParticle,Particles);
             if CurrentParticle > TwinParticle
                 CurrentParticle = CurrentParticle -1;
             end
        elseif sum(FrameFilterTwin) & ~sum(~FrameFilterTwin)
            Particles=JoinParticleTraces(CurrentParticle,TwinParticle,Particles);
            if TwinParticle < CurrentParticle
                CurrentParticle = CurrentParticle -1;
            end
            TwinParticle = [];
        else
            disp('Cannot switch twin particle at current frame.')
            
        end
        
    end
    


   
end
Particles(CurrentParticle).Nucleus = Schnitz;
Particles(CurrentParticle).Schnitz = Schnitz;
if ~isempty(TwinParticle)
    Particles(TwinParticle).Nucleus = Schnitz;
    Particles(TwinParticle).Schnitz = Schnitz;
end
cptState.Particles{cptState.CurrentChannelIndex} = Particles;
cptState.CurrentParticle = CurrentParticle;
end