function Particles=JoinParticleTraces(OriginalParticle,ClickedParticle,Particles)

%This function joins two particle traces and renumbers all particles in the
%Particles structure accordingly

%Transfer the information to the original particle
Particles(OriginalParticle).Frame=[Particles(OriginalParticle).Frame,Particles(ClickedParticle).Frame];
if isfield(Particles,'xPos')
    Particles(OriginalParticle).xPos=[Particles(OriginalParticle).xPos,Particles(ClickedParticle).xPos];
end
if isfield(Particles,'yPos')
    Particles(OriginalParticle).yPos=[Particles(OriginalParticle).yPos,Particles(ClickedParticle).yPos];
end
Particles(OriginalParticle).Index=[Particles(OriginalParticle).Index,Particles(ClickedParticle).Index];
Particles(OriginalParticle).Approved=0;
%Particles(OriginalParticle).nc=[Particles(OriginalParticle).nc,Particles(ClickedParticle).nc];
if isfield(Particles,'FrameApproved')
    Particles(OriginalParticle).FrameApproved=logical([Particles(OriginalParticle).FrameApproved,Particles(ClickedParticle).FrameApproved]);
end



%Now, get rid of the clicked particle
Particles=Particles([1:ClickedParticle-1,ClickedParticle+1:end]);

%Deals with the indexing changing because of the removal of
%the old particle.
 if ClickedParticle<OriginalParticle
     OriginalParticle=OriginalParticle-1;
 end

%Sort the frames within the particle. This is useful if we
%connected to a particle that came before.
[SortedFrame,Permutations]=sort(Particles(OriginalParticle).Frame);
Particles(OriginalParticle).Frame=Particles(OriginalParticle).Frame(Permutations);
try
    if isfield(Particles,'xPos')
        Particles(OriginalParticle).xPos=Particles(OriginalParticle).xPos(Permutations);
    end
    if isfield(Particles,'yPos')
        Particles(OriginalParticle).yPos=Particles(OriginalParticle).yPos(Permutations);
    end
end
Particles(OriginalParticle).Index=Particles(OriginalParticle).Index(Permutations);
% 9/4 EL: Added the if statement to remove error of the index (Permutations)
% exceeding matrix dimension. Did this always create an error when the
% frames of the original particle was not approved before?
if ~isempty(Particles(OriginalParticle).FrameApproved)
    Particles(OriginalParticle).FrameApproved=Particles(OriginalParticle).FrameApproved(Permutations);
end

end