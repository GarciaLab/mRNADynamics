function Particles=JoinParticleTraces(OriginalParticle,ClickedParticle,Particles)

%This function joins two particle traces and renumbers all particles in the
%Particles structure accordingly

%Transfer the information to the original particle
Particles(OriginalParticle).Frame=[Particles(OriginalParticle).Frame,Particles(ClickedParticle).Frame];
Particles(OriginalParticle).Index=[Particles(OriginalParticle).Index,Particles(ClickedParticle).Index];
Particles(OriginalParticle).Approved=0;
%Particles(OriginalParticle).nc=[Particles(OriginalParticle).nc,Particles(ClickedParticle).nc];
if isfield(Particles,'FrameApproved')
    Particles(OriginalParticle).FrameApproved=[Particles(OriginalParticle).FrameApproved,Particles(ClickedParticle).FrameApproved];
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
Particles(OriginalParticle).Index=Particles(OriginalParticle).Index(Permutations);
Particles(OriginalParticle).FrameApproved=Particles(OriginalParticle).FrameApproved(Permutations);    
