for i=1:length(Particles)
    index=find(abs(Particles(i).Frame-33*ones(1,length(Particles(i).Frame)))<4);
    if length(index)>1
    Particles=SeparateParticleTraces(i,index(1),Particles);
    end
end