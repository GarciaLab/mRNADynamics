
Prefix='2020-10-07-hbBAC-MS2xHisRFP_MCP-GFP-17_5C-Anterior-Embryo1';
lE = LiveExperiment(Prefix);
FrameInfo = lE.getFrameInfo;
PixelSize = lE.pixelSize_um;
Particles = lE.getParticles;
schnitzcells = lE.getSchnitzcells;
Distances = [];
Velocities = [];
Xpos = [];
Ypos  = [];
Frames = [];
ParticleIDs = [];
NCs = [];
j = 1;
for p=1:length(Particles)
    Particle = Particles(p);
    if (Particle.Approved ~= -1) & (length(Particle.Frame) > 1)
        if ~isempty(Particle.Nucleus)
            NuclearCycle = schnitzcells(Particle.Nucleus).cycle;
        else
            NuclearCycle = NaN;
        end
        for f = 1:(length(Particle.Frame)-1)
            Frames(j) = Particle.Frame(f+1);
            Xpos(j) = Particle.xPos(f+1);
            Ypos(j) = Particle.yPos(f+1);
            ParticleIDs(j) = p;
            
            NCs(j) = NuclearCycle;
            Distances(j) = sqrt((Particle.xPos(f+1)-Particle.xPos(f))^2 + (Particle.yPos(f+1)-Particle.yPos(f))^2);
            Velocities(j) = Distances(j)*PixelSize/(FrameInfo(f+1).Time-FrameInfo(f).Time);
            j = j+ 1;
        end
    end
end