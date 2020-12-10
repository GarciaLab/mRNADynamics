Prefix='2020-10-30-hbBAC-MS2xHisRFP_MCP-GFP-27_5C-Anterior-Embryo1';
liveExperiment = LiveExperiment(Prefix);
schnitzcells = liveExperiment.getSchnitzcells;
particles = liveExperiment.getParticles;
spots = liveExperiment.getSpots;
FrameInfo = liveExperiment.getFrameInfo;
CompiledParticles = GetFluoZInfo(liveExperiment);
CompiledParticles = CompiledParticles{1};
%compiledspots = main01_compile_traces('hbBAC-MS2');
%%
pixelsize = liveExperiment.pixelSize_um;
nucleusDiameters = zeros(1, 6);
for nc=9:14
    nucleusDiameters(nc-8) = getDefaultParameters(FrameInfo,[['d', num2str(nc)]])/pixelsize; % in pixels
end
xDim = liveExperiment.xDim;
yDim = liveExperiment.yDim;
zDim = liveExperiment.zDim;
%%

if ~isfield('CompiledParticles', 'NuclearXYFlaggedFrames')
    for j = 1:length(CompiledParticles)
        CompiledParticles(j).NuclearXYFlaggedFrames = [];
    end
end
if ~isfield('CompiledParticles', 'SpotZFlaggedFrames')
    for j = 1:length(CompiledParticles)
        CompiledParticles(j).SpotZFlaggedFrames = [];
    end
end
if ~isfield('CompiledParticles', 'NoNucleusFlag')
    for j = 1:length(CompiledParticles)
        CompiledParticles(j).NoNucleusFlag = 0;
    end
end
for i=1:length(CompiledParticles)
    particle = CompiledParticles(i);
    if particle.Approved <= 0
        continue
    end
    NuclearXYFlaggedFrames = [];
    SpotZFlaggedFrames = [];
    NoNucleusFlag = 0;
    zFlag = 0;
    MinStartTime = 2; % minutes - need to come up with a formula for this.
    if ~isempty(particle.Nucleus)
        schnitz_index = particle.Nucleus;
        nc = schnitzcells(schnitz_index).cycle;
        diam = nucleusDiameters(nc-8);
        FrameVector = schnitzcells(schnitz_index).frames;
        TimeVector = horzcat(FrameInfo(FrameVector).Time)/60;
        TimeVector = TimeVector-min(TimeVector);
        xPosVector = schnitzcells(schnitz_index).cenx;
        yPosVector = schnitzcells(schnitz_index).ceny;
        MinFrame = find(TimeVector >= MinStartTime, 1);
        for f = MinFrame:length(FrameVector)
            if min([xPosVector(f), xDim-xPosVector(f),yPosVector(f),yDim-yPosVector(f)]) < diam
                NuclearXYFlaggedFrames = [NuclearXYFlaggedFrames, FrameVector(f)];
            end
            
        end
        for idx = 1:length(particle.Frame)
            f = particle.Frame(idx);
            if (particle.FluoZInfo(idx,1) == max(particle.FluoZInfo(idx,:))) | (particle.FluoZInfo(idx,zDim) == max(particle.FluoZInfo(idx,:))) 
                SpotZFlaggedFrames = [SpotZFlaggedFrames, f];
            end
        end

    else
        NoNucleusFlag = 1;
    end
    CompiledParticles(i).NoNucleusFlag = NoNucleusFlag;
    CompiledParticles(i).SpotZFlaggedFrames  = SpotZFlaggedFrames;
    CompiledParticles(i).NuclearXYFlaggedFrames = NuclearXYFlaggedFrames ;
end













