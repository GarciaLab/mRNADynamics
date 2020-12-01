function [ParticleFrames,gxVec,gyVec,gzVec,gFitsCell,f3Vec,f3RawVec]=...
                            getGauss3DFitInfo(CurrentParticle,Particles,Spots)

% NL 5/2/2019

%First, get the different intensity values corresponding to this particle.
ParticleFrames = Particles(CurrentParticle).Frame;
gxVec = NaN(size(ParticleFrames));
gyVec = NaN(size(ParticleFrames));
gzVec = NaN(size(ParticleFrames));
f3Vec = NaN(size(ParticleFrames));
f3RawVec = NaN(size(ParticleFrames));

gFitsCell = cell(size(ParticleFrames));
for i=1:length(ParticleFrames)    
    spot = Spots(ParticleFrames(i)).Fits(Particles(CurrentParticle).Index(i));    
    if isfield(spot,'gauss3DIntensity') && ~isempty(spot.gauss3DIntensity)
        % position
        GaussPos = spot.GaussPos3D;
        gxVec(i) = GaussPos(1);
        gyVec(i) = GaussPos(2);
        gzVec(i) = GaussPos(3);    
        f3Vec(i) = spot.gauss3DIntensity;
        f3RawVec(i) = spot.gauss3DIntensityRaw;        
        % general fit info
        gFitsCell{i} = spot.SpotFits3D;  
    end
end