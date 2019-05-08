function [particle_frames,gx_vec,gy_vec,gz_vec,g_fits_cell]=...
    getGauss3DFitInfo(CurrentParticle,Particles,Spots)

% NL 5/2/2019

%First, get the different intensity values corresponding to this particle.
particle_frames = Particles(CurrentParticle).Frame;
gx_vec = NaN(size(particle_frames));
gy_vec = NaN(size(particle_frames));
gz_vec = NaN(size(particle_frames));
g_fits_cell = cell(size(particle_frames));
for i=1:length(particle_frames)    
    spot = Spots(particle_frames(i)).Fits(Particles(CurrentParticle).Index(i));    
    if isfield(spot,'gauss3DIntensity') 
        % position
        GaussPos = spot.GaussPos;
        gx_vec(i) = GaussPos(1);
        gy_vec(i) = GaussPos(2);
        gz_vec(i) = GaussPos(3);        
        % general fit info
        g_fits_cell{i} = spot.fits3D;  
    end
end