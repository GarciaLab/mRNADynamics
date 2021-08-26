function [Particles, Spots, NChannels] = loadParticlesIfExists(DropboxFolder)
	if exist([DropboxFolder,filesep,Prefix,filesep,'Particles.mat'], 'file')
	    load([DropboxFolder,filesep,Prefix,filesep,'Particles.mat'], 'Particles', 'SpotFilter')
	    load([DropboxFolder,filesep,Prefix,filesep,'Spots.mat'], 'Spots')
	    
	    % Create the particle array. This is done so that we can support multiple
	    % channels. Also figure out the number of channels
	    if iscell(Particles)
	        NChannels = length(Particles);
	    else
	        Particles = {Particles};
	        Spots = {Spots};
	        NChannels = 1;
	    end
	    
	    % Now, get the particle positions (if they're not there already).
	    for ChN = 1:NChannels
	        Particles = addPositionsToParticles(Particles, Spots, ChN);
	    end
	    
	    if isfield(Particles{ChN}, 'DVpos')
	        warning('Particles.mat already has DV positions stored. They will be rewritten')
	    end
	    if isfield(Particles{ChN}, 'APpos')
	        warning('Particles.mat already has AP positions stored. They will be rewritten')
	    end
	else
		Particles = []
		Spots = []
		NChannels = []
	    warning('No Particles.mat found. Just updating APDetection.mat')
	end

end