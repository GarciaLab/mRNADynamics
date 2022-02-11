function [CompiledParticles] =...
    GetFluoZInfo(liveExperiment, CompiledParticles, Particles, Spots) %plotTraceSettings, noSpline)
% This function uses the total intensity mask to calculate the particles
% intensity and subtracts the background obtained from the fit.
if ~exist('CompiledParticles', 'var')
    CompiledParticlesAll = getCompiledParticles(liveExperiment);
    CompiledParticles = CompiledParticlesAll.CompiledParticles;
end
if ~exist('Particles', 'var')
    Particles = getParticles(liveExperiment);
end
if ~exist('Spots', 'var')
    Spots = getSpots(liveExperiment);
end
%% 

% First, get the different intensity values corresponding to this particle.
for ChN=1:length(liveExperiment.spotChannels)
if ~isfield('CompiledParticles', 'FluoZInfo')
    for j=1:length(CompiledParticles{ChN})
        CompiledParticles{ChN}(j).FluoZInfo = zeros(length(CompiledParticles{ChN}(j).Fluo), liveExperiment.zDim);
    end
end

defaultArea = 109; %109 pixels is the default area when the pixels are assumed to be 212nm x 212 nm AR 9/3/18

warning('off', 'MATLAB:rankDeficientMatrix'); %suppress the spline fitting warnings
for j=1:length(CompiledParticles{ChN})
    CurrentParticle = CompiledParticles{ChN}(j).OriginalParticle;
    
    for i=1:length(Particles{ChN}(CurrentParticle).Frame)
        
        spot = Spots(Particles{ChN}(CurrentParticle).Frame(i)).Fits(Particles{ChN}(CurrentParticle).Index(i));
        
        %Determine the brightest Z plane of this particle
        for zIndex = 1:length(spot.z)
            z = spot.z(zIndex);
            CompiledParticles{ChN}(j).FluoZInfo(i,z) = spot.FixedAreaIntensity(zIndex);
        end
    end
    
end
%% 
end

DropboxFolder = liveExperiment.userResultsFolder;
CompiledParticlesFile = [DropboxFolder,filesep,liveExperiment.Prefix,filesep,'CompiledParticles.mat'];

try
    save(CompiledParticlesFile, 'CompiledParticles','-append','-v6');
catch
    %save as 7.3 only if we really need to
    save(CompiledParticlesFile, 'CompiledParticles','-append','-v7.3', '-nocompression');
end

CompiledParticlesToken = now;
save([DropboxFolder,filesep,liveExperiment.Prefix,filesep,'CompiledParticlesToken.mat'],'CompiledParticlesToken', '-v6')

disp('CompiledParticles.mat saved.');

