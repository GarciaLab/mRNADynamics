function CompiledParticles = addCycle(Prefix, varargin)

if ischar(Prefix)
    
%What type of experiment are we dealing with? Get this out of MovieDatabase
[rawDataPath,ProcPath,DropboxFolder,MS2CodePath, PreProcPath,...
    rawDataFolder, Prefix, ExperimentType,Channel1,Channel2,OutputFolder,...
    Channel3, spotChannels, MovieDataBaseFolder, movieDatabase]...
    = readMovieDatabase(Prefix, optionalResults);


load([DropboxFolder,filesep,Prefix,filesep,'CompiledParticles.mat']);
liveExperiment = LiveExperiment(Prefix);

anaphaseFrames = getAnaphaseFrames(liveExperiment); 

% ncFrames = [zeros(1, 8), nc9, nc10, nc11, nc12, nc13, nc14];
ncFrames = [zeros(1, 8), anaphaseFrames];
else
    %optionally just give the function compiledparticles straight
    CompiledParticles = Prefix;
    ncFrames = varargin{1};
end
    


if iscell(CompiledParticles)
    nCh = length(CompiledParticles);
else
    nCh = 1;
end


for ch = 1:nCh
    
    for p = 1:length(CompiledParticles{ch})
        
        inds = find(CompiledParticles{ch}(p).Frame(1) < ncFrames); %find anaphases occuring after first frame of trace
        if isempty(inds)
            CompiledParticles{ch}(p).cycle = 13; %edge case for 14th cycle
        else 
            CompiledParticles{ch}(p).cycle = inds(1) - 1;
        end
    end
    
end

CompiledParticles{ch}(p).nc = CompiledParticles{ch}(p).cycle;



end