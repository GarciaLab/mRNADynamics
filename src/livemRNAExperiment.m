classdef liveExperiment
    %livemRNAExperiment Object to organize data related to a live imaging
    %experiment
    
    properties
        
        Prefix = '';
        isUnhealthy = false;
        
        anaphaseFrames = [0; 0; 0; 0; 0; 0];
        project = '';
        
        rawExportedDirectory = '';
        processedDirectory = '';
        resultsDirectory = '';
        
        preFolder = '';
        procFolder = '';
        resultsFolder = '';
        
        hasCompiledParticlesFile = false;
        hasSchnitzcellsFile = false;
        hasSpotsFile = false;
        hasParticlesFile = false;
        hasDoGs = false;
        hasRawStacks = false;
        hasMovieMatFile = false;
        hasHisMatFile = false;
        hasEllipsesFile = false;
        
        hasChannelsFile = false;
        hasAnaphaseFile = false;
        
    end
    
    methods
        
        
        %%Constructors
        
        function obj = livemRNAExperiment(Prefix)
            %livemRNAExperiment Construct an instance of this class
            
            obj.Prefix = Prefix;
            
            [rawDataPath, ProcPath, DropboxFolder, ~, PreProcPath,...
                ~, ~, ~, ~]= DetermineLocalFolders(obj.Prefix);
            
            obj.preFolder = [PreProcPath, filesep, Prefix, filesep];
            obj.procFolder = [ProcPath, filesep, Prefix, '_', filesep];
            obj.resultsFolder = [DropboxFolder, filesep, Prefix, filesep];
            
            [rawDir, procDir, resultsDir] = browseExperiment(obj.Prefix);
            obj.rawExportedDirectory = rawDir;
            obj.processedDirectory = procDir;
            obj.resultsDirectory = resultsDir;
            
            
            isUnhealthyFile = [DropboxFolder,filesep,obj.Prefix,filesep, 'isUnhealthy.mat'];
            if exist(isUnhealthyFile, 'file')
                load(isUnhealthyFile, 'isUnhealthy');
            else
                isUnhealthy = NaN;
            end
            obj.isUnhealthy = isUnhealthy;
            
            
            obj.anaphaseFrames = retrieveAnaphaseFrames(obj.Prefix);
            obj.hasAnaphaseFile=sum(contains(obj.resultsDirectory{:, 1}, 'anaphaseFrames'));
            
            obj.project = '';
            
            obj.hasCompiledParticlesFile = sum(contains(obj.resultsDirectory{:, 1}, 'CompiledParticles'));
            obj.hasSchnitzcellsFile = sum(contains(obj.resultsDirectory{:, 1}, '_lin'));
            obj.hasSpotsFile = sum(contains(obj.resultsDirectory{:, 1}, 'Spots'));
            obj.hasParticlesFile = sum(contains(obj.resultsDirectory{:, 1}, 'Particles'));
            obj.hasEllipsesFile = sum(contains(obj.resultsDirectory{:, 1}, 'Ellipses'));
            
            obj.hasDoGs = sum(contains(obj.processedDirectory{:, 1}, 'dog', 'IgnoreCase', true));
            obj.hasRawStacks = sum(contains(obj.rawExportedDirectory{:, 1}, 'stacks', 'IgnoreCase', true));
            obj.hasMovieMatFile = sum(contains(obj.rawExportedDirectory{:, 1}, 'movieMat', 'IgnoreCase', true));
            obj.hasHisMatFile = sum(contains(obj.rawExportedDirectory{:, 1}, 'hisMat', 'IgnoreCase', true));
            
            
            obj.hasChannelsFile =sum(contains(obj.resultsDirectory{:, 1}, 'Channels'));
            
            
        end
        
        
        
        %%
        
        
        %Methods
        
        
        function movieMat = getMovieMat(obj)
            
            movieMat = loadMovieMat(obj.Prefix);
            
        end
        
        function hisMat = getHisMat(obj)
            
            hisMat = loadHisMat(obj.Prefix);
            
        end
        
        function schnitzcells = getSchnitzcells(obj)
            
            schnitzcellsFile = [obj.resultsFolder, obj.Prefix, '_lin.mat'];
            if obj.hasSchnitzcellsFile
                schnitzcells = load(schnitzcellsFile);
            end
            
        end
        
        function Ellipses = getEllipses(obj)
            
            ellipsesFile = [obj.resultsFolder, 'Ellipses.mat'];
            if obj.hasEllipsesFile
                Ellipses = load(ellipsesFile);
            end
            
        end
        
    end
    
    
    
end
