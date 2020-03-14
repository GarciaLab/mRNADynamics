classdef livemRNAExperiment
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
        
        hasCompiledParticlesFile = false;
        hasSchnitzcellsFile = false;
        hasSpotsFile = false;
        hasParticlesFile = false;
        hasDoGs = false;
        hasRawStacks = false;
        hasMovieMatFile = false;
        hasHisMatFile = false;
        
        hasChannelsFile = false;
        hasAnaphaseFile = false;
        
    end
    
    methods
        
               
        %%Constructors
        
        function obj = livemRNAExperiment(Prefix)
            %livemRNAExperiment Construct an instance of this class
            
            obj.Prefix = Prefix;
            
            [~, ~, DropboxFolder, ~, ~] = DetermineLocalFolders(obj.Prefix);

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
            obj.hasSchnitzcellsFile = sum(contains(obj.resultsDirectory{:, 1}, 'schnitzcells'));
            obj.hasSpotsFile = sum(contains(obj.resultsDirectory{:, 1}, 'Spots'));
            obj.hasParticlesFile = sum(contains(obj.resultsDirectory{:, 1}, 'Particles'));
            
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
        
        
    end
    
    
    
end

