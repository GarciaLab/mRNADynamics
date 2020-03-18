classdef liveExperiment
    %livemRNAExperiment Object to organize data related to a live imaging
    %experiment
    
    properties
        
        Prefix = '';
        preFolder = '';
        procFolder = '';
        resultsFolder = '';
        MLFolder = '';
        project = '';
        Channels = {};
        spotChannel = [];
        
        isUnhealthy = false;
        
        anaphaseFrames = [0; 0; 0; 0; 0; 0];
        
        
    end
    
    properties (Hidden)
        
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
        hasEllipsesFile = false;
        
        hasChannelsFile = false;
        hasAnaphaseFile = false;
        
        yDim = 0;
        xDim = 0;
        zDim = 0;
        nFrames = 0;
        
        zStep = 0;
        nDigits = 0;
        pixelSize_nm = 0;
        
        
    end
    
    methods
        
        
        %%Constructors
        
        function obj = liveExperiment(Prefix)
            %livemRNAExperiment Construct an instance of this class
            
            obj.Prefix = Prefix;
            
            [rawDataPath, ProcPath, DropboxFolder, ~, PreProcPath,...
                ~, ~, ~, ~]= DetermineLocalFolders(obj.Prefix);
            
            obj.preFolder = [PreProcPath, filesep, Prefix, filesep];
            obj.procFolder = [ProcPath, filesep, Prefix, '_', filesep];
            obj.resultsFolder = [DropboxFolder, filesep, Prefix, filesep];
            obj.MLFolder = [DropboxFolder, filesep, 'training_data_and_classifiers', filesep];
            
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
            
            
            [~,~,~,~, ~,...
                ~, ~, ~,~,~,~,...
                ~, ~, ~, movieDatabase]...
                = readMovieDatabase(Prefix);
            
            [~, ~, ~, ~, ~, ~,...
                Channel1, Channel2,~, ~,  ~, ~, ~,...
                ~, ~, ~, ~, ~, ~, ~, Channel3,~,~, ~, ~]...
                = getExperimentDataFromMovieDatabase(Prefix, movieDatabase);
            
            obj.Channels = {Channel1{1}, Channel2{1}, Channel3{1}};
            
            FrameInfo = getFrameInfo(obj);
            
            [xSize, ySize, pixelSize_nm, zStep, snippet_size,...
                nFrames, nSlices, nDigits] = getFrameInfoParams(FrameInfo);
            
            obj.xDim = xSize;
            obj.yDim = ySize;
            obj.zDim = nSlices;
            obj.nFrames = nFrames;
            obj.zStep = zStep;
            obj.nDigits = nDigits;
            obj.pixelSize_nm = pixelSize_nm;
            
            obj.spotChannel = getCoatChannel(Channel1, Channel2, Channel3);
            
            
        end
        
        
        
        %%
        
        
        %Methods
        
        
        function movieMat = getMovieMat(obj)
            
            movieMat = loadMovieMat(obj.Prefix);
            
        end
        
        function hisMat = getHisMat(obj)
            
            hisMat = loadHisMat(obj.Prefix);
            
        end
        
        function FrameInfo = getFrameInfo(obj)
            
            load([obj.resultsFolder,filesep,'FrameInfo.mat'], 'FrameInfo');
            
        end
        
        function schnitzcells = getSchnitzcells(obj)
            
            schnitzcellsFile = [obj.resultsFolder, obj.Prefix, '_lin.mat'];
            if obj.hasSchnitzcellsFile
                load(schnitzcellsFile, 'schnitzcells');
            end
            
        end
        
        function Ellipses = getEllipses(obj)
            
            ellipsesFile = [obj.resultsFolder, 'Ellipses.mat'];
            if obj.hasEllipsesFile
                load(ellipsesFile, 'Ellipses');
            end
            
        end
        
        function CompiledParticles = getCompiledParticles(obj)
            
            compiledParticlesFile = [obj.resultsFolder, 'CompiledParticles.mat'];
            if obj.hasCompiledParticlesFile
                CompiledParticles = load(compiledParticlesFile);
            end
            
        end
        
        function Spots = getSpots(obj)
            
            spotsFile = [obj.resultsFolder, 'Spots.mat'];
            if obj.hasSpotsFile
               load(spotsFile, 'Spots');
            end
            
        end
        
        function [Particles, SpotFilter] = getParticles(obj)
            
            particlesFile = [obj.resultsFolder, 'Particles.mat'];
            if obj.hasParticlesFile
                load(particlesFile, 'Particles', 'SpotFilter');
            end
            
        end
        
        
    end
    
    
    
end