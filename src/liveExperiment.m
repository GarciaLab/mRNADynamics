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
        spotChannels = [];
        
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
        
        zStep_um = 0;
        nDigits = 0;
        pixelSize_nm = 0;
        pixelSize_um = 0;
        snippetSize_px = 0;
        
        nc9 = 0;
        nc10 = 0;
        nc11 = 0;
        nc12 = 0;
        nc13 = 0;
        nc14 = 0;
        
        experimentType = '';
        experimentAxis = '';
        APResolution = '';
        DVResolution = '';
        
        
        
    end
    
    methods
        
        
        %%Constructors
        
        function obj = liveExperiment(Prefix)
            %livemRNAExperiment Construct an instance of this class
            
            obj.Prefix = Prefix;
            
            [~, ProcPath, DropboxFolder, ~, PreProcPath,...
                ~, ~, ~, movieDatabase]= DetermineLocalFolders(obj.Prefix);
            
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
            
            [~, obj.experimentType, obj.experimentAxis, ~, ~, obj.APResolution,...
                Channel1, Channel2,~, ~,  ~, ~, ~,...
                ~, ~, ~, ~, ~, ~, ~, Channel3,~,~, ~, obj.DVResolution]...
                = getExperimentDataFromMovieDatabase(Prefix, movieDatabase);
                        
            obj.Channels = {Channel1{1}, Channel2{1}, Channel3{1}};
            
            try
                FrameInfo = getFrameInfo(obj);
                [obj.xDim, obj.yDim, obj.pixelSize_nm, obj.zStep_um, obj.snippetSize_px,...
                    obj.nFrames, obj.zDim, obj.nDigits] = getFrameInfoParams(FrameInfo);
                obj.pixelSize_um = obj.pixelSize_nm/1000;
            catch
                warning('FrameInfo not found.')
            end
            
            obj.spotChannels = getCoatChannel(Channel1, Channel2, Channel3);
            
            obj.anaphaseFrames = retrieveAnaphaseFrames(obj.Prefix);
            obj.hasAnaphaseFile=sum(contains(obj.resultsDirectory{:, 1}, 'anaphaseFrames'));
            if numel(obj.anaphaseFrames) < 6
                obj.anaphaseFrames = vertcat(obj.anaphaseFrames, nan(6-numel(obj.anaphaseFrames), 1));
            end
            obj.nc9 = obj.anaphaseFrames(1);
            obj.nc10 = obj.anaphaseFrames(2);
            obj.nc11 = obj.anaphaseFrames(3);
            obj.nc12 = obj.anaphaseFrames(4);
            obj.nc13 = obj.anaphaseFrames(5);
            obj.nc14 = obj.anaphaseFrames(6);
            
            
            
        end
        
        
        
        %%
        
        
        %Methods
        
        
        function out = getMovieMat(obj)
            
            persistent movieMat;
            %load movie only if it hasn't been loaded or if we've switched
            %Prefixes (determined by num frames)
            if isempty(movieMat) || ~isequal( size(movieMat, 4), obj.nFrames)
                if obj.hasMovieMatFile
                    movieMat = loadMovieMat(obj.Prefix);
                else
                    movieMat = makeMovieMats(obj.Prefix, [], [], [], 'loadHis', false, 'makeMovie', true, 'loadMovie', false);
                end
            end
            out = movieMat;
            
        end
        
        function out = getHisMat(obj)
            
            persistent hisMat;
            %load histone movie only if it hasn't been loaded or if we've switched
            %Prefixes (determined by num frames)
            if isempty(hisMat) || ~isequal( size(hisMat, 3), obj.nFrames)
                if obj.hasHisMatFile
                    hisMat = loadHisMat(obj.Prefix);
                else
                    [~,hisMat] = makeMovieMats(obj.Prefix, [], [], [], 'loadMovie', false,  'loadHis', false, 'makeMovie', false, 'makeHis', true);
                end
            end
            out = hisMat;
            
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
        
        function APDetection = getAPDetection(obj)
            
            APDetectionFile = [obj.resultsFolder, 'APDetection.mat'];
            if exist(APDetectionFile, 'file')
                APDetection =  load(APDetectionFile);
            end
            
        end
        
    end
    
    methods(Static)
        
        function movieDatabase = getMovieDatabase
            [~, ~, ~, ~, ~,...
                ~, ~, ~, movieDatabase] = DetermineLocalFolders;
        end
        
    end
    
    
end