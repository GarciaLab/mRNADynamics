classdef liveExperiment
    %livemRNAExperiment object to organize data related to a live imaging
    %experiment
    
    properties
        
        Prefix = '';
        preFolder = '';
        procFolder = '';
        resultsFolder = '';
        MLFolder = '';
        project = '';
        Channels = {};
        spotChannelIndices = [];
        inputChannelIndices = [];
        nuclearChannels = {};
        
        isUnhealthy = false;
        
        anaphaseFrames = [0; 0; 0; 0; 0; 0];
        
        
    end
    
    properties (Hidden)
        
        
        rawFolder = '';
        
        rawExportedDirectory = '';
        processedDirectory = '';
        resultsDirectory = '';
        
        userRawFolder = '';
        userPreFolder = '';
        userProcFolder = '';
        userResultsFolder = '';
        
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
        
        Channel1 = '';
        Channel2 = '';
        Channel3 = '';
        
        MS2CodePath = '';
        
        
        
    end
    
    methods
        
        
        %%Constructors
        
        function obj = liveExperiment(Prefix)
            %livemRNAExperiment Construct an instance of this class
            
            obj.Prefix = Prefix;
            
            %caching the results of this function since
            %its csv2cell call is the
            %most time-consuming part of
            %project initialization and the output
            %is not memory intensive
            memoizedDetermineLocalFolders = memoize(@DetermineLocalFolders);
            
            [rawPath, ProcPath, DropboxFolder, obj.MS2CodePath, PreProcPath,...
                ~, ~, ~, movieDatabase]= memoizedDetermineLocalFolders(obj.Prefix);
            
            
            obj.userPreFolder = PreProcPath;
            obj.userProcFolder = ProcPath;
            obj.userResultsFolder = DropboxFolder;
            obj.userRawFolder =  rawPath;
            Subfolder = [obj.Prefix(1:10),filesep,obj.Prefix(12:length(obj.Prefix))];
            
            obj.rawFolder = strcat(obj.userRawFolder,filesep,Subfolder);
            obj.preFolder = [PreProcPath, filesep, Prefix, filesep];
            obj.procFolder = [ProcPath, filesep, Prefix, '_', filesep];
            obj.resultsFolder = [DropboxFolder, filesep, Prefix, filesep];
            obj.MLFolder = [DropboxFolder, filesep, 'training_data_and_classifiers', filesep];
            
            %             [rawDir, procDir, resultsDir] = browseExperiment(obj.Prefix);
            %             obj.rawExportedDirectory = rawDir;
            %             obj.processedDirectory = procDir;
            %             obj.resultsDirectory = resultsDir;
            
            
            isUnhealthyFile = [DropboxFolder,filesep,obj.Prefix,filesep, 'isUnhealthy.mat'];
            if exist(isUnhealthyFile, 'file')
                load(isUnhealthyFile, 'isUnhealthy');
            else isUnhealthy = NaN; end
            
            obj.isUnhealthy = isUnhealthy;
            
            obj.project = '';
            
            obj.hasCompiledParticlesFile = exist([obj.resultsFolder, 'CompiledParticles.mat'] , 'file');
            obj.hasSchnitzcellsFile = exist([obj.resultsFolder,Prefix, '_lin.mat'] , 'file');
            obj.hasSpotsFile = exist([obj.resultsFolder, 'Spots.mat'] , 'file');
            obj.hasParticlesFile = exist([obj.resultsFolder, 'Particles.mat'] , 'file');
            obj.hasEllipsesFile = exist([obj.resultsFolder, 'Ellipses.mat'] , 'file');
            obj.hasChannelsFile =exist([obj.resultsFolder, 'Channels.mat'] , 'file');
            obj.hasAnaphaseFile=exist([obj.resultsFolder, 'anaphaseFrames.mat'] , 'file');
            
            obj.hasDoGs = exist([obj.procFolder, 'dogs'], 'dir');
            
            obj.hasRawStacks = exist([obj.preFolder, 'stacks'], 'dir');
            obj.hasMovieMatFile = exist([obj.preFolder, filesep, Prefix, '_movieMatCh1.mat'], 'file');
            obj.hasHisMatFile = exist([obj.preFolder, filesep, Prefix, '_hisMat.mat'], 'file');
            
            [~, obj.experimentType, obj.experimentAxis, ~, ~, obj.APResolution,...
                Channel1, Channel2,~, ~,  ~, ~, ~,...
                ~, ~, ~, ~, ~, ~, ~, Channel3,~,~, ~, obj.DVResolution]...
                = getExperimentDataFromMovieDatabase(Prefix, movieDatabase, obj.userResultsFolder);
            
            obj.Channels = {Channel1{1}, Channel2{1}, Channel3{1}};
            obj.Channel1 = Channel1{1};
            obj.Channel2 = Channel2{1};
            obj.Channel3 = Channel3{1};
            
            try
                FrameInfo = getFrameInfo(obj);
                [obj.xDim, obj.yDim, obj.pixelSize_nm, obj.zStep_um, obj.snippetSize_px,...
                    obj.nFrames, obj.zDim, obj.nDigits] = getFrameInfoParams(FrameInfo);
                obj.pixelSize_um = obj.pixelSize_nm/1000;
            catch, warning('FrameInfo not found.'); end
            
            
            obj.inputChannelIndices = find(contains(obj.Channels, 'input', 'IgnoreCase', true));
            
            obj.spotChannelIndices = getCoatChannel(Channel1, Channel2, Channel3);
            
            obj.anaphaseFrames = retrieveAnaphaseFrames(obj.Prefix, obj.userResultsFolder);
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
        
        function Channels = getChannels(this)
            % JP: This is filtering out NaN, empty Strings and Channels
            % beginning with 'His'. As discussed with armando, it's
            % probably better to get the channels from PreProcess data. If
            % we do that, this method would be the place to implement that
            filteredChannels = this.Channels(~strcmp(this.Channels, 'NaN'));
            filteredChannels = filteredChannels(~strcmp(filteredChannels, ''));
            filteredChannels = filteredChannels(~startsWith(filteredChannels, 'His'));
            Channels = filteredChannels;
        end
        
        function out = getMovieMat(this)
            
            persistent preTifDir;
            if isempty(preTifDir) ||...
                    ~isequal( length(preTifDir), this.nFrames)
            preTifDir = dir([this.preFolder, '*_ch0*.tif']);
            end
            
            exportedChannels = [];
            % find what channels were exported           
            for k = 1:5  %i don't see channel number going beyond 6 any time soon. 
                exportedChannels(k) =  any(cellfun(@(x) contains(x, ['_ch0',num2str(k)]),...
                    {preTifDir.name}));
            end
            channelsToRead = find(exportedChannels);
     
            % this is for backwards compatibility, exported tiffs used to be one per z slice.
            haveTifStacks = any(cellfun(@(x) ~contains(x, '_z'),...
                {preTifDir.name})); 
            
            
            persistent movieMat;
            %load movie only if it hasn't been loaded or if we've switched
            %Prefixes (determined by num frames)
            if isempty(movieMat) ||...
                    ~isequal( size(movieMat, 4), this.nFrames)
                
                if this.hasMovieMatFile
%                 if false
                    movieMat = loadMovieMat(this.Prefix);
                elseif ~haveTifStacks
                    movieMat = makeMovieMats(this.Prefix, [], [], [],...
                        'loadHis', false, 'makeMovie', true, 'loadMovie', false);
                elseif haveTifStacks
                    nPadding = 2;
                    
                    movieMat = zeros(this.yDim, this.xDim,...
                        this.zDim+nPadding, this.nFrames, length(channelsToRead), 'uint8'); % y x z t ch
                    
                    chIndex = 0;
                   
                    for ch = channelsToRead
                        chIndex = chIndex + 1;
                        
                        preChDir = preTifDir(cellfun(@(x) contains(x,...
                            ['ch0', num2str(ch)]), {preTifDir.name}));
                        
                        for f = 1:this.nFrames
                            movieMat(:, :, :, f, chIndex) =...
                                imreadStack2([this.preFolder, filesep, preChDir(f).name],...
                                this.yDim, this.xDim, this.zDim+nPadding);     
                        end
                        
                    end
                else
                    error('can''t load movie.')
                end
                
            end
            out = movieMat;
            
        end
        
        function out = getHisMat(obj)
            
            if exist([obj.preFolder, filesep,obj.Prefix, '-His.tif'], 'file')
                haveHisTifStack = true;
            end
            persistent hisMat;
            %load histone movie only if it hasn't been loaded or if we've switched
            %Prefixes (determined by num frames)
            if isempty(hisMat) || ~isequal( size(hisMat, 3), obj.nFrames)
                if obj.hasHisMatFile
                    hisMat = loadHisMat(obj.Prefix);
                elseif haveHisTifStack
                    
                    hisMat = imreadStack2([obj.preFolder, filesep,obj.Prefix, '-His.tif'], obj.yDim, obj.xDim, obj.nFrames);
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