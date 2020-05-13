classdef LiveExperiment
    %LivemRNAExperiment object to organize data related to a live imaging
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
        inputChannels = [];
        nuclearChannels = {};
        
        isUnhealthy = false;
        
        anaphaseFrames (:, 1) uint16 = [0; 0; 0; 0; 0; 0];
        
        
    end
    
    properties (Access = private)
                
    end
    
    properties (Hidden)
        
        
        rawFolder = '';
        
        rawExportedDirectory = '';
        processedDirectory = '';
        resultsDirectory = '';
        
        experimentFolder = '';
        
        userRawFolder = '';
        userPreFolder = '';
        userProcFolder = '';
        userResultsFolder = '';
        userLivemRNAFolder = '';
        userExperimentsFolder = '';
        
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
        
        function this = LiveExperiment(Prefix)
            %livemRNAExperiment Construct an instance of this class
            
            this.Prefix = Prefix;
            
            %caching the results of this function since
            %its csv2cell call is the
            %most time-consuming part of
            %project initialization and the output
            %is not memory intensive
            %             DetermineLocalFolders = memoize(@DetermineLocalFolders);
            
            [this.userRawFolder, this.userProcFolder, this.userResultsFolder,...
                this.MS2CodePath, this.userPreFolder,...
                ~, ~, ~, movieDatabase]= DetermineLocalFolders(this.Prefix);

            dateString = this.Prefix(1:10);
            experimentName = this.Prefix(12:length(this.Prefix));
            rawSubFolder = [dateString,filesep,experimentName];
            
            liveFolderIndex = strfind(lower(this.userPreFolder), lower('LivemRNA'))...
                +  length('livemRNA') - 1;
            
            this.userLivemRNAFolder = this.userPreFolder(1:liveFolderIndex); 
            
            this.userExperimentsFolder = [this.userLivemRNAFolder, filesep, 'Experiments']; 
            
           
            if exist( this.userExperimentsFolder, 'dir')
                this.experimentFolder = [this.userExperimentsFolder, filesep, Prefix];
            end
            
            this = setExperimentFolders(this);
            
            if isempty(this.preFolder) || isempty(this.procFolder) || isempty(this.rawFolder)
                
                this.rawFolder = strcat(this.userRawFolder,filesep,rawSubFolder);
                this.preFolder = [this.userPreFolder, filesep, this.Prefix, filesep];
                this.procFolder = [this.userProcFolder, filesep, this.Prefix, '_', filesep];
                
            end
            
            this.resultsFolder = [this.userResultsFolder, filesep, this.Prefix, filesep];
            this.MLFolder = [this.userResultsFolder, filesep, 'training_data_and_classifiers', filesep];
            
            
            isUnhealthyFile = [this.userResultsFolder,filesep,this.Prefix,filesep, 'isUnhealthy.mat'];
            if exist(isUnhealthyFile, 'file')
                load(isUnhealthyFile, 'isUnhealthy');
            else, isUnhealthy = NaN;
            end
            
            this.isUnhealthy = isUnhealthy;
            
            this.project = '';
            
            this.hasCompiledParticlesFile = exist([this.resultsFolder, 'CompiledParticles.mat'] , 'file');
            this.hasSchnitzcellsFile = exist([this.resultsFolder,this.Prefix, '_lin.mat'] , 'file');
            this.hasSpotsFile = exist([this.resultsFolder, 'Spots.mat'] , 'file');
            this.hasParticlesFile = exist([this.resultsFolder, 'Particles.mat'] , 'file');
            this.hasEllipsesFile = exist([this.resultsFolder, 'Ellipses.mat'] , 'file');
            this.hasChannelsFile =exist([this.resultsFolder, 'Channels.mat'] , 'file');
            this.hasAnaphaseFile=exist([this.resultsFolder, 'anaphaseFrames.mat'] , 'file');
            
            this.hasDoGs = exist([this.procFolder, 'dogs'], 'dir');
            
            this.hasRawStacks = exist([this.preFolder, 'stacks'], 'dir');
            this.hasMovieMatFile = exist([this.preFolder, filesep, Prefix, '_movieMatCh1.mat'], 'file');
            this.hasHisMatFile = exist([this.preFolder, filesep, Prefix, '_hisMat.mat'], 'file');
            
            [~, this.experimentType, this.experimentAxis, ~, ~, this.APResolution,...
                Channel1, Channel2,~, ~,  ~, ~, ~,...
                ~, ~, ~, ~, ~, ~, ~, Channel3,~,~, ~, this.DVResolution]...
                = getExperimentDataFromMovieDatabase(this.Prefix, movieDatabase, this.userResultsFolder);
            
            this.Channels = {Channel1{1}, Channel2{1}, Channel3{1}};
            this.Channel1 = Channel1{1};
            this.Channel2 = Channel2{1};
            this.Channel3 = Channel3{1};
            
            try
                [this.xDim, this.yDim, this.pixelSize_nm, this.zStep_um, this.snippetSize_px,...
                    this.nFrames, this.zDim, this.nDigits] = getFrameInfoParams(getFrameInfo(this));
                this.pixelSize_um = this.pixelSize_nm/1000;
            catch %nothing to see here
            end
            
            
            this.inputChannels = find(contains(this.Channels, 'input', 'IgnoreCase', true));
            
            this.spotChannels = getCoatChannel(Channel1, Channel2, Channel3);
            
            this.anaphaseFrames = retrieveAnaphaseFrames(this.Prefix, this.userResultsFolder);
            
            this.nc9 = this.anaphaseFrames(1);
            this.nc10 = this.anaphaseFrames(2);
            this.nc11 = this.anaphaseFrames(3);
            this.nc12 = this.anaphaseFrames(4);
            this.nc13 = this.anaphaseFrames(5);
            this.nc14 = this.anaphaseFrames(6);
            
            
            
        end
        
        
        
        %%
        
        
        %Methods
        
        function this = setExperimentFolders(this)

                expFolder = [this.userExperimentsFolder, filesep, this.Prefix];
                
                if isfolder(expFolder)

                    this.preFolder = [expFolder, filesep, 'PreProcessedData'];
                    this.procFolder = [expFolder, filesep, 'ProcessedData'];
                    this.rawFolder = [expFolder, filesep, 'RawDynamicsData'];

                end
                
        end
        
        
        function Channels = getChannels(this)
            % JP: This is filtering out NaN, empty Strings and Channels
            % beginning with 'His'. As discussed with armando, it's
            % probably better to get the channels from PreProcess data. If
            % we do that, this method would be the place to implement that
            filteredChannels = this.Channels(~strcmpi(this.Channels, 'NaN'));
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
                exportedChannels(k) =  any(contains(...
                    string({preTifDir.name}), ['_ch0',num2str(k)]));
            end
            channelsToRead = find(exportedChannels);
            
            % this is for backwards compatibility,
            %exported tiffs used to be one per z slice.
            haveTifStacks = any(~contains(...
                string({preTifDir.name}), '_z'));
            
            
            persistent movieMat;
            %load movie only if it hasn't been loaded or if we've switched
            %Prefixes (determined by num frames)
            if isempty(movieMat) ||...
                    ~isequal( size(movieMat, 4), this.nFrames)
                if haveTifStacks
                    movieMat = makeMovieMatFromTifStacks(this, preTifDir, channelsToRead);
                elseif this.hasMovieMatFile
                    %load in .mat file
                    movieMat = loadMovieMat(this.Prefix);
                elseif ~haveTifStacks
                    %load movie from individual tif slices
                    movieMat = makeMovieMats(this.Prefix, [], [], [],...
                        'loadHis', false, 'makeMovie', true, 'loadMovie', false);
                else
                    error('can''t load movie.')
                end
                
            end
            out = movieMat;
            
            if max(movieMat(:)) < 255
                movieMat = uint8(movieMat);
            end
            
        end
        
        
        
        function movieMat = makeMovieMatFromTifStacks(this, preTifDir, channelsToRead)
            
            nPadding = 2;
            
            moviePrecision = 'uint16';
            movieMat = zeros(this.yDim, this.xDim,...
                this.zDim+nPadding, this.nFrames,...
                length(channelsToRead), moviePrecision); % y x z t ch
            
            chIndex = 0;
            
            for ch = channelsToRead
                
                chIndex = chIndex + 1;
                
                preChDir = preTifDir( ...
                    contains(...
                    string({preTifDir.name}), ['ch0', num2str(ch)]) &...
                    ~contains(string({preTifDir.name}), '_z') );
                
                %making these temporary variables to avoid passing all
                %of
                %liveExperiment to the parforloop
                this_nFrames = this.nFrames;
                this_yDim = this.yDim;
                this_preFolder = this.preFolder;
                this_xDim = this.xDim;
                this_zDim = this.zDim;
                
                parfor f = 1:this_nFrames
                    movieMat(:, :, :, f, chIndex) =...
                        imreadStack2([this_preFolder, filesep, preChDir(f).name],...
                        this_yDim, this_xDim, this_zDim+nPadding);
                end
                
                
            end
            
        end
        
        function out = getHisMat(this)
            
            if exist([this.preFolder, filesep,this.Prefix, '-His.tif'], 'file')
                haveHisTifStack = true;
            else
                haveHisTifStack = false;
            end
            persistent hisMat;
            %load histone movie only if it hasn't been loaded or if we've switched
            %Prefixes (determined by num frames)
            if isempty(hisMat) || ~isequal( size(hisMat, 3), this.nFrames)
                if haveHisTifStack
                    %load in sequential tif stacks
                    hisMat = imreadStack2([this.preFolder, filesep,...
                        this.Prefix, '-His.tif'], this.yDim, this.xDim, this.nFrames);
                
                %deprecated filetype. here for backwards compatibility
                elseif this.hasHisMatFile
                    %load up .mat histone file
                    hisMat = loadHisMat(this.Prefix);
                    
                %deprecated filetype. here for backwards compatibility
                else
                    %load in individual tif slices
                    [~,hisMat] = makeMovieMats(this.Prefix, [], [], [],...
                        'loadMovie', false,  'loadHis', false,...
                        'makeMovie', false, 'makeHis', true);
                    
                end
            end
            out = hisMat;
            
            if max(hisMat(:)) < 255
                hisMat = uint8(hisMat);
            end
            
        end
        
        function FrameInfo = getFrameInfo(this)
            
            load([this.resultsFolder,filesep,'FrameInfo.mat'], 'FrameInfo');
            
        end
        
        function schnitzcells = getSchnitzcells(this)
            
            schnitzcellsFile = [this.resultsFolder, this.Prefix, '_lin.mat'];
            if this.hasSchnitzcellsFile
                load(schnitzcellsFile, 'schnitzcells');
            end
            
        end
        
        function Ellipses = getEllipses(this)
            
            ellipsesFile = [this.resultsFolder, 'Ellipses.mat'];
            if this.hasEllipsesFile
                load(ellipsesFile, 'Ellipses');
            end
            
        end
        
        function CompiledParticles = getCompiledParticles(this)
            
            compiledParticlesFile = [this.resultsFolder, 'CompiledParticles.mat'];
            if this.hasCompiledParticlesFile
                CompiledParticles = load(compiledParticlesFile);
            end
            
        end
        
        function Spots = getSpots(this)
            
            spotsFile = [this.resultsFolder, 'Spots.mat'];
            if this.hasSpotsFile
                load(spotsFile, 'Spots');
            end
            
        end
        
        function [Particles, SpotFilter] = getParticles(this)
            
            particlesFile = [this.resultsFolder, 'Particles.mat'];
            if this.hasParticlesFile
                load(particlesFile, 'Particles', 'SpotFilter');
            end
            
        end
        
        function APDetection = getAPDetection(this)
            
            APDetectionFile = [this.resultsFolder, 'APDetection.mat'];
            if exist(APDetectionFile, 'file')
                APDetection =  load(APDetectionFile);
            end
            
        end
        
        
        function anaphaseFrames = getAnaphaseFrames(this)
            
           anaphaseFrames = this.anaphaseFrames;
            
        end
        
    end
    
    methods(Static)
        
        function movieDatabase = getMovieDatabase
            [~, ~, ~, ~, ~,...
                ~, ~, ~, movieDatabase] = DetermineLocalFolders;
        end
        
    end
    
    
end