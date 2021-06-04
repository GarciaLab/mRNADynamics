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
        
        anaphaseFrames = [0; 0; 0; 0; 0; 0];

        experimentName = '';

        hasCustomMasks = false;
        
        
    end
    
    properties (Access = private)
        
        preLoadMovie = false;
        
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
        
        fileMode = '';
        
        experimentType = '';
        experimentAxis = '';
        APResolution = '';
        DVResolution = '';
        
        Channel1 = '';
        Channel2 = '';
        Channel3 = '';
        
        MS2CodePath = '';
        
        %% Added by GM 9/29/20
        
        Temp_set = 0;
        Temp_obs = 0;
        

        
        %% 
        
        
    end
    
    methods
        
        
        %%Constructors
        
        function this = LiveExperiment(Prefix, preLoadMovie)
            %livemRNAExperiment Construct an instance of this class
            
            this.Prefix = Prefix;
            
            if nargin > 1
                this.preLoadMovie = preLoadMovie;
            end
            
            [this.userRawFolder, this.userProcFolder, this.userResultsFolder,...
                this.MS2CodePath, this.userPreFolder,...
                ~, ~, ~, movieDatabase]= DetermineLocalFolders(this.Prefix);
            
            dateString = this.Prefix(1:10);
            this.experimentName = this.Prefix(12:length(this.Prefix));
            rawSubFolder = [dateString,filesep,this.experimentName];
            
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
            
            
            try
                isUnhealthyFile = [this.userResultsFolder,filesep,this.Prefix,filesep, 'isUnhealthy.mat'];
                if exist(isUnhealthyFile, 'file')
                    load(isUnhealthyFile, 'isUnhealthy');
                else, isUnhealthy = NaN;
                end
                this.isUnhealthy = isUnhealthy;
            catch        
                %not important for most things if this fails. 
            end
            
            
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

            this.hasCustomMasks = isfolder([this.rawFolder, filesep, 'masks']);

            
            [~, this.experimentType, this.experimentAxis, ~, ~, this.APResolution,...
                Channel1, Channel2,~, ~,  ~, ~, ~,...
                ~, ~, ~, ~, ~, ~, ~, Channel3,~,~, ~, this.DVResolution, this.Temp_set, this.Temp_obs]...
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
            
            try
                evalc('[~, this.fileMode] = DetermineFileMode(this.rawFolder)');    %Using evalc to supress displays to the command window from the function DetermineFileMode
            catch
                %let's not do anything and let the caller handle this
                %exception. 
            end
            
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
        
        function nuclearMask = getNuclearMask(this, frameNumber, zIndex)
            T = sprintf('%02d', frameNumber);
            Z = sprintf('%03d', zIndex);
            maskFilePath = [this.rawFolder, filesep, 'masks', filesep, this.experimentName, '_T', T, '_Z', Z, '.tif'];
            
            disp(['Reading nuclear mask from file: ', maskFilePath]);

            nuclearMask = imread(maskFilePath);
            
            nuclearMask = cast(nuclearMask, 'logical');
        end
        
        function out = getMovieMat(this)
            
            %we're going to check if this is a new prefix by verifying the
            %exact equality of frame times in frameinfo. i don't think
            %these should ever be precisely identical in different movies
            persistent FrameInfo_movie;
            persistent movieMat;
            
            tempInfo = load([this.resultsFolder,filesep,'FrameInfo.mat'], 'FrameInfo');
            
            isNewMovie = isempty(FrameInfo_movie) ||...
                length([tempInfo.FrameInfo.Time]) ~= length([FrameInfo_movie.Time]) || ...
                any([tempInfo.FrameInfo.Time] ~= [FrameInfo_movie.Time]) ||...
                size(movieMat, 4) ~= this.nFrames;
            
            persistent preTifDir;
            if isempty(preTifDir) || isNewMovie
                FrameInfo_movie = tempInfo.FrameInfo;
                preTifDir = dir([this.preFolder, '*_ch0*.tif']);
            end
                        
            %just return an empty array if we can't load the movie.
            %leave the handling to the caller, presumably by enabling
            %sequential file loading.
            if ~haveSufficientMemory(preTifDir) || ~this.preLoadMovie

                out = [];
                return;
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
            
            
            %load movie only if it hasn't been loaded or if we've switched
            %Prefixes (determined by num frames) or if the old FrameInfo doesn't match
            %the new FrameInfo
            if isempty(movieMat) || isNewMovie
                
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
            
            %let's reduce the memory footprint of the movie if we can
%             if max(movieMat(:)) < 255
%                 movieMat = uint8(movieMat);
%             end
            movieMat = double(movieMat);
            
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
                
                for f = 1:this_nFrames
                    movieMat(:, :, :, f, chIndex) =...
                        imreadStack2([this_preFolder, filesep, preChDir(f).name],...
                        this_yDim, this_xDim, this_zDim+nPadding);
                end
                
                
            end
            
        end
        
        function out = getHisMat(this)
            
            hisFile = [this.preFolder, filesep,this.Prefix, '-His.tif'];
            if exist(hisFile, 'file')
                haveHisTifStack = true;
            else
                haveHisTifStack = false;
            end
            
                       
            persistent hisMat;
            persistent FrameInfo_His;
            tempInfo = load([this.resultsFolder,filesep,'FrameInfo.mat'], 'FrameInfo');
            
            isNewMovie = isempty(FrameInfo_His) ||...
                size([tempInfo.FrameInfo.Time],2) ~= size([FrameInfo_His.Time],2) || ...
                any([tempInfo.FrameInfo.Time] ~= [FrameInfo_His.Time]) ||...
                size(hisMat, 3) ~= this.nFrames;
            
            
            %load histone movie only if it hasn't been loaded or if we've switched
            %Prefixes (determined by num frames)
            if isempty(hisMat) || isNewMovie
                
                FrameInfo_His = tempInfo.FrameInfo;
                
                if haveHisTifStack
                    %load in sequential tif stacks
                    hisMat = imreadStack2([this.preFolder, filesep,...
                        this.Prefix, '-His.tif'], this.yDim, this.xDim, this.nFrames);
                    
                    %deprecated filetype. here for backwards compatibility
                elseif this.hasHisMatFile
                    %load up .mat histone file
                    hisMat = loadHisMat(this.Prefix);
                else
                    %load in individual tif slices
                    [~,hisMat] = makeMovieMats(this.Prefix, [], [], [],...
                        'loadMovie', false,  'loadHis', false,...
                        'makeMovie', false, 'makeHis', true);
                    
                end
                
            end
            out = hisMat;
            
            %let's reduce the memory footprint of the movie if we can
%             if max(hisMat(:)) < 255
%                 hisMat = uint8(hisMat);
%             end
            hisMat = double(hisMat);
            
        end
        
        %this function retrieves or generates a max projection of the full
        %movie
        function out = getMaxMat(this)
            
            persistent FrameInfo_max;
            persistent maxMat;
            
            tempInfo = load([this.resultsFolder,filesep,'FrameInfo.mat'], 'FrameInfo');
            
            isNewMovie = isempty(FrameInfo_max) ||...
                length([tempInfo.FrameInfo.Time]) ~= length([FrameInfo_max.Time]) || ...
                any([tempInfo.FrameInfo.Time] ~= [FrameInfo_max.Time]) ||...
                size(maxMat, 3) ~= this.nFrames;
            
            if isNewMovie
                FrameInfo_max = tempInfo.FrameInfo;
                movieMat = getMovieMat(this);
                if ~isempty(movieMat)
                    maxMat = max(movieMat(:,:,:,:,:),[],3);
                else
                    maxMat = [];
                end
            end
            out = squeeze(maxMat);
            out = double(out);
            
        end
        
        function out = getMovieFrame(this, frame, channel)
            
            stackFile = [this.preFolder, filesep, this.Prefix, '_',...
                iIndex(frame, this.nDigits),...
                '_ch', iIndex(channel, 2), '.tif'];
            
            if exist(stackFile, 'file')
                out = imreadStack(stackFile);
            elseif exist([this.preFolder, filesep, this.Prefix, '_',...
                iIndex(frame, this.nDigits),...
                '_z', iIndex(1, 2),...
                '_ch', iIndex(channel, 2), '.tif'], 'file')
                
                out = zeros(this.yDim, this.xDim, this.zDim);
                for z = 1:this.zDim
                    out(:, :, z) = getMovieSlice(this, frame, channel, z);
                end
            end
            out = double(out);
        end
        
        function out = getMovieSlice(this, frame, channel, slice)
            
            imFile = [this.preFolder, filesep, this.Prefix, '_',...
                iIndex(frame, this.nDigits),...
                '_z', iIndex(slice, 2),...
                '_ch', iIndex(channel, 2), '.tif'];
            
            if exist(imFile, 'file')
                out = imread(imFile);
            elseif exist([this.preFolder, filesep, this.Prefix, '_',...
                iIndex(frame, this.nDigits),...
                '_ch', iIndex(channel, 2), '.tif'], 'file')
            
                imStack = getMovieFrame(this, frame, channel);
                out = imStack(:, :, slice);
            end
            out = double(out);
            
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
            else
                warning('No Ellipses structure found for this LiveExperiment.')
                Ellipses = {};
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
        
        function rawDataFile = getRawDataFile(this)
            dashPositions = strfind(this.Prefix,'-');
            rawFileName= this.Prefix((dashPositions(3)+1):end);
            
            if strcmpi(this.fileMode, 'LIFExport')
                rawDataFile = [this.rawFolder, filesep, rawFileName, '.lif'];
            else
                error('getRawDataFile function only supports LIF files right now')
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
