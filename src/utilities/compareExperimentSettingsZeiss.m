% function [ComparedSettings,RawSettings] = compareExperimentSettingsZeiss(DynamicsResultsPath,RawDataPath,DataType,varargin)
%
% DESCRIPTION
% Compares the experimental settings (aka metadata) for multiple datasets
% and determines whether or not they are the identical. Only works for
% Leica datasets at the moment.
%
% PARAMETERS
% DynamicsResultsPath: Full or relative path to DynamicsResults folder, where the
%             DataStatus.xlsx file is saved
% RawDataPath: Full path to RawDynamicsData folder, where all the sets to
%              be compared are saved
% DataType: This is a cell array of char variable(s) which are the exact 
%           name(s) of the tab(s) in DataStatus.xlsx that you wish to 
%           analyze.
%           E.g. {'DataType'} for only one tab
%           E.g. {'DataType1', 'DataType2', 'DataType3'} for multiple
%           tabsst
% 
%
% OPTIONS
% 'CompareDecimals': Should be followed by an integer between 1 and 4,
%                    indicating the number of decimal places you want
%                    compared to determine if the settings match
% 'LaserTolerance': Should be followed by a number between 0 and 1
%                   indicating the percent tolerance that is allowed when 
%                   comparing the laser power between datasets. The default
%                   tolerance is 0.10.
% 'JustReady': Only compares the settings of Prefixes that have "Ready" or
%              "ApproveAll" in the Compile
%
% OUTPUT
% ComparedSettings: Returns a structure containing booleans indicating 
%                   which settings match.                 
% RawSettings: Returns a structure containing all the microsocpe settings 
%           for the datasets found in the indicated DataStatus tab,   
%           allowing for easy visual comparison.
% 
%
% Author (contact): Meghan Turner (meghan_turner@berkeley.edu)
% Created: 2019-04-30
% Last Updated: 2019-10-31
%
% Documented by: Meghan Turner (meghan_turner@berkeley.edu)

function [ComparedSettings,RawSettings] = compareExperimentSettingsZeiss(DynamicsResultsPath,RawDataPath,DataType,varargin)

%Tolerances - number of decimals to which to round off any setting to consider settings the same
CompareDecimals = 2;
pinholeDecimals = 4; %This needs to be fixed

LaserTolerance = 0.1;   % Percent tolerance (as opposed to decimal tolerance above)
justReady = false;      

for i= 1:length(varargin)
    if strcmpi(varargin{i},'CompareDecimals')
        CompareDecimals = varargin{i+1};
    elseif strcmpi(varargin{i},'LaserTolerance')
        LaserTolerance = varargin{i+1};
        if LaserTolerance < 0 || LaserTolerance > 1
            error('LaserTolerance must be between 0 and 1.')
        end
    
    elseif strcmpi(varargin{i}, 'JustReady')
        justReady = true;
    end
end

% Look in DataStatus.xlsx in folderPath and find the tab given by
% the input variable dataType.
if isfile([DynamicsResultsPath,'\','DataStatus.xlsx']) || ...
                        isfile([DynamicsResultsPath,'\','DataStatus.csv'])  % NB: isfile doesn't let you use * to check for multiple file extensions
    DataStatusDir = dir([DynamicsResultsPath,filesep,'DataStatus.*']);
    
    if isa(DataType, 'char') || ~isa(DataType, 'cell')
        error('Wrong data type for DataType. It must be a cell array of single quote string(s), e.g. {''DataType''} or {''DataType1'', ''DataType2''}')
    end
    
    [~,DataTab]=xlsread([DynamicsResultsPath,filesep,DataStatusDir(1).name],DataType{1,1});
    for i = 2:numel(DataType)
        [~,TempDataTab]=xlsread([DynamicsResultsPath,filesep,DataStatusDir(1).name],DataType{1,i});
        numNewSets = size(TempDataTab,2) - 1;
        DataTab(:,end+1:end+numNewSets) = TempDataTab(:,2:end);
    end
        
else
    error(['DataStatus.xlsx does not exist in the indicated folder path, ' ...
            DynamicsResultsPath, '.'])
end

%Find and load the different prefixes
PrefixRow = find(strcmpi(DataTab(:,1),'Prefix:'));
SizeDataTab = size(DataTab);
%Check if we only need to load the datasets marked as "ready"
if justReady
    CompileRow=find(strcmpi(DataTab(:,1),'AnalyzeLiveData Compile Particles')|...
                    strcmpi(DataTab(:,1),'CompileParticles')|...
                    strcmpi(DataTab(:,1),'CompileNuclearProtein'));
    prefixFilter=find(strcmpi(DataTab(CompileRow,:),'READY')|strcmpi(DataTab(CompileRow,:),'ApproveAll'));
    NumDatasets = length(prefixFilter);
else
    NumDatasets = SizeDataTab(2) - 1;
end

Prefixes = cell(1,NumDatasets);

if ~justReady
    for i = 2:(NumDatasets+1)
        PrefixCell = DataTab{PrefixRow,i};
        QuotePositions = strfind(PrefixCell,'''');
        PrefixName = PrefixCell((QuotePositions(1)+1):(QuotePositions(end)-1));
        Prefixes(i-1) = {PrefixName};
    end 
else
    for i = 1:length(prefixFilter)
        PrefixCell = DataTab{PrefixRow,prefixFilter(i)};
        QuotePositions = strfind(PrefixCell,'''');
        PrefixName = PrefixCell((QuotePositions(1)+1):(QuotePositions(end)-1));
        Prefixes(i) = {PrefixName};
    end 
end

% Extract the MetaData from each Prefix
PrefixFolders = {};
RawSettings = struct('Prefix','',  'SizeX',[], 'SizeY',[],... 
                    'SizeZ',[], 'PixelSizeX',[], 'PixelSizeY',[], ...
                    'PixelSizeZ',[], 'firstLaser', [], 'secondLaser', [], 'thirdLaser', [], 'pinholeSize', [],...
                    'gain1', [], 'gain2', [], 'zoom', []);
DatasetCounter = 1; %Counter only advances when dataset exists. Prevents empty entries in RawSettings, which can produce false negatives in ComparedSettings

for i = 1:NumDatasets
    % Need to regenerate the folder paths from the Prefixes
    CurrPrefix = char(Prefixes(i));
    DashPositions = strfind(CurrPrefix,'-');
    if numel(DashPositions) < 3
        error('Check that the Prefix entry for this dataset contains a dash (and not a slash) after the date.')
    end
    CurrPrefixFolder = [CurrPrefix(1:(DashPositions(3)-1)),filesep,...
            CurrPrefix((DashPositions(3)+1):end)];
    PrefixFolders(i,1) = {CurrPrefixFolder};
    CurrLSM = CurrPrefix((DashPositions(3)+1):end);
    
    % Check which type of microscopy data we have
    CurrRawDataPath = [RawDataPath, filesep, CurrPrefixFolder];
    try
        evalc('[~, FileMode] = DetermineFileMode(CurrRawDataPath)');    %Using evalc to supress displays to the command window from the function DetermineFileMode
    catch
        warning(['Either dataset does not exist or no tifs found for dataset ' num2str(i) '. Skipping settings extraction for this dataset.'])
        continue
    end
    
    if strcmpi(FileMode, 'LSM')
        CurrSettingStruct = [];

        % Read in and add to CurrSettingsStruct any settings available in
        % BioFormats
        % Read in only the metadata without having to open the .lif files
        evalc('MetaReader = bfGetReader([CurrRawDataPath, filesep, CurrLSM, ''_A.lsm''])');   %Using evalc to repress displays to the command window from the function bfGetReader
        MetaData = MetaReader.getMetadataStore(); 
        % Access the desired settings
        SeriesIndex = 0;
        CurrSettingStruct.Prefix = CurrPrefix;   % # pixels per frame

        CurrSettingStruct.SizeX = str2double(MetaData.getPixelsSizeX(SeriesIndex));   % # pixels per frame
        CurrSettingStruct.SizeY = str2double(MetaData.getPixelsSizeY(SeriesIndex));   % # pixels per frame
        CurrSettingStruct.SizeZ = str2double(MetaData.getPixelsSizeZ(SeriesIndex));   % # steps per stack
        CurrSettingStruct.PixelSizeX = str2double(MetaData.getPixelsPhysicalSizeX(SeriesIndex).value);   % um
        CurrSettingStruct.PixelSizeY = str2double(MetaData.getPixelsPhysicalSizeY(SeriesIndex).value);   % um
        CurrSettingStruct.PixelSizeZ = str2double(MetaData.getPixelsPhysicalSizeZ(SeriesIndex).value);   % um
        CurrSettingStruct.firstLaser = str2double(MetaData.getLaserWavelength(0, 0).value);
        try
            CurrSettingStruct.secondLaser = str2double(MetaData.getLaserWavelength(0, 1).value);
        catch
            CurrSettingStruct.secondLaser = [];
        end
        try
            CurrSettingStruct.thirdLaser = str2double(MetaData.getLaserWavelength(0, 2).value);
        catch
             CurrSettingStruct.thirdLaser = [];
        end
        CurrSettingStruct.pinholeSize = str2double(MetaData.getChannelPinholeSize(0,0).value);
        CurrSettingStruct.gain1 = str2double(MetaData.getDetectorGain(0,0));
        try
            CurrSettingStruct.gain2 = str2double(MetaData.getDetectorGain(0,1));
        catch
            CurrSettingStruct.gain2 = [];
        end
        CurrSettingStruct.zoom = str2double(MetaData.getDetectorZoom(0,0));

        % Close the reader when finished with it
        
        MetaReader.close;
        
        % Add current settings to the master RawSettings structure
        RawSettings(DatasetCounter) = CurrSettingStruct;
    else
        error('Only LSM datasets currently supported by this script.')
    end
    DatasetCounter = DatasetCounter + 1;
    disp(['Settings for dataset ' num2str(i) ' of ' num2str(NumDatasets) ' extracted.'])
end

%% Compare each setting across datasets

% Find the difference in settings between each pair of datasets, check for any
% nonzero entries indicating nonmatching settings, then save the logical 
% % result in the corresponding field
% ComparedSettings.pixelDwellTime = isempty(nonzeros(diff([RawSettings.pixelDwellTime])));
% ComparedSettings.lineAccumulation = isempty(nonzeros(diff([RawSettings.lineAccumulation])));
% ComparedSettings.frameAccumulation = isempty(nonzeros(diff([RawSettings.frameAccumulation])));
% ComparedSettings.frameAverage = isempty(nonzeros(diff([RawSettings.frameAverage])));
% ComparedSettings.emWaveForPinAiryCalc = isempty(nonzeros(diff([RawSettings.emWaveForPinAiryCalc])));
% ComparedSettings.pinholeAiry = isempty(nonzeros(diff(round([RawSettings.pinholeAiry],CompareDecimals))));
% ComparedSettings.pinhole = isempty(nonzeros(diff(round([RawSettings.pinhole],pinholeDecimals))));  %Need to round to prevent incorrect falses
% ComparedSettings.rotatorAngle = isempty(nonzeros(diff([RawSettings.rotatorAngle])));
% ComparedSettings.phaseX = isempty(nonzeros(diff([RawSettings.phaseX])));

% Compare the strings between settings, check for any zero entries 
% indicating nonmatching settings, then save the logical result in the 
% % corresponding field
% ComparedSettings.scanDirectionXName = ... 
%     isempty(nonzeros(~strcmp({RawSettings.scanDirectionXName}, ...
%                                 RawSettings(1).scanDirectionXName)));
% 
% ComparedSettings.scanDirectionX = isempty(nonzeros(diff([RawSettings.scanDirectionX])));
% ComparedSettings.zoomSetting = isempty(nonzeros(diff([RawSettings.zoomSetting])));
% ComparedSettings.scanSpeed = isempty(nonzeros(diff([RawSettings.scanSpeed])));
% ComparedSettings.refractionIndex = isempty(nonzeros(diff([RawSettings.refractionIndex])));
% ComparedSettings.numericalAperture = isempty(nonzeros(diff([RawSettings.numericalAperture])));
% 
% ComparedSettings.immersion = isempty(nonzeros(~strcmp({RawSettings.immersion}, ...
%                                                         RawSettings(1).immersion)));
% 
% ComparedSettings.objectiveNumber = isempty(nonzeros(diff([RawSettings.objectiveNumber])));
% ComparedSettings.magnification = isempty(nonzeros(diff([RawSettings.magnification])));
% ComparedSettings.frameTime = isempty(nonzeros(diff(round([RawSettings.frameTime],CompareDecimals))));
% % ComparedSettings.completeTime = isempty(nonzeros(diff([RawSettings.completeTime])));
% ComparedSettings.cycleTime = isempty(nonzeros(diff(round([RawSettings.cycleTime],CompareDecimals))));    %Need to round to prevent incorrect falses
% 
% ComparedSettings.zStackDirectionModeName = ... 
%     isempty(nonzeros(~strcmp({RawSettings.zStackDirectionModeName}, ...
%                                 RawSettings(1).zStackDirectionModeName)));
% ComparedSettings.zUseModeName = ... 
%     isempty(nonzeros(~strcmp({RawSettings.zUseModeName}, ...
%                                 RawSettings(1).zUseModeName)));
% ComparedSettings.xGalvoMovementModeName = ... 
%     isempty(nonzeros(~strcmp({RawSettings.xGalvoMovementModeName}, ...
%                                 RawSettings(1).xGalvoMovementModeName)));

ComparedSettings.SizeX = isempty(nonzeros(diff([RawSettings.SizeX])));
ComparedSettings.SizeY = isempty(nonzeros(diff([RawSettings.SizeY])));
ComparedSettings.SizeZ = isempty(nonzeros(diff([RawSettings.SizeZ])));
ComparedSettings.PixelSizeX = isempty(nonzeros(diff([RawSettings.PixelSizeX])));
ComparedSettings.PixelSizeY = isempty(nonzeros(diff([RawSettings.PixelSizeY])));
ComparedSettings.PixelSizeZ = isempty(nonzeros(diff(round([RawSettings.PixelSizeZ],CompareDecimals)))); %Need to round to prevent incorrect falses
ComparedSettings.firstLaser = isempty(nonzeros(diff([RawSettings.firstLaser])));
ComparedSettings.secondLaser = isempty(nonzeros(diff([RawSettings.secondLaser])));
ComparedSettings.thirdLaser = isempty(nonzeros(diff([RawSettings.thirdLaser])));
ComparedSettings.gain1 = isempty(nonzeros(diff([RawSettings.gain1])));
ComparedSettings.gain2 = isempty(nonzeros(diff([RawSettings.gain2])));
ComparedSettings.zoom = isempty(nonzeros(diff([RawSettings.zoom])));
ComparedSettings.pinholeSize = isempty(nonzeros(diff([RawSettings.pinholeSize])));






% Display comparison structure and final warnings
ComparedSettings
disp('Settings compared.')
disp('If any settings display a zero, check the RawSettings structure to identify the nonmatching dataset.')
disp('For some settings, nonmatching does not necessary mean a faulty dataset. Use your own judgement.')