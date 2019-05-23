% function [ComparedSettings,RawSettings] = compareExperimentSettings(DynamicsResultsPath,RawDataPath,DataType,varargin)
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
% DataType: This is a string that is identical to the name of the tab in
%           DataStatus.xlsx that you wish to analyze.
% 
%
% OPTIONS
% 'LaserTolerance': Should be followed by a number between 0 and 1
%                   indicating the percent tolerance that is allowed when 
%                   comparing the laser power between datasets. The default
%                   tolerance is 0.10.
%
% OUTPUT
% ComparedSettings: Returns a ?? containing booleans indicating whether all
%                   settings match.
% RawSettings: Returns a table containing all the microsocpe settings for
%           the datasets found in the indicated DataStatus tab, allowing  
%           for easy visual comparison.
% 
%
% Author (contact): Meghan Turner (meghan_turner@berkeley.edu)
% Created: 2019-04-30
% Last Updated: 2019-05-07
%
% Documented by: Meghan Turner (meghan_turner@berkeley.edu)

function [ComparedSettings,RawSettings] = compareExperimentSettings(DynamicsResultsPath,RawDataPath,DataType,varargin)
 
LaserTolerance = 0.1;
pinholeTolerance = 4;   % Number of decimals to which to round off the pinhole to consider settings the same
cycleTimeTolerance = 2;   % Number of decimals to which to round off the cycleTime to consider settings the same
pixelSizeZTolerance = 2;    % Number of decimals to which to round off the pixelSizeZ to consider settings the same

for i= 1:length(varargin)
    if strcmpi(varargin{i},'LaserTolerance')
        LaserTolerance = varargin{i+1};
        if LaserTolerance < 0 || LaserTolerance > 1
            error('LaserTolerance must be between 0 and 1.')
        end
    end
end

% Look in DataStatus.xlsx in folderPath and find the tab given by
% the input variable dataType.
if isfile([DynamicsResultsPath,'\','DataStatus.xlsx']) || ...
                        isfile([DynamicsResultsPath,'\','DataStatus.csv'])  % NB: isfile doesn't let you use * to check for multiple file extensions
    DataStatusDir = dir([DynamicsResultsPath,filesep,'DataStatus.*']);
    [~,DataTab]=xlsread([DynamicsResultsPath,filesep,DataStatusDir(1).name],DataType);
else
    error(['DataStatus.xlsx does not exist in the indicated folder path, ' ...
            DynamicsResultsPath, '.'])
end

%Find and load the different prefixes
PrefixRow = find(strcmpi(DataTab(:,1),'Prefix:'));
SizeDataTab = size(DataTab);
NumDatasets = SizeDataTab(2) - 1;
Prefixes = cell(1,NumDatasets);
for i = 2:(NumDatasets+1)
    PrefixCell = DataTab{PrefixRow,i};
    QuotePositions = strfind(PrefixCell,'''');
    PrefixName = PrefixCell((QuotePositions(1)+1):(QuotePositions(end)-1));
    Prefixes(i-1) = {PrefixName};
end 

% Extract the MetaData from each Prefix
PrefixFolders = {};
RawSettings = struct('Prefix','', 'pixelDwellTime',[], ...                  % Initialize the struct
                    'lineAccumulation',[], 'frameAccumulation',[],...
                    'lineAverage',[], 'frameAverage',[], ... 
                    'emWaveForPinAiryCalc',[], 'pinholeAiry', [], ...
                    'pinhole',[], 'rotatorAngle',[], 'phaseX',[],...
                    'scanDirectionXName','',  'scanDirectionX',[], ...
                    'zoomSetting',[], 'scanSpeed',[], 'refractionIndex',[],...
                    'numericalAperture',[], 'immersion','', ... 
                    'objectiveNumber',[], 'magnification',[], 'frameTime',[],...
                    'cycleTime',[], ...%'completeTime',[],
                    'zStackDirectionModeName','', 'zUseModeName', '',...
                    'xGalvoMovementModeName','', 'SizeX',[], 'SizeY',[],... 
                    'SizeZ',[], 'PixelSizeX',[], 'PixelSizeY',[], ...
                    'PixelSizeZ',[]);
for i = 1:NumDatasets
    % Need to regenerate the folder paths from the Prefixes
    CurrPrefix = char(Prefixes(i));
    DashPositions = strfind(CurrPrefix,'-');
    CurrPrefixFolder = [CurrPrefix(1:(DashPositions(3)-1)),filesep,...
            CurrPrefix((DashPositions(3)+1):end)];
    PrefixFolders(i,1) = {CurrPrefixFolder};
    CurrLIF = CurrPrefix((DashPositions(3)+1):end);
    
    % Check which type of microscopy data we have
    % MT: only supports LIF files at the moment
    CurrRawDataPath = [RawDataPath, filesep, CurrPrefixFolder];
    [~, FileMode] = DetermineFileMode(CurrRawDataPath);
    
    if strcmpi(FileMode, 'LIFEXport')
        CurrSettingStruct = [];
        % Read in the .xml file for those settings not supported by
        % Bio-Formats
        if isfolder([CurrRawDataPath, filesep, 'MetaData'])
            tic
            xml_file_path = dir([CurrRawDataPath, filesep, 'MetaData', filesep, '*.xml']);
            xml_file = xml_file_path(1).name;   % By default, selects the first .xml file, probably not the best idea
            [~, CurrSettingStruct] = readSettingsFromXML(CurrPrefix, [CurrRawDataPath,...
                                             filesep, 'MetaData', filesep, xml_file]);
        else 
            warning('No MetaData folder found. Have you exportedDataForLivemRNA yet?')
        end
        
        % Read in and add to CurrSettingsStruct any settings available in
        % BioFormats
        % Read in only the metadata without having to open the .lif files
        MetaReader = bfGetReader([CurrRawDataPath, filesep, CurrLIF, '.lif']);
        MetaData = MetaReader.getMetadataStore(); 
        % Access the desired settings
        SeriesIndex = 0;
        CurrSettingStruct.SizeX = str2double(MetaData.getPixelsSizeX(SeriesIndex));   % # pixels per frame
        CurrSettingStruct.SizeY = str2double(MetaData.getPixelsSizeY(SeriesIndex));   % # pixels per frame
        CurrSettingStruct.SizeZ = str2double(MetaData.getPixelsSizeZ(SeriesIndex));   % # steps per stack
        CurrSettingStruct.PixelSizeX = str2double(MetaData.getPixelsPhysicalSizeX(SeriesIndex).value);   % um
        CurrSettingStruct.PixelSizeY = str2double(MetaData.getPixelsPhysicalSizeY(SeriesIndex).value);   % um
        CurrSettingStruct.PixelSizeZ = str2double(MetaData.getPixelsPhysicalSizeZ(SeriesIndex).value);   % um
        % Close the reader when finished with it
        MetaReader.close;
        
        % Add current settings to the master RawSettings structure
        RawSettings(i) = CurrSettingStruct;
        
    else
        error('Only Leica datasets currently supported by this script.')
    end
    disp(['Settings for dataset ' num2str(i) ' of ' num2str(NumDatasets) ' extracted.'])
end

%% Compare each setting across datasets

% Find the difference in settings between each pair of datasets, check for any
% nonzero entries indicating nonmatching settings, then save the logical 
% result in the corresponding field
ComparedSettings.pixelDwellTime = isempty(nonzeros(diff([RawSettings.pixelDwellTime])));
ComparedSettings.lineAccumulation = isempty(nonzeros(diff([RawSettings.lineAccumulation])));
ComparedSettings.frameAccumulation = isempty(nonzeros(diff([RawSettings.frameAccumulation])));
ComparedSettings.frameAverage = isempty(nonzeros(diff([RawSettings.frameAverage])));
ComparedSettings.emWaveForPinAiryCalc = isempty(nonzeros(diff([RawSettings.emWaveForPinAiryCalc])));
ComparedSettings.pinholeAiry = isempty(nonzeros(diff([RawSettings.pinholeAiry])));
ComparedSettings.pinhole = isempty(nonzeros(diff(round([RawSettings.pinhole],pinholeTolerance))));  %Need to round to prevent incorrect falses
ComparedSettings.rotatorAngle = isempty(nonzeros(diff([RawSettings.rotatorAngle])));
ComparedSettings.phaseX = isempty(nonzeros(diff([RawSettings.phaseX])));

% Compare the strings between settings, check for any zero entries 
% indicating nonmatching settings, then save the logical result in the 
% corresponding field
ComparedSettings.scanDirectionXName = ... 
    isempty(nonzeros(~strcmp({RawSettings.scanDirectionXName}, ...
                                RawSettings(1).scanDirectionXName)));

ComparedSettings.scanDirectionX = isempty(nonzeros(diff([RawSettings.scanDirectionX])));
ComparedSettings.zoomSetting = isempty(nonzeros(diff([RawSettings.zoomSetting])));
ComparedSettings.scanSpeed = isempty(nonzeros(diff([RawSettings.scanSpeed])));
ComparedSettings.refractionIndex = isempty(nonzeros(diff([RawSettings.refractionIndex])));
ComparedSettings.numericalAperture = isempty(nonzeros(diff([RawSettings.numericalAperture])));

ComparedSettings.immersion = isempty(nonzeros(~strcmp({RawSettings.immersion}, ...
                                                        RawSettings(1).immersion)));

ComparedSettings.objectiveNumber = isempty(nonzeros(diff([RawSettings.objectiveNumber])));
ComparedSettings.magnification = isempty(nonzeros(diff([RawSettings.magnification])));
ComparedSettings.frameTime = isempty(nonzeros(diff([RawSettings.frameTime])));
% ComparedSettings.completeTime = isempty(nonzeros(diff([RawSettings.completeTime])));
ComparedSettings.cycleTime = isempty(nonzeros(diff(round([RawSettings.cycleTime],cycleTimeTolerance))));    %Need to round to prevent incorrect falses

ComparedSettings.zStackDirectionModeName = ... 
    isempty(nonzeros(~strcmp({RawSettings.zStackDirectionModeName}, ...
                                RawSettings(1).zStackDirectionModeName)));
ComparedSettings.zUseModeName = ... 
    isempty(nonzeros(~strcmp({RawSettings.zUseModeName}, ...
                                RawSettings(1).zUseModeName)));
ComparedSettings.xGalvoMovementModeName = ... 
    isempty(nonzeros(~strcmp({RawSettings.xGalvoMovementModeName}, ...
                                RawSettings(1).xGalvoMovementModeName)));

ComparedSettings.SizeX = isempty(nonzeros(diff([RawSettings.SizeX])));
ComparedSettings.SizeY = isempty(nonzeros(diff([RawSettings.SizeY])));
ComparedSettings.SizeZ = isempty(nonzeros(diff([RawSettings.SizeZ])));
ComparedSettings.PixelSizeX = isempty(nonzeros(diff([RawSettings.PixelSizeX])));
ComparedSettings.PixelSizeY = isempty(nonzeros(diff([RawSettings.PixelSizeY])));
ComparedSettings.PixelSizeZ = isempty(nonzeros(diff(round([RawSettings.PixelSizeZ],pixelSizeZTolerance)))); %Need to round to prevent incorrect falses

% Print of comparison and final warnings
ComparedSettings
disp('Settings compared.')
disp('If any settings display a zero, check the RawSettings structure to identify the nonmatching dataset.')
disp('For some settings, nonmatching does not necessary mean a faulty dataset. Use your own judgement.')