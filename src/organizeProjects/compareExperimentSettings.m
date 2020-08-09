function [comparedSettings,rawSettings] = compareExperimentSettings(dataTypes,varargin)

% function [ComparedSettings,RawSettings] = compareExperimentSettings(dataTypes,varargin)
%
% DESCRIPTION
% Compares the experimental settings (aka metadata) for multiple datasets
% and determines whether or not they are the identical. Only works for
% Leica datasets at the moment.
%
% PARAMETERS
% dataTypes: This is a cell array of char variable(s) which are the exact 
%           name(s) of the tab(s) in DataStatus.xlsx that you wish to 
%           analyze. Single dataType can also be passed as a string (char
%           array).
%           E.g. dataTypes = {'dataType1'} OR dataType = 'dataType1' for 
%                only one tab
%           E.g. dataTypes = {'dataType1', 'dataType2', 'dataType3'} for 
%                multiple tabs
% 
%
% OPTIONS
% 'onlyApproved': Only compares the settings of Prefixes that have "Ready" or
%                 "ApproveAll" in the Compile
% 'defaultTolerance': Should be followed by an integer between 1 and 4,
%                     indicating the number of decimal places you want
%                     compared to determine if the settings match
% 'laserTolerance': Should be followed by a number between 0 and 1
%                   indicating the percent tolerance that is allowed when 
%                   comparing the laser power between datasets. The default
%                   tolerance is 0.10.
%
%
% OUTPUT
% comparedSettings: Returns a structure containing booleans indicating 
%                   which settings match.                 
% rawSettings: Returns a structure containing all the microsocpe settings 
%           for the datasets found in the indicated DataStatus tab,   
%           allowing for easy visual comparison. 
%           See getSettingsFromLIF.m for which settings are currently 
%           included.
% comparedSettings_[...].mat: MAT file containing comparedSettings, 
%                             rawSettings, and dataType. Saved to the same
%                             results folder where the DataStatus.xlsx file
%                             containing your dataTypes exists
% 
%
% Author (contact): Meghan Turner (meghan_turner@berkeley.edu)
% Created: 2019-04-30
% Last Updated: 2020-05-21 MT
%

%Initialize variables
onlyApproved = false; 

%Default tolerances to consider settings the same - these defaults were 
% selected by MT after looking at typical variations in LIF settings

%Decimal place tolerances
defaultDecimals = 2;
pixelXYDecimals = 3;
pinholeDecimals = 4;
%Percent tolerances
laserPercent = 0.1;


%Determine user inputs
for i= 1:length(varargin)
    if strcmpi(varargin{i},'defaultTolerance')
        defaultDecimals = varargin{i+1};
    elseif strcmpi(varargin{i},'laserTolerance')
        laserTolerance = varargin{i+1};
        if laserTolerance < 0 || laserTolerance > 1
            error('laserTolerance must be a percentage, expressed as a number between 0 and 1.')
        end
    elseif strcmpi(varargin{i}, 'onlyApproved')
        onlyApproved = true;
    end
end

%Store tolerances in a table to easily pass to other functions
settingsTolerances = table(defaultDecimals,pinholeDecimals,pixelXYDecimals,laserPercent);

%Convert to cell array if the user passed a single dataType as a char array
if ischar(dataTypes)
    projectList = {dataTypes};
else
    projectList = dataTypes;
end

% Get the experiment names (formerly prefixes) for the project(s)
experimentNames = {};
for i = 1:length(projectList)
    projectName = projectList{i};
    
    %Grab only the experiments the user wants
    if onlyApproved    
        newExperimentNames = getProjectPrefixes(projectName,'onlyApproved');
    else
        newExperimentNames = getProjectPrefixes(projectName);
    end
    
    %Add new experiments into the master list
    numNew = length(newExperimentNames);
    experimentNames(end+1:end+numNew,1) = newExperimentNames;
end


%% Extract the MetaData from each experiment in the project

numExperiments = length(experimentNames);
experimentExistsCount = 1; %Counter only advances when dataset exists. Prevents empty entries in RawSettings, which can produce false negatives in ComparedSettings

for i = 1:numExperiments
    %Get raw data folder and file info for this experiment
    currExperimentName = experimentNames{i};
    currExperiment = LiveExperiment(currExperimentName);
    
    %Only proceed with settings extraction if we're dealing with a LIF file
    %(See compareExperimentSettingsZeiss for Zeiss support)
    FileMode = currExperiment.fileMode; 
    if isempty(FileMode)
        warning(['Dataset ' currExperimentName ' not detected. Skipping settings extraction for this dataset.'])
        continue    
    elseif strcmpi(FileMode, 'LIFEXport')
        %Extract settings for this experiment
        try
            currSettingsStruct = getSettingsFromLIF(currExperiment);

            % Add current settings to the master RawSettings structure
            rawSettings(experimentExistsCount) = currSettingsStruct;  
        catch
            warning(['Failed to extract settings for experiment ', i]);
        end
    else
        error('Only Leica datasets currently supported by this script.')
    end
    
    experimentExistsCount = experimentExistsCount + 1;
    disp(['Settings for dataset ' num2str(i) ' of ' num2str(numExperiments) ' extracted.'])
end

%% Compare each setting across datasets
comparedSettings = compareLIFSettings(rawSettings,settingsTolerances)  %Not supressed so it will display in command window

disp('Settings compared and saved in your specified DynamicResults folder. ')
disp('If any settings display a zero, check the RawSettings structure to identify the nonmatching dataset.')
disp('For some settings, nonmatching does not necessary mean a faulty dataset. Use your own judgement.')


%% Save structures for later reference
%Use the last experiment's specific results folder to find the root results
% folder
resultsFolder = currExperiment.resultsFolder;
filesepPositions = strfind(resultsFolder,filesep);
resultsFolder = resultsFolder(1:filesepPositions(end-1));

%Name is the first project from projectList plus an indicator of whether 
%other projects were also compared
if numel(projectList) == 1
    resultsFolder = [resultsFolder, projectList{1}, filesep];
    if ~exist(resultsFolder,'dir')
        mkdir(resultsFolder)
    end
    saveName = ['comparedSettings_', projectList{1}];
else
    saveName = ['comparedSettings_', projectList{1}, '_plusOtherProjects'];
end

save([resultsFolder,saveName,'.mat'],'rawSettings','comparedSettings','dataTypes');