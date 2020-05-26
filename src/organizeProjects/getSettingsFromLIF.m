function settingsStruct = getSettingsFromLIF(liveExperiment)

% function settingsStruct = getSettingsFromLIF(liveExperiment)
%
% DESCRIPTION
% Compares the experimental settings (aka metadata) for multiple datasets
% and determines whether or not they are the identical. Only works for
% Leica datasets at the moment.
%
% PARAMETERS
% liveExperiment: LiveExperiment object for the Leica (LIF) experiment for 
%                 which you want to extract settings 
%
% OPTIONS
% N/A
%
% OUTPUT
% settingsStruct: Returns a structure containing all the microsocpe  
%                 settings for this experiment
%
%                 Includes these settings:
%                 experimentName, pixelDwellTime, lineAccumulation,
%                 frameAccumulation, lineAverage, frameAverage,
%                 emWaveForPinAiryCalc, pinholeAiry, pinhole,
%                 rotatorAngle, phaseX, scanDirectionXName, scanDirectionX,
%                 zoomSetting, scanSpeed, refractionIndex,
%                 numericalAperture, immersion, objectiveNumber, 
%                 magnification, frameTime, cycleTime, completeTime,
%                 zStackDirectionModeName, zUseModeName,
%                 xGalvoMovementModeName, laserLines, isLaserLineChecked,
%                 laserIntensityDev, laserAOBSIntensityDev, 
%                 laserAOBSIntensityLowDev, laserOutCheckedIntensity,
%                 detectorName, isDetectorActive, 
%                 detectorAcquisitionModeName, detectorAcquisitionMode,
%                 detectorTimeGateWavelength, detectorTimeGatePulseEnd,
%                 detectorTimeGatePulseStart, detectorIsTimeGateActivated,
%                 multibandTargetWaveLengthEnd, 
%                 multibandTargetWaveLengthBegin, multibandChannelName,
%                 multibandChannel, SizeX, SizeY, SizeZ, PixelSizeX, 
%                 PixelSizeY, PixelSizeZ
% 
%
% Author (contact): Meghan Turner (meghan_turner@berkeley.edu)
% Created: 2020-05-21
% Last Updated:
%
experimentName = liveExperiment.Prefix;
rawDataPath = liveExperiment.rawFolder;
rawDataFile = getRawDataFile(liveExperiment);   %used in an evalc statement in line 59

% settingsStruct.experimentName = liveExperiment.Prefix;

%% Read XML file to extract those settings not supported by Bio-Formats

metaDataDir = dir([rawDataPath, filesep, 'MetaData', filesep, '*.xml']); %the old LASX xml file(s)
xmlFile = [rawDataPath, filesep, 'lifMeta.xml'];   %the new, generated xml file

%If the dataset has a MetaData folder from LASX export, use that
if ~isempty(metaDataDir)
    disp('Using MetaData folder for settings extraction ...')
    xmlFileName = metaDataDir(1).name;   % By default, selects the first .xml file, probably not the best idea
    xmlFile = [rawDataPath, filesep, 'MetaData', filesep, xmlFileName];
    [~, settingsStruct] = readSettingsFromLIFMetaDataXML(experimentName, xmlFile);
    
%If a lifMeta.xml file already exists, use that
elseif exist(xmlFile, 'file')
    disp('Using lifMeta file for settings extraction...')
    [~, settingsStruct] = readSettingsFromLIFMetaDataXML(experimentName, xmlFile);
    
%Otherwise, generate and use a new lifMeta.xml file
else
    disp('Generating new lifMeta file for settings extraction...')
    generateLIFMetaDataXML(currExperiment, xmlFile);
    [~, settingsStruct] = readSettingsFromLIFMetaDataXML(experimentName, xmlFile);
end


%% Read in and add to settingsStruct all settings available in BioFormats

% Read in only the metadata without having to open the .lif files
evalc('MetaReader = bfGetReader(rawDataFile)');   %Using evalc to repress displays to the command window from the function bfGetReader
% MetaReader = bfGetReader(rawDataFile);
MetaData = MetaReader.getMetadataStore();

%Grab the desired settings
SeriesIndex = 0;
settingsStruct.SizeX = str2double(MetaData.getPixelsSizeX(SeriesIndex));   % # pixels per frame
settingsStruct.SizeY = str2double(MetaData.getPixelsSizeY(SeriesIndex));   % # pixels per frame
settingsStruct.SizeZ = str2double(MetaData.getPixelsSizeZ(SeriesIndex));   % # steps per stack
settingsStruct.PixelSizeX = str2double(MetaData.getPixelsPhysicalSizeX(SeriesIndex).value);   % um
settingsStruct.PixelSizeY = str2double(MetaData.getPixelsPhysicalSizeY(SeriesIndex).value);   % um
settingsStruct.PixelSizeZ = str2double(MetaData.getPixelsPhysicalSizeZ(SeriesIndex).value);   % um

%Close the reader when finished with it
MetaReader.close;
