function [theStruct, settingStruct] = readSettingsFromLIFMetaDataXML(Prefix, filename)

% function [theStruct, settingStruct] = readSettingsFromLIFMetaDataXML(Prefix, filename)
%
% DESCRIPTION
% Searches Leica XML files for the following microscope settings and  
% saves them as fields in a structure to be accessed by other scripts. This
% allows for direct access to settings that are not supported by
% Bio-Formats.
% 
% ARGUMENTS
% filename: Full or relative path to the .xml file to be read
% 
% OUTPUT
% theStruct: Structure created by recursing over the nodes of the .xml
%            document
% settingStruct: Structure containing fields corresponding to the
%                settings read from teh .xml document. Numbers have
%                already been converted to doubles.
%
% CALLED BY
% compareExperimentSettings.m
% 
% Author (contact): Meghan Turner (meghan_turner@berkeley.edu)
% Created: 2019-04-30
% Last Updated: 2019-05-07
%
% Based on the searchXML.m script written by Armando Reimer. Uses nested
%   functions and structures instead of local functions & assignin/evalin.
%
% Documented by: Meghan Turner (meghan_turner@berkeley.edu)


    % Initilaize variables so they're accessible in the nested functions
    pixelDwellTime = [];
    lineAccumulation = [];
    frameAccumulation = [];
    lineAverage = [];
    frameAverage = [];
    emWaveForPinAiryCalc = [];
    pinholeAiry = [];
    pinhole = [];
    rotatorAngle = [];
    phaseX = [];
    scanDirectionXName = '';
    scanDirectionX = [];
    zoomSetting = [];
    scanSpeed = [];
    refractionIndex = [];
    numericalAperture = [];
    immersion = '';
    objectiveNumber = [];
    magnification = [];
    frameTime = [];
%     completeTime = [];
    cycleTime = [];
    zStackDirectionModeName = '';
    zUseModeName = '';
    xGalvoMovementModeName = '';
    laserLines = [];
    isLaserLineChecked = [];
    laserIntensityDev = [];
    laserAOBSIntensityDev = [];
    laserAOBSIntensityLowDev = [];
    laserOutCheckedIntensity = [];
    detectorName = {};
    isDetectorActive = [];
    detectorAcquisitionModeName = {};
    detectorAcquisitionMode = [];
    detectorTimeGateWavelength = [];
    detectorTimeGatePulseEnd = [];
    detectorTimeGatePulseStart = [];
    detectorIsTimeGateActivated = [];
    multibandTargetWaveLengthEnd = [];
    multibandTargetWaveLengthBegin = [];
    multibandChannelName = {};
    multibandChannel = [];

    % Read in the .xml file
    try
       tree = xmlread(filename);
    catch
       error('Failed to read XML file %s.',filename);
    end

    % Recurse over child nodes. This could run into problems 
    % with very deeply nested trees.
%     try
       theStruct = parseChildNodes(tree);
%     catch
%        error('Unable to parse XML file %s.',filename);
%     end
    
    % Save settings into a structure
    % hardware settings
    settingStruct.experimentName = Prefix;
    settingStruct.pixelDwellTime = str2double(pixelDwellTime);
    settingStruct.lineAccumulation = str2double(lineAccumulation);
    settingStruct.frameAccumulation = str2double(frameAccumulation);
    settingStruct.lineAverage = str2double(lineAverage);
    settingStruct.frameAverage = str2double(frameAverage);
    settingStruct.emWaveForPinAiryCalc = str2double(emWaveForPinAiryCalc);  % in nm; emission wavelength for the pinhole airy calculation
    settingStruct.pinholeAiry = str2double(pinholeAiry);    % in AU, should be around 1 AU
    settingStruct.pinhole = str2double(pinhole);    % in um
    settingStruct.rotatorAngle = str2double(rotatorAngle);  % in degrees
    settingStruct.phaseX = str2double(phaseX);
    settingStruct.scanDirectionXName = char(scanDirectionXName);
    settingStruct.scanDirectionX = str2double(scanDirectionX);
    settingStruct.zoomSetting = str2double(zoomSetting);
    settingStruct.scanSpeed = str2double(scanSpeed);
    settingStruct.refractionIndex = str2double(refractionIndex);
    settingStruct.numericalAperture = str2double(numericalAperture);
    settingStruct.immersion = char(immersion);
    settingStruct.objectiveNumber = str2double(objectiveNumber);
    settingStruct.magnification = str2double(magnification);
    settingStruct.frameTime = str2double(frameTime);
%     settingStruct.completeTime = str2double(completeTime);
    settingStruct.cycleTime = str2double(cycleTime);
    settingStruct.zStackDirectionModeName = char(zStackDirectionModeName);
    settingStruct.zUseModeName = char(zUseModeName);
    settingStruct.xGalvoMovementModeName = char(xGalvoMovementModeName);
    % laser line settings
    settingStruct.laserLines = laserLines;
    settingStruct.isLaserLineChecked = isLaserLineChecked;
    settingStruct.laserIntensityDev = laserIntensityDev; 
    settingStruct.laserAOBSIntensityDev = laserAOBSIntensityDev;
    settingStruct.laserAOBSIntensityLowDev = laserAOBSIntensityLowDev;
    settingStruct.laserOutCheckedIntensity = laserOutCheckedIntensity;
    % detector settings
    settingStruct.detectorName = detectorName;
    settingStruct.isDetectorActive = isDetectorActive;
    settingStruct.detectorAcquisitionModeName = detectorAcquisitionModeName;
    settingStruct.detectorAcquisitionMode = detectorAcquisitionMode;
    settingStruct.detectorTimeGateWavelength = detectorTimeGateWavelength;
    settingStruct.detectorTimeGatePulseEnd = detectorTimeGatePulseEnd;
    settingStruct.detectorTimeGatePulseStart = detectorTimeGatePulseStart;
    settingStruct.detectorIsTimeGateActivated = detectorIsTimeGateActivated;
    settingStruct.multibandTargetWaveLengthEnd = multibandTargetWaveLengthEnd;
    settingStruct.multibandTargetWaveLengthBegin  = multibandTargetWaveLengthBegin;
    settingStruct.multibandChannelName = multibandChannelName;
    settingStruct.multibandChannel = multibandChannel;


    % ----- Nested function PARSECHILDNODES -----
    function children = parseChildNodes(theNode)
        % Recurse over node children.
        children = [];
        if theNode.hasChildNodes
           childNodes = theNode.getChildNodes;
           numChildNodes = childNodes.getLength;
           allocCell = cell(1, numChildNodes);

           children = struct(             ...
              'Name', allocCell, 'Attributes', allocCell,    ...
              'Data', allocCell, 'Children', allocCell);

            for count = 1:numChildNodes
                theChild = childNodes.item(count-1);
                children(count) = makeStructFromNode(theChild);
            end
        end
    end

    % ----- Nested function MAKESTRUCTFROMNODE -----
    function nodeStruct = makeStructFromNode(theNode)
        % Create structure of node info.

        nodeStruct = struct(                        ...
           'Name', char(theNode.getNodeName),       ...
           'Attributes', parseAttributes(theNode),  ...
           'Data', '',                              ...
           'Children', parseChildNodes(theNode));

        if any(strcmp(methods(theNode), 'getData'))
           nodeStruct.Data = char(theNode.getData); 
        else
           nodeStruct.Data = '';
        end
    end

    % ----- Nested function PARSEATTRIBUTES -----
    function attributes = parseAttributes(theNode)
        % Create attributes structure.

        attributes = [];
        if theNode.hasAttributes
           theAttributes = theNode.getAttributes;
           numAttributes = theAttributes.getLength;
           allocCell = cell(1, numAttributes);
           attributes = struct('Name', allocCell, 'Value', ...
                               allocCell);
           nodeName = char(theNode.getNodeName);
           
           % Create attributes structure for all nodes and extract
           % attribute values for the settings we care about
           
           % extract most hardware settings
           if strcmp(nodeName,'ATLConfocalSettingDefinition')
               for count = 1:numAttributes
                  attrib = theAttributes.item(count-1);
                  attributes(count).Name = char(attrib.getName);
                  attributes(count).Value = char(attrib.getValue);
                  
                  if strcmp(char(attrib.getName),'PixelDwellTime')
                      pixelDwellTime = attrib.getValue;
                  elseif strcmp(char(attrib.getName),'Line_Accumulation')
                      lineAccumulation = attrib.getValue;
                  elseif strcmp(char(attrib.getName),'FrameAccumulation')
                      assignin('base', 'frameAccumulation', attrib.getValue);
                      frameAccumulation = attrib.getValue;
                  elseif strcmp(char(attrib.getName),'LineAverage')
                      lineAverage = attrib.getValue;
                  elseif strcmp(char(attrib.getName),'FrameAverage')
                      frameAverage = attrib.getValue;
                  elseif strcmp(char(attrib.getName),'EmissionWavelengthForPinholeAiryCalculation')
                      emWaveForPinAiryCalc = attrib.getValue;
                  elseif strcmp(char(attrib.getName),'PinholeAiry')
                      pinholeAiry = attrib.getValue;
                  elseif strcmp(char(attrib.getName),'Pinhole')
                      pinhole = attrib.getValue;
                  elseif strcmp(char(attrib.getName),'RotatorAngle')
                      rotatorAngle = attrib.getValue;
                  elseif strcmp(char(attrib.getName),'PhaseX')
                      phaseX = attrib.getValue;
                  elseif strcmp(char(attrib.getName),'ScanDirectionXName')
                      scanDirectionXName = attrib.getValue;
                  elseif strcmp(char(attrib.getName),'ScanDirectionX')
                      scanDirectionX = attrib.getValue;
                  elseif strcmp(char(attrib.getName),'Zoom')
                      zoomSetting = attrib.getValue;
                  elseif strcmp(char(attrib.getName),'ScanSpeed')
                      scanSpeed = attrib.getValue;
                  elseif strcmp(char(attrib.getName),'RefractionIndex')
                      refractionIndex = attrib.getValue;
                  elseif strcmp(char(attrib.getName),'NumericalAperture')
                      numericalAperture = attrib.getValue;
                  elseif strcmp(char(attrib.getName),'Immersion')
                      immersion = attrib.getValue;
                  elseif strcmp(char(attrib.getName),'ObjectiveNumber')
                      objectiveNumber = attrib.getValue;
                  elseif strcmp(char(attrib.getName),'Magnification')
                      magnification = attrib.getValue;
                  elseif strcmp(char(attrib.getName),'FrameTime')
                      frameTime = attrib.getValue;
%                   elseif strcmp(char(attrib.getName),'CompleteTime')
%                       completeTime = attrib.getValue;
                  elseif strcmp(char(attrib.getName),'CycleTime')
                      cycleTime = attrib.getValue;
                  elseif strcmp(char(attrib.getName),'ZStackDirectionModeName')
                      zStackDirectionModeName = attrib.getValue;
                  elseif strcmp(char(attrib.getName),'ZUseModeName')
                    zUseModeName = attrib.getValue;
                  elseif strcmp(char(attrib.getName),'XGalvoMovementModeName')
                    xGalvoMovementModeName = attrib.getValue;
                  end
               end
           % Get information about the lasers
           elseif strcmp(nodeName, 'LaserLineSetting')
               for count = 1:numAttributes
                  attrib = theAttributes.item(count-1);
                  attributes(count).Name = char(attrib.getName);
                  attributes(count).Value = char(attrib.getValue);
                  
                  if strcmp(char(attrib.getName),'LaserLine')
                      laserLines(end+1) = str2double(attrib.getValue);
                  elseif strcmp(char(attrib.getName),'IsLineChecked')
                      isLaserLineChecked(end+1) = str2double(attrib.getValue);
                  elseif strcmp(char(attrib.getName),'IntensityDev')
                      laserIntensityDev(end+1) = str2double(attrib.getValue);
                  elseif strcmp(char(attrib.getName),'AOBSIntensityDev')
                      laserAOBSIntensityDev(end+1) = str2double(attrib.getValue);
                  elseif strcmp(char(attrib.getName),'AOBSIntensityLowDev')
                      laserAOBSIntensityLowDev(end+1) = str2double(attrib.getValue);   
                  elseif strcmp(char(attrib.getName),'OutCheckedIntensity')
                      laserOutCheckedIntensity(end+1) = str2double(attrib.getValue);
                  end
               end
           % extract detector settings
           elseif strcmp(nodeName, 'Detector')
               for count = 1:numAttributes
                  attrib = theAttributes.item(count-1);
                  attributes(count).Name = char(attrib.getName);
                  attributes(count).Value = char(attrib.getValue);
                  
                  if strcmp(char(attrib.getName),'Name')
                      detectorName{end+1} = char(attrib.getValue);
                  elseif strcmp(char(attrib.getName),'IsActive')
                      isDetectorActive(end+1) = str2double(attrib.getValue);   
                  elseif strcmp(char(attrib.getName),'AcquisitionModeName')
                      detectorAcquisitionModeName{end+1} = char(attrib.getValue); 
                  elseif strcmp(char(attrib.getName),'AcquisitionMode')
                      detectorAcquisitionMode(end+1) = str2double(attrib.getValue); 
                  elseif strcmp(char(attrib.getName),'TimeGateWavelength')
                      detectorTimeGateWavelength(end+1) = str2double(attrib.getValue); 
                  elseif strcmp(char(attrib.getName),'TimeGatePulseEnd')
                      detectorTimeGatePulseEnd(end+1) = str2double(attrib.getValue); 
                  elseif strcmp(char(attrib.getName),'TimeGatePulseStart')
                      detectorTimeGatePulseStart(end+1) = str2double(attrib.getValue);    
                  elseif strcmp(char(attrib.getName),'IsTimeGateActivated')
                      detectorIsTimeGateActivated(end+1) = str2double(attrib.getValue);
                  end
               end
           % extract detector ranges, which are stored in the MultiBand
           % element
           elseif strcmp(nodeName, 'MultiBand')
               for count = 1:numAttributes
                  attrib = theAttributes.item(count-1);
                  attributes(count).Name = char(attrib.getName);
                  attributes(count).Value = char(attrib.getValue);
                  
                  if strcmp(char(attrib.getName),'RightWorld')
                      multibandTargetWaveLengthEnd(end+1) = str2double(attrib.getValue);
                  elseif strcmp(char(attrib.getName),'LeftWorld')
                      multibandTargetWaveLengthBegin(end+1) = str2double(attrib.getValue);
                  elseif strcmp(char(attrib.getName),'ChannelName')
                      multibandChannelName{end+1} = char(attrib.getValue);
                  elseif strcmp(char(attrib.getName),'Channel')
                      multibandChannel(end+1) = str2double(attrib.getValue);
                  end
               end
           % don't do any extraction for any other elements    
           else
               for count = 1:numAttributes
                  attrib = theAttributes.item(count-1);
                  attributes(count).Name = char(attrib.getName);
                  attributes(count).Value = char(attrib.getValue);
               end
           end
        end
    end
end