% function [theStruct, settingStruct] = readSettingsFromXML(filename)
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

function [theStruct, settingStruct] = readSettingsFromXML(Prefix, filename)
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

    % Read in the .xml file
    try
       tree = xmlread(filename);
    catch
       error('Failed to read XML file %s.',filename);
    end

    % Recurse over child nodes. This could run into problems 
    % with very deeply nested trees.
    try
       theStruct = parseChildNodes(tree);
    catch
       error('Unable to parse XML file %s.',filename);
    end
    
    % Save settings into a structure
    settingStruct.Prefix = Prefix;
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

           for count = 1:numAttributes
              attrib = theAttributes.item(count-1);
              % Settings that only occur once or are always the same
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
              elseif strcmp(char(attrib.getName),'CompleteTime')
                  completeTime = attrib.getValue;
              elseif strcmp(char(attrib.getName),'CycleTime')
                  cycleTime = attrib.getValue;
              elseif strcmp(char(attrib.getName),'ZStackDirectionModeName')
                  zStackDirectionModeName = attrib.getValue;
              elseif strcmp(char(attrib.getName),'ZUseModeName')
                zUseModeName = attrib.getValue;
              elseif strcmp(char(attrib.getName),'XGalvoMovementModeName')
                xGalvoMovementModeName = attrib.getValue;
              end

              % Settings that occur multiple times with different values
        %       if strcmp(char(attrib.getName),'LaserLine')
        %           if exist('laserLine','var')
        %               newLaserLine = attrib.getValue;
        %               newLaserLines = {laserLine, newLaserLine}
        %               assignin('base', 'laserLine', newLaserLines);
        %           else
        %               assignin('base', 'laserLine', attrib.getValue);
        %           end
        %       elseif
        %       end

              attributes(count).Name = char(attrib.getName);
              attributes(count).Value = char(attrib.getValue);
           end
        end
    end
end