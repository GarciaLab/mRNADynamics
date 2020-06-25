function comparedSettings = compareLIFSettings(rawSettings,settingsTolerances)

% function settingsStruct = getSettingsFromLIF(liveExperiment)
%
% DESCRIPTION
% Compares the experimental settings (aka metadata) for multiple LIF 
% (Leica) datasets and determines whether or not they are the identical.
%
% PARAMETERS
% rawSettings: Structure containing all microscope settings for multiple
%              experiments that you would like to compare
% settingsTolerances: Table containing tolerances for various settings, 
%                     e.g. number of decimal places that need to match for 
%                     two settings to be considered the same
%
% OPTIONS
% N/A
%
% OUTPUT
% comparedSettings: Returns a structure of booleans that indicate whether
%                   or not all each setting matches across all experiments
%
%                   Compares only the following settings:
%                   experimentName, pixelDwellTime, lineAccumulation,
%                   frameAccumulation, lineAverage, frameAverage,
%                   emWaveForPinAiryCalc, pinholeAiry, pinhole,
%                   rotatorAngle, phaseX, scanDirectionXName, scanDirectionX,
%                   zoomSetting, scanSpeed, refractionIndex,
%                   numericalAperture, immersion, objectiveNumber, 
%                   magnification, frameTime, cycleTime, completeTime,
%                   zStackDirectionModeName, zUseModeName,
%                   xGalvoMovementModeName, SizeX, SizeY, SizeZ, PixelSizeX, 
%                   PixelSizeY, PixelSizeZ
% 
%
% Author (contact): Meghan Turner (meghan_turner@berkeley.edu)
% Created: 2020-05-21
% Last Updated:
%

%Initializing structure so these are displayed in a particular order
comparedSettings = struct('pixelDwellTime',0, 'lineAccumulation',0, ...
                          'frameAccumulation',0, 'frameAverage',0,...
                          'emWaveForPinAiryCalc',0, 'pinholeAiry', 0,...
                          'pinhole',0, 'rotatorAngle',0, 'phaseX',0, ...
                          'scanDirectionXName',0, 'scanDirectionX',0, ...
                          'zoomSetting',0, 'scanSpeed',0, ...
                          'refractionIndex',0,'numericalAperture',0, ...
                          'immersion',0, 'objectiveNumber',0, ...
                          'magnification',0, 'frameTime',0, 'cycleTime',0, ...
                          'zStackDirectionModeName',0, 'zUseModeName',0, ...
                          'xGalvoMovementModeName',0, 'SizeX',0, 'SizeY',0,...
                          'SizeZ',0, 'PixelSizeX',0, 'PixelSizeY',0, ...
                          'PixelSizeZ',0);

%Compare the settings that are numbers
% Method description: Find the difference in settings between each pair of 
% datasets. Nonzero entries indicate nonmatching settings, so if the
% nonzero array is empty, all our settings match.
comparedSettings.pixelDwellTime = isempty(nonzeros(diff([rawSettings.pixelDwellTime])));
comparedSettings.lineAccumulation = isempty(nonzeros(diff([rawSettings.lineAccumulation])));
comparedSettings.frameAccumulation = isempty(nonzeros(diff([rawSettings.frameAccumulation])));
comparedSettings.frameAverage = isempty(nonzeros(diff([rawSettings.frameAverage])));
comparedSettings.emWaveForPinAiryCalc = isempty(nonzeros(diff([rawSettings.emWaveForPinAiryCalc])));
comparedSettings.pinholeAiry = isempty(nonzeros(diff(round([rawSettings.pinholeAiry],settingsTolerances.defaultDecimals))));
comparedSettings.pinhole = isempty(nonzeros(diff(round([rawSettings.pinhole],settingsTolerances.pinholeDecimals))));  %Round by default tolerance to prevent false negatives
comparedSettings.rotatorAngle = isempty(nonzeros(diff([rawSettings.rotatorAngle])));
comparedSettings.phaseX = isempty(nonzeros(diff([rawSettings.phaseX])));
comparedSettings.scanDirectionX = isempty(nonzeros(diff([rawSettings.scanDirectionX])));
comparedSettings.zoomSetting = isempty(nonzeros(diff(round([rawSettings.zoomSetting],settingsTolerances.defaultDecimals))));%Round by default tolerance to prevent false negatives
comparedSettings.scanSpeed = isempty(nonzeros(diff([rawSettings.scanSpeed])));
comparedSettings.refractionIndex = isempty(nonzeros(diff([rawSettings.refractionIndex])));
comparedSettings.numericalAperture = isempty(nonzeros(diff([rawSettings.numericalAperture])));
comparedSettings.objectiveNumber = isempty(nonzeros(diff([rawSettings.objectiveNumber])));
comparedSettings.magnification = isempty(nonzeros(diff([rawSettings.magnification])));
comparedSettings.frameTime = isempty(nonzeros(diff(round([rawSettings.frameTime],settingsTolerances.defaultDecimals))));
% ComparedSettings.completeTime = isempty(nonzeros(diff([RawSettings.completeTime])));
comparedSettings.cycleTime = isempty(nonzeros(diff(round([rawSettings.cycleTime],settingsTolerances.defaultDecimals))));    %Round by default tolerance to prevent false negatives
comparedSettings.SizeX = isempty(nonzeros(diff([rawSettings.SizeX])));
comparedSettings.SizeY = isempty(nonzeros(diff([rawSettings.SizeY])));
comparedSettings.SizeZ = isempty(nonzeros(diff([rawSettings.SizeZ])));
comparedSettings.PixelSizeX = isempty(nonzeros(diff(round([rawSettings.PixelSizeX],settingsTolerances.pixelXYDecimals)))); %Round by default tolerance to prevent false negatives
comparedSettings.PixelSizeY = isempty(nonzeros(diff(round([rawSettings.PixelSizeY],settingsTolerances.pixelXYDecimals)))); %Round by default tolerance to prevent false negatives
comparedSettings.PixelSizeZ = isempty(nonzeros(diff(round([rawSettings.PixelSizeZ],settingsTolerances.defaultDecimals)))); %Round by default tolerance to prevent false negatives

%Compare the settings that are strings
% Method desctiption: Compare all settings to the first setting to find the
% the strings that DO NOT match. Nonzero entries indicate that one or more
% settings don't match the first setting string, so if the nonzero array is
% empty, all our settings match.
comparedSettings.scanDirectionXName = isempty(nonzeros(~strcmp({rawSettings.scanDirectionXName},rawSettings(1).scanDirectionXName)));
comparedSettings.immersion = isempty(nonzeros(~strcmp({rawSettings.immersion}, rawSettings(1).immersion)));
comparedSettings.zStackDirectionModeName = ... 
    isempty(nonzeros(~strcmp({rawSettings.zStackDirectionModeName},rawSettings(1).zStackDirectionModeName)));
comparedSettings.zUseModeName = ... 
    isempty(nonzeros(~strcmp({rawSettings.zUseModeName},rawSettings(1).zUseModeName)));
comparedSettings.xGalvoMovementModeName = ... 
    isempty(nonzeros(~strcmp({rawSettings.xGalvoMovementModeName},rawSettings(1).xGalvoMovementModeName)));

