function schnitzcells = CalculateNuclearMovement(Prefix)
%% Load necessary experiment parameters and schnitzcells 
liveExperiment = LiveExperiment(Prefix);

FrameInfo = getFrameInfo(liveExperiment);
FrameTimes = [FrameInfo(:).Time]; % in seconds
nc_info = [liveExperiment.nc9, liveExperiment.nc10, liveExperiment.nc11,...
    liveExperiment.nc12, liveExperiment.nc13, liveExperiment.nc14, length(FrameInfo)];
PixelSize = liveExperiment.pixelSize_um;
nucleusDiameters = zeros(1, 6);
for nc=9:14
    nucleusDiameters(nc-8) = getDefaultParameters(FrameInfo,[['d', num2str(nc)]])/PixelSize; % in pixels
end
xDim = liveExperiment.xDim;
yDim = liveExperiment.yDim;
zDim = liveExperiment.zDim;
schnitzcells = DetermineNucleiEndFrames(Prefix);


%% Calculate movement of each schnitz cell over each nuclear cycle and each
% frame for every nucleus and plot velocities and 

% Add Velocity and Displacement Info for each schnitz cell with the
% following fields:
% TotalDistanceTraveled (sum of displacements values)in microns
% MeanDistanceTraveledPerSecond in microns/s
% TotalDisplacement (Distance between first and last points of cycle) in
% microns
% MeanDisplacementPerSecond in microns/s
% TotalXDistanceTraveled (sum of x displacement values) in microns
% TotalYDistanceTraveled (sum of y displacement values) in microns
% TotalXDisplacement (sum of x displacement values) in microns
% TotalYDisplacement(sum of y displacement values) in microns
% xVelocity (microns/s)
% yVelocity (micron/s)
% Velocity (microns/s) 

if ~isfield(schnitzcells, 'anaphaseFrame')
    for sc = 1:length(schnitzcells)
        schnitzcells(sc).anaphaseFrame = [];
        schnitzcells(sc).inferredAnaphaseFrame = false;
    end
end


for schnitz_index=1:length(schnitzcells)
    CurrentSchnitz = schnitzcells(schnitz_index);
    NFrames = length(CurrentSchnitz.frames);
    % Initialize VelocityInfo  
    VelocityInfo = {};
    VelocityInfo.TotalXDistanceTraveled = NaN;
    VelocityInfo.TotalYDistanceTraveled = NaN;
    VelocityInfo.TotalDistanceTraveled = NaN;
    VelocityInfo.MeanDistanceTraveledPerSecond = NaN;
    VelocityInfo.TotalXDisplacement = NaN;
    VelocityInfo.TotalYDisplacement = NaN;
    VelocityInfo.TotalDisplacement = NaN;
    VelocityInfo.MeanDisplacementPerSecond = NaN;
    VelocityInfo.xVelocity = NaN(1, NFrames-1);
    VelocityInfo.yVelocity = NaN(1, NFrames-1);
    VelocityInfo.Velocity = NaN(1, NFrames-1);
    VelocityInfo.SchnitzCellTimes = FrameTimes(CurrentSchnitz.frames);
    VelocityInfo.SchnitzHasAllFrames = false;
    VelocityInfo.SchnitzHasFirstAndLastCycleFrames = false;
    VelocityInfo.SchnitzCycleDuration = NaN;
    if isempty(CurrentSchnitz.cycle)
        CurrentNC = NaN;
        schnitzcells(schnitz_index).VelocityInfo = VelocityInfo;
        continue
    else
        CurrentNC = CurrentSchnitz.cycle;
    end
    
    
    
    % First add total info for cells with info for the first and last franes
    % of the cycle. Note that for NC14, the cell must be present in the
    % last frame of the movie rather than the last frame of the cycle. 

    
    SchnitzHasFirstAndLastCycleFrames = false;
    if ~isempty(CurrentSchnitz.anaphaseFrame) & ~isempty(CurrentSchnitz.inferredAnaphaseFrame)
        if ~CurrentSchnitz.inferredAnaphaseFrame & ~isnan(CurrentSchnitz.anaphaseFrame)
            schnitzcells(schnitz_index).containsFirstFrameOfCycle = 1;
            if CurrentSchnitz.containsLastFrameOfCycle == 1
                SchnitzHasFirstAndLastCycleFrames = true;
            end
        end     
    end
    
    % If only one frame exists, it is not possible to calculate any
    % velocity info  
    if length(CurrentSchnitz.frames) <= 1
        schnitzcells(schnitz_index).VelocityInfo = VelocityInfo;
        continue
    end
    VelocityInfo.SchnitzHasFirstAndLastCycleFrames =SchnitzHasFirstAndLastCycleFrames;
    
    SchnitzHasAllFrames = false;
    FrameDifferences = diff(CurrentSchnitz.frames);
    if SchnitzHasFirstAndLastCycleFrames
        if isempty(find(FrameDifferences ~= 1, 1))
            SchnitzHasAllFrames = true;
        end
    end
    VelocityInfo.SchnitzHasAllFrames = SchnitzHasAllFrames;
    % First calculate TotalXDistanceTraveled, TotalYDistanceTraveled,
    % TotalDistanceTraveled and MeanDistanceTraveledPerSecond for schnitz
    % cells with all frames in cycle accounted for 
    DeltaXValuesInPixels = diff(CurrentSchnitz.cenx);
    DeltaYValuesInPixels = diff(CurrentSchnitz.ceny);
    DeltaSValuesInPixels = sqrt(DeltaXValuesInPixels.^2+ DeltaYValuesInPixels.^2);
    DeltaTs = diff(VelocityInfo.SchnitzCellTimes);
    
    if SchnitzHasAllFrames
        VelocityInfo.TotalXDistanceTraveled = sum(abs(DeltaXValuesInPixels))*PixelSize; % in microns
        VelocityInfo.TotalYDistanceTraveled = sum(abs(DeltaYValuesInPixels))*PixelSize; % in microns
        VelocityInfo.TotalDistanceTraveled = sum(DeltaSValuesInPixels)*PixelSize; % in microns
        VelocityInfo.MeanDistanceTraveledPerSecond = VelocityInfo.TotalDistanceTraveled/(VelocityInfo.SchnitzCellTimes(end)-VelocityInfo.SchnitzCellTimes(1)); % in microns per second
    end
    
    % Add TotalDisplacementInfo for schnitz cells with first and last cycle
    % frames accounted for
    if SchnitzHasFirstAndLastCycleFrames
        VelocityInfo.TotalXDisplacement = (CurrentSchnitz.cenx(end)-CurrentSchnitz.cenx(1))*PixelSize;  % in microns
        VelocityInfo.TotalYDisplacement = (CurrentSchnitz.ceny(end)-CurrentSchnitz.ceny(1))*PixelSize;  % in microns
        VelocityInfo.TotalDisplacement = sqrt(VelocityInfo.TotalXDisplacement^2 + VelocityInfo.TotalYDisplacement^2);  % in microns
        VelocityInfo.MeanDisplacementPerSecond = VelocityInfo.TotalDisplacement/(VelocityInfo.SchnitzCellTimes(end)-VelocityInfo.SchnitzCellTimes(1));% in microns per second
        if CurrentNC ~= 14
            VelocityInfo.SchnitzCycleDuration = VelocityInfo.SchnitzCellTimes(end)-VelocityInfo.SchnitzCellTimes(1); % in seconds
        end
    end


    % Finally add info per frame for all adjacent frames 
    
    xVelocitiesInPixels = DeltaXValuesInPixels./DeltaTs;
    yVelocitiesInPixels = DeltaYValuesInPixels./DeltaTs;
    for frame_index = 1:(NFrames-1)
        if FrameDifferences(frame_index) == 1
            VelocityInfo.xVelocity(frame_index) = PixelSize*xVelocitiesInPixels(frame_index);
            VelocityInfo.yVelocity(frame_index) = PixelSize*yVelocitiesInPixels(frame_index);
            VelocityInfo.Velocity(frame_index) = sqrt(VelocityInfo.xVelocity(frame_index)^2 + VelocityInfo.yVelocity(frame_index)^2);
        end   
    end
    
    schnitzcells(schnitz_index).VelocityInfo = VelocityInfo;
end

%% Save schnitz changes
schnitzPath = [liveExperiment.resultsFolder, filesep, Prefix, '_lin.mat'];
if whos(var2str(schnitzcells)).bytes < 2E9
    save(schnitzPath, 'schnitzcells', '-v6');
else
    save(schnitzPath, 'schnitzcells', '-v7.3', '-nocompression');
end


