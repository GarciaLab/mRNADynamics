function [NCDivisionInfo,DivisionStdInfo] = CalculateSchnitzDivisionCycleTimes(Prefix)
%% Load necessary info into memory
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
DetermineNucleiEndFrames(Prefix);
schnitzcells = CalculateNuclearMovement(Prefix);

%% Calculate schnitz cell division times
%  EmbryosWithNuclearTrackingChecked = [];

%APResolution = liveExperiment.APResolution;
%numAPbins = uint16(1/APResolution+1);

%Temp_obs = liveExperiment.Temp_obs;
EmbryoNCInfo = cell(1, 4);
NCDivisionInfo = NaN(1,4);
DivisionStdInfo  = NaN(1,4);

for schnitz_index = 1:length(schnitzcells)
    CurrentSchnitz = schnitzcells(schnitz_index);
    if CurrentSchnitz.VelocityInfo.SchnitzHasFirstAndLastCycleFrames & ...
            ~isnan(CurrentSchnitz.VelocityInfo.SchnitzCycleDuration) & ...
            ~isempty(CurrentSchnitz.cycle)
        NucleusCycleDuration = CurrentSchnitz.VelocityInfo.SchnitzCycleDuration/60; % in minutes
        if ~isempty(EmbryoNCInfo{CurrentSchnitz.cycle-9})
            EmbryoNCInfo{CurrentSchnitz.cycle-9}(length(EmbryoNCInfo{CurrentSchnitz.cycle-9}) + 1) = ...
                NucleusCycleDuration;
        else
            EmbryoNCInfo{CurrentSchnitz.cycle-9} = [NucleusCycleDuration];
        end
    end
    
end

%%
for NC=10:13
    if ~isempty(EmbryoNCInfo{NC-9})
        if length(EmbryoNCInfo{NC-9}) >= 10
            NCDivisionInfo(NC-9) = mean(EmbryoNCInfo{NC-9});
            DivisionStdInfo(NC-9)= std(EmbryoNCInfo{NC-9});
        end
        
    end
end


