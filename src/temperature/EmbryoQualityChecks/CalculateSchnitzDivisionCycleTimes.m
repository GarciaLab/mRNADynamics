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
for i = 1:4
    EmbryoNCInfo{i} = [];
end
NCDivisionInfo = NaN(1,4);
DivisionStdInfo  = NaN(1,4);

EmbryoNCInfo2 = cell(1, 4);
for i = 1:4
    EmbryoNCInfo2{i} = [];
end
NCDivisionInfo2 = NaN(1,4);
DivisionStdInfo2  = NaN(1,4);
GoodPrimary = false(1, length(schnitzcells));
GoodSecondary = false(1, length(schnitzcells));

for schnitz_index = 1:length(schnitzcells)
    CurrentSchnitz = schnitzcells(schnitz_index);
    if isfield(CurrentSchnitz, 'FlaggingInfo')
        if all(CurrentSchnitz.FlaggingInfo.SchnitzPresent) & ~isempty(CurrentSchnitz.anaphaseFrame) & ~isempty(CurrentSchnitz.cycle)...
                & ~CurrentSchnitz.inferredAnaphaseFrame
            if CurrentSchnitz.cycle < 14 & CurrentSchnitz.anaphaseFrame >= nc_info(CurrentSchnitz.cycle-8) & ...
                    CurrentSchnitz.anaphaseFrame <= nc_info(CurrentSchnitz.cycle-8)+ 7 &...
                    CurrentSchnitz.FlaggingInfo.AllSchnitzFrames(end) <= nc_info(CurrentSchnitz.cycle-7)+7
                NucleusCycleDuration = (max(FrameTimes(CurrentSchnitz.FlaggingInfo.AllSchnitzFrames)) - ...
                    min(FrameTimes(CurrentSchnitz.FlaggingInfo.AllSchnitzFrames)))/60;
                EmbryoNCInfo{CurrentSchnitz.cycle-9}(end + 1) = ...
                    NucleusCycleDuration;
                GoodPrimary(schnitz_index) = true;
            end
        end
    end
        if CurrentSchnitz.VelocityInfo.SchnitzHasFirstAndLastCycleFrames & ...
                ~isnan(CurrentSchnitz.VelocityInfo.SchnitzCycleDuration) & ...
                ~isempty(CurrentSchnitz.cycle)
            NucleusCycleDuration = CurrentSchnitz.VelocityInfo.SchnitzCycleDuration/60; % in minutes
           GoodSecondary(schnitz_index) = true;
           EmbryoNCInfo2{CurrentSchnitz.cycle-9}(end + 1) = ...
                    NucleusCycleDuration;
%             EmbryoNCInfo{CurrentSchnitz.cycle-9}(end + 1) = ...
%                     NucleusCycleDuration;
     
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


for NC=10:13
    if ~isempty(EmbryoNCInfo2{NC-9})
        if length(EmbryoNCInfo2{NC-9}) >= 10
            NCDivisionInfo2(NC-9) = mean(EmbryoNCInfo2{NC-9});
            DivisionStdInfo2(NC-9)= std(EmbryoNCInfo2{NC-9});
        end
        
    end
end

