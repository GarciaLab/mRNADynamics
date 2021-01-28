function schnitzcells = DetermineNucleiEndFrames(Prefix)


%% Load experiment information and schnitz cells
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
schnitzcells = getSchnitzcells(liveExperiment);

%% Distribution of nuclear velocities
% Grabs information about all schnitz cells movement from frame to frame
Frames = [];
cycles = [];
SchnitzIDs = [];
Velocities = [];

for schnitz_index=1:length(schnitzcells)
    CurrentSchnitz = schnitzcells(schnitz_index);
    if length(CurrentSchnitz.frames) <= 1
        continue
    end
    if isempty(CurrentSchnitz.cycle)
        CurrentNC = NaN;
    else
        CurrentNC = CurrentSchnitz.cycle;
    end
    for frame_index = 1:(length(CurrentSchnitz.frames)-1)
        CurrentFrame = CurrentSchnitz.frames(frame_index);
        NextFrame =  CurrentSchnitz.frames(frame_index+1);
        if NextFrame-CurrentFrame ~= 1
            continue
        end
        DeltaX = CurrentSchnitz.cenx(frame_index+1)-CurrentSchnitz.cenx(frame_index);
        DeltaY = CurrentSchnitz.ceny(frame_index+1)-CurrentSchnitz.ceny(frame_index);
        DistanceInPixels = sqrt(DeltaX^2 + DeltaY^2);
        DeltaT = FrameTimes(NextFrame)-FrameTimes(CurrentFrame);
        Velocity = DistanceInPixels*PixelSize/DeltaT;
        Frames = [Frames, uint16(CurrentFrame)];
        cycles = [cycles, uint16(CurrentNC)];
        SchnitzIDs = [SchnitzIDs,uint16(schnitz_index)];
        Velocities = [Velocities, Velocity];
        
        
    end
end


ReasonableVelocityCutoffs = zeros(length(FrameInfo)-1, 6);
for nc =9:14
    for i =1:(length(FrameInfo)-1)
        VelocitySubset = Velocities(Frames == i & cycles == nc);
        if ~isempty(VelocitySubset)
            ReasonableVelocityCutoffs(i, nc-8) = prctile(VelocitySubset, 99);
        end
    end
end


%% Gets information about anaphase frames and positions at which they occur
AnaphaseFrames = [];
Xpositions = [];
Ypositions = [];
CycleInfoForAnaphase = [];

for schnitz_index=1:length(schnitzcells)
    if isfield(schnitzcells, 'inferredAnaphaseFrame')
        if isempty(schnitzcells(schnitz_index).inferredAnaphaseFrame)
            schnitzcells(schnitz_index).inferredAnaphaseFrame = 0;
        end
        
    end
    CurrentSchnitz = schnitzcells(schnitz_index);
    if isempty(CurrentSchnitz.cycle)
        continue
    else
        CurrentNC = CurrentSchnitz.cycle;
    end
    
    
    
    if isfield(schnitzcells, 'anaphaseFrame')
        if ~isempty(CurrentSchnitz.anaphaseFrame)
            if ~isempty(CurrentSchnitz.inferredAnaphaseFrame)
                if CurrentSchnitz.inferredAnaphaseFrame == 0
                    if ~isnan(CurrentSchnitz.anaphaseFrame)
                        CycleInfoForAnaphase = [CycleInfoForAnaphase, CurrentNC];
                        AnaphaseFrames = [AnaphaseFrames, uint16(CurrentSchnitz.anaphaseFrame)];
                        Xpositions = [Xpositions, CurrentSchnitz.cenx(1)];
                        Ypositions = [Ypositions, CurrentSchnitz.ceny(1)];
                    end
                end
            elseif (~isnan(CurrentSchnitz.anaphaseFrame) & (CurrentSchnitz.anaphaseFrame ~= 0))
                CycleInfoForAnaphase = [CycleInfoForAnaphase, CurrentNC];
                AnaphaseFrames = [AnaphaseFrames, uint16(CurrentSchnitz.anaphaseFrame)];
                Xpositions = [Xpositions, CurrentSchnitz.cenx(1)];
                Ypositions = [Ypositions, CurrentSchnitz.ceny(1)];
                
            end
            
        end
        
    end
    
end

%%

for schnitz_index=1:length(schnitzcells)
    schnitzcells(schnitz_index).containsLastFrameOfCycle = 0;
end



for schnitz_index = 1:length(schnitzcells)
    CurrentSchnitz = schnitzcells(schnitz_index);
    CurrentNC = CurrentSchnitz.cycle;
    if isempty(CurrentNC)
        continue
    elseif (CurrentNC == 14) &(CurrentSchnitz.frames(end) == nc_info(end))
        schnitzcells(schnitz_index).containsLastFrameOfCycle = 1;
        continue
    elseif (CurrentNC == 14)
        schnitzcells(schnitz_index).containsLastFrameOfCycle = -1;
        continue
    end
    
    SubsequentCycleAnaphaseFrames = unique(AnaphaseFrames(CycleInfoForAnaphase == CurrentNC+1));
    PossibleEndCycleFrames = SubsequentCycleAnaphaseFrames-1;
    
    SchnitzCellLastFrame = CurrentSchnitz.frames(end);
    if ~ismember(SchnitzCellLastFrame, PossibleEndCycleFrames)
        schnitzcells(schnitz_index).containsLastFrameOfCycle = -1;
    elseif (CurrentSchnitz.Flag == 6) | (CurrentSchnitz.Approved == 0)
        schnitzcells(schnitz_index).containsLastFrameOfCycle = -1;
    elseif SchnitzCellLastFrame == max(PossibleEndCycleFrames)
        schnitzcells(schnitz_index).containsLastFrameOfCycle = 1;
    else
        VelocityCutoff = ReasonableVelocityCutoffs(SchnitzCellLastFrame, CurrentNC-8);
        DeltaT = FrameTimes(SchnitzCellLastFrame+1)-FrameTimes(SchnitzCellLastFrame);
        MaxPixelDistance = DeltaT*VelocityCutoff/PixelSize;
        CurrentSchnitzXpos = CurrentSchnitz.cenx(end);
        CurrentSchnitzYpos = CurrentSchnitz.ceny(end);
        Distances = sqrt((Xpositions-CurrentSchnitzXpos).^2+(Ypositions-CurrentSchnitzYpos).^2);
        SortedDistances = sort(Distances(CycleInfoForAnaphase == CurrentNC +1));
        NearbyAnaphases = AnaphaseFrames((CycleInfoForAnaphase == CurrentNC +1) & ...
            (Distances <= SortedDistances(min(10, length(SortedDistances)))));
        if ismember(SchnitzCellLastFrame+1, NearbyAnaphases)
            if (CurrentSchnitzXpos+MaxPixelDistance <= xDim) & ...
                    (CurrentSchnitzXpos-MaxPixelDistance >= 0) & ...
                    (CurrentSchnitzYpos+MaxPixelDistance <= yDim) & ...
                    (CurrentSchnitzYpos-MaxPixelDistance >= 0)
                schnitzcells(schnitz_index).containsLastFrameOfCycle = 1;
            end
        end
        
    end
    
end


%% Save schnitz changes
schnitzPath = [liveExperiment.resultsFolder, filesep, Prefix, '_lin.mat'];
if whos(var2str(schnitzcells)).bytes < 2E9
    save(schnitzPath, 'schnitzcells', '-v6');
else
    save(schnitzPath, 'schnitzcells', '-v7.3', '-nocompression');
end


