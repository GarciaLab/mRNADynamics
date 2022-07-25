function [CompiledParticles] = AddQCInfoToCompiledParticles(Prefix, CompiledParticles, FluoString, varargin)

UseZInfo = true;
UseHistoneInfo = false;
UsePositionInfoForApprovedParticles = true;
UsePositionInfo = false;
UseFluoInfo = false;


x = 1;
while x <= length(varargin)
    if strcmp(lower(varargin{x}), 'nozinfo')
        UseZInfo = false;
    elseif strcmp(lower(varargin{x}), 'nohistone')
        UseHistoneInfo = false;
    elseif strcmp(lower(varargin{x}), 'noposition')
        UsePositionInfoForApprovedParticles = false;
    elseif strcmp(lower(varargin{x}), 'nofluo')
        UseFluoInfo = false;
    end
    x = x+1;
    
end
%% Add histone fluorescence information for each schnitz cell z-stack
liveExperiment = LiveExperiment(Prefix);
DropboxFolder = liveExperiment.userResultsFolder;
schnitzcells = getSchnitzcells(liveExperiment);
if ~strcmpi(liveExperiment.experimentType, '2spot')
    Use2SpotPositionInfo = false;
end

NoHistone = true;
for ch = 1:length(liveExperiment.Channels)
    if contains(lower(liveExperiment.Channels{ch}), 'his')
        NoHistone = false;
    end
end

if NoHistone
    UseHistoneInfo = false;
end


FrameInfo = getFrameInfo(liveExperiment);
FrameTimes = [FrameInfo(:).Time];
nc_info = [liveExperiment.nc9, liveExperiment.nc10, liveExperiment.nc11,...
    liveExperiment.nc12, liveExperiment.nc13, liveExperiment.nc14, length(FrameInfo)];
ncFrames = [0 0 0 0 0 0 0 0 nc_info(1:6)];
PixelSize = liveExperiment.pixelSize_um;
SnippetSize = liveExperiment.snippetSize_px;
nucleusDiameters = zeros(1, 6);
for nc=9:14
    nucleusDiameters(nc-8) = getDefaultParameters(FrameInfo,[['d', num2str(nc)]])/PixelSize; % in pixels
end
xDim = liveExperiment.xDim;
yDim = liveExperiment.yDim;
zDim = liveExperiment.zDim;

if ~isfield(schnitzcells, 'HistoneFluo')
    schnitzcells = integrateHistoneFluo(Prefix, schnitzcells, FrameInfo);
end

%  Add spot fluorescence information for each particle z-stack in CompiledParticles

ParticleSchnitzDistances = GetParticleSchnitzDistanceDist(schnitzcells, CompiledParticles, xDim, yDim);
ParticleFluos = GetParticleFluoDist(CompiledParticles);
MaxParticleFluo = max(ParticleFluos);
MinSchnitzDistanceToBoundary = prctile(ParticleSchnitzDistances, 95);
[CompiledParticles] = GetFluoZInfo(liveExperiment, CompiledParticles);





%%
his_zpadding = 1;
spot_zpadding = 3;
ChN = 1;
EarliestTurnOnTimeCells = cell(1, 6);
EarliestTurnOnTimes = NaN(1, 6);
for NC = 9:14
    EarliestTurnOnTimeCells{NC-8} = [];
    for i =1:length(CompiledParticles{ChN})
        if schnitzcells(CompiledParticles{ChN}(i).schnitz).cycle == NC & ...
                min(CompiledParticles{ChN}(i).Frame) >= nc_info(NC-8)-1
            EarliestTurnOnTimeCells{NC-8}(end+1) = min(CompiledParticles{ChN}(i).Frame);
        end
    end
    if ~isempty(EarliestTurnOnTimeCells{NC-8})
        EarliestTurnOnTimes(NC-8) = min(EarliestTurnOnTimeCells{NC-8});
    end
end


for ChN=1:length(liveExperiment.spotChannels)
   
    for i =1:length(CompiledParticles{ChN})
        CompiledParticles{ChN}(i).schnitzcell = schnitzcells(CompiledParticles{ChN}(i).schnitz);
        CompiledParticles{ChN}(i).cycle = CompiledParticles{ChN}(i).schnitzcell.cycle;
        CompiledParticles{ChN}(i).FlaggingInfo = {};
        if isfield(CompiledParticles{ChN}, 'ManualFrameApproved')
            CompiledParticles{ChN}(i).FrameApproved = CompiledParticles{ChN}(i).ManualFrameApproved;
        end
        if isfield(CompiledParticles{ChN}, 'ManualApproved')
            CompiledParticles{ChN}(i).Approved = CompiledParticles{ChN}(i).ManualApproved;
        end
    end
    %%
    
    for i=1:length(CompiledParticles{ChN})%[test_idx]
        p = CompiledParticles{ChN}(i);
        
        sc = p.schnitzcell;
        if sc.Flag == 6
            p.FlaggingInfo.SickNucleus = true;
        else
            p.FlaggingInfo.SickNucleus = false;
        end
        
        if ~isfield(p, 'Approved')
            p.Approved = 0;
        end
        if ~isfield(sc, 'Approved')
            sc.Approved = 0;
        end
        if p.Approved == 1 & sc.Approved == 1
            p.FlaggingInfo.ApprovedNucleus = true;
        else
            p.FlaggingInfo.ApprovedNucleus = false;
        end
        
        
        
        % Find True Nuclear start frame by inference or direct observation
        if ~isempty(sc.inferredAnaphaseFrame)
            if sc.inferredAnaphaseFrame == 0
                if ~isempty(sc.anaphaseFrame) & sc.anaphaseFrame ~= 0
                    if sc.anaphaseFrame >= nc_info(sc.cycle-8)-2
                        p.FlaggingInfo.FirstNucFrame = min(sc.anaphaseFrame, min(p.Frame));
                        p.FlaggingInfo.FNFinferred = false;
                    else
                        p.FlaggingInfo.FirstNucFrame=  min([sc.frames.' p.Frame,  max([1, nc_info(sc.cycle-8)])]);
                        p.FlaggingInfo.FNFinferred = true;
                    end
                else
                    p.FlaggingInfo.FirstNucFrame =  min([sc.frames.' p.Frame,  max([1, nc_info(sc.cycle-8)])]);
                    p.FlaggingInfo.FNFinferred = true;
                end
            elseif sc.anaphaseFrame >= nc_info(sc.cycle-8)-2 & (nc_info(sc.cycle-8) > 0)
                p.FlaggingInfo.FirstNucFrame =  sc.anaphaseFrame;
                p.FlaggingInfo.FNFinferred = true;
            else
                p.FlaggingInfo.FirstNucFrame =   min([sc.frames.' p.Frame max([1, nc_info(sc.cycle-8)])]);
                p.FlaggingInfo.FNFinferred = true;
            end
        elseif ~isempty(sc.anaphaseFrame)
            if sc.anaphaseFrame >= nc_info(sc.cycle-8)-2
                p.FlaggingInfo.FirstNucFrame = min([max([1, sc.anaphaseFrame]) sc.frames.' p.Frame, max([1, nc_info(sc.cycle-8)])]);
                p.FlaggingInfo.FNFinferred = true;
            else
                p.FlaggingInfo.FirstNucFrame =  min([sc.frames.' p.Frame  max([1, nc_info(sc.cycle-8)])]);
                p.FlaggingInfo.FNFinferred = true;
            end
        else
            p.FlaggingInfo.FirstNucFrame =  min([sc.frames.' p.Frame  max([1, nc_info(sc.cycle-8)])]);
            p.FlaggingInfo.FNFinferred = true;
        end
        % Find True Nuclear end frame by inference or direct observation
        
        p.FlaggingInfo.LastNucFrame = max([1 p.NucEnd sc.frames.' p.Frame max([1 nc_info(sc.cycle-7)])]);
        if p.FlaggingInfo.LastNucFrame ~= max(sc.frames.')
            p.FlaggingInfo.FNEinferred = true;
        else
            p.FlaggingInfo.FNEinferred = false;
        end
        
        p.FlaggingInfo.AllSchnitzFrames =...
            p.FlaggingInfo.FirstNucFrame:p.FlaggingInfo.LastNucFrame;
        
        p.FlaggingInfo.FluoVector = NaN(1, length(p.FlaggingInfo.AllSchnitzFrames));
        p.FlaggingInfo.SpotStateCanBeCalled = false;
        p.FlaggingInfo.SpotStateDefinitive = zeros(1, length(p.FlaggingInfo.AllSchnitzFrames), 'logical');
        p.FlaggingInfo.ParticleFrameApproved = -1*ones(1, length(p.FlaggingInfo.AllSchnitzFrames));
        p.FlaggingInfo.SchnitzFrameApproved = -1*ones(1, length(p.FlaggingInfo.AllSchnitzFrames));
        p.FlaggingInfo.SchnitzPresent = zeros(1, length(p.FlaggingInfo.AllSchnitzFrames), 'logical');
        p.FlaggingInfo.SchnitzAwayFromBoundary = zeros(1, length(p.FlaggingInfo.AllSchnitzFrames), 'logical');
        p.FlaggingInfo.SchnitzOn = zeros(1, length(p.FlaggingInfo.AllSchnitzFrames), 'logical');
        p.FlaggingInfo.SchnitzOff = zeros(1, length(p.FlaggingInfo.AllSchnitzFrames), 'logical');
        p.FlaggingInfo.SchnitzTemporarilyOff = zeros(1, length(p.FlaggingInfo.AllSchnitzFrames), 'logical');
        p.FlaggingInfo.SchnitzFinishedTranscribing = zeros(1, length(p.FlaggingInfo.AllSchnitzFrames), 'logical');
        p.FlaggingInfo.UseTraceFluo = zeros(1, length(p.FlaggingInfo.AllSchnitzFrames), 'logical');
        
        
        NFrames = length(p.FlaggingInfo.AllSchnitzFrames);
        
        if all(ismember(p.Frame, p.FlaggingInfo.AllSchnitzFrames))
            p.FlaggingInfo.FluoVector(find( ismember(p.FlaggingInfo.AllSchnitzFrames, p.Frame))) = p.(FluoString);
            try
                p.FlaggingInfo.ParticleFrameApproved(find( ismember(p.FlaggingInfo.AllSchnitzFrames, p.Frame))) = p.FrameApproved;
            catch
                p.FlaggingInfo.ParticleFrameApproved(~isnan(p.FlaggingInfo.FluoVector)) = true;
            end
        else
            SubFluoVector =  p.(FluoString);
            SubFluoVector = SubFluoVector(ismember(p.Frame, p.FlaggingInfo.AllSchnitzFrames));
            p.FlaggingInfo.FluoVector(find( ismember(p.FlaggingInfo.AllSchnitzFrames, p.Frame))) = SubFluoVector;
            SubFrameApprovedVector = p.FrameApproved;
            SubFrameApprovedVector = SubFrameApprovedVector(ismember(p.Frame, p.FlaggingInfo.AllSchnitzFrames));
            p.FlaggingInfo.ParticleFrameApproved(find( ismember(p.FlaggingInfo.AllSchnitzFrames, p.Frame))) = SubFrameApprovedVector;
        end
        
        if all(ismember(sc.frames.', p.FlaggingInfo.AllSchnitzFrames))
            p.FlaggingInfo.SchnitzPresent( ismember(p.FlaggingInfo.AllSchnitzFrames, sc.frames.')) = true;
            p.FlaggingInfo.SchnitzFrameApproved(ismember(p.FlaggingInfo.AllSchnitzFrames, sc.frames.')) = sc.FrameApproved;
            
        else
            p.FlaggingInfo.SchnitzPresent(find( ismember(p.FlaggingInfo.AllSchnitzFrames,  sc.frames.'))) = true;
            SchnitzFrameApproved = sc.FrameApproved;
            SubFrameApprovedVector = SchnitzFrameApproved(ismember(sc.frames.', p.FlaggingInfo.AllSchnitzFrames));
            p.FlaggingInfo.SchnitzFrameApproved(find( ismember(p.FlaggingInfo.AllSchnitzFrames, sc.frames.'))) = SubFrameApprovedVector;
        end
        
        p.FlaggingInfo.UseTraceFluo = ~isnan(p.FlaggingInfo.FluoVector) & p.FlaggingInfo.ParticleFrameApproved;
        
        MaxSpotFluoLevel = NaN(1, NFrames);
        MaxSpotZPos= NaN(1, NFrames);
        SpotZTraceApproved = true(1, NFrames, 'logical');
        
        
        for j = 1:length(p.Frame)
            p_frame_index = find(p.FlaggingInfo.AllSchnitzFrames == p.Frame(j), 1);
            if ~p.FlaggingInfo.UseTraceFluo(p_frame_index)
                continue
            end
            SpotZSliceVector = p.FluoZInfo(j,:);
            MaxF = max(SpotZSliceVector);
            if ~isempty(MaxF)
                MaxSpotFluoLevel(p_frame_index) = MaxF;
                MaxZ = find(round(SpotZSliceVector(spot_zpadding+1:zDim-spot_zpadding),4) == round(MaxF,4), 1);
                if ~isempty(MaxZ)
                    MaxSpotZPos(p_frame_index) = MaxZ+spot_zpadding;
                    SpotZTraceApproved(p_frame_index) = true;
                else
                    MaxZ = find(SpotZSliceVector == MaxF, 1);
                    if ~isempty(MaxZ)
                        MaxSpotZPos(p_frame_index) = MaxZ;
                        SpotZTraceApproved(p_frame_index) = false;
                    end
                end
                
            end
            
        end
        p.FlaggingInfo.MaxSpotFluoLevel = MaxSpotFluoLevel;
        p.FlaggingInfo.MaxSpotZPos = MaxSpotZPos;
        p.FlaggingInfo.SpotZTraceApproved = SpotZTraceApproved;
        
        % Add spot position deviation velocity flags
        
        if UseZInfo
            p.FlaggingInfo.UseTraceFluo =  p.FlaggingInfo.UseTraceFluo  &  p.FlaggingInfo.SpotZTraceApproved;
        end
        
        
        
        
        for frame_index = 1:length(p.FlaggingInfo.AllSchnitzFrames)
            if p.FlaggingInfo.SchnitzPresent(frame_index)
                sc_frame_index = find(sc.frames.' == p.FlaggingInfo.AllSchnitzFrames(frame_index));
                xpos = sc.cenx(sc_frame_index);
                ypos = sc.ceny(sc_frame_index);
                if (xpos > MinSchnitzDistanceToBoundary) & (xpos < xDim- MinSchnitzDistanceToBoundary) & ...
                        (ypos > MinSchnitzDistanceToBoundary) & (ypos < yDim- MinSchnitzDistanceToBoundary)
                    p.FlaggingInfo.SchnitzAwayFromBoundary(frame_index) = true;
                end
            end
        end
        
        % Fill in missing Frames to trace
        ParticleIsOn = ~isnan(p.FlaggingInfo.FluoVector) & p.FlaggingInfo.ParticleFrameApproved;
        NumFrames = length(p.FlaggingInfo.AllSchnitzFrames);
        FrameIsApproved = p.FlaggingInfo.ParticleFrameApproved;
        
        if sc.cycle >= 12
            for p_frame_index = 2:NumFrames-1
                if ParticleIsOn(p_frame_index) & ~p.FlaggingInfo.SpotZTraceApproved(p_frame_index) &(p.FlaggingInfo.FluoVector(p_frame_index)/MaxParticleFluo < 0.2)
                    nearest_off_low = find((isnan(p.FlaggingInfo.FluoVector) |~FrameIsApproved) & (1:NumFrames < p_frame_index), 1, 'last');
                    nearest_off_high = find((isnan(p.FlaggingInfo.FluoVector) |~FrameIsApproved) &  (1:NumFrames > p_frame_index), 1);
                    if isempty(nearest_off_low)
                        nearest_off_low = 0;
                    end
                    if isempty(nearest_off_high)
                        nearest_off_high = length(p.FlaggingInfo.AllSchnitzFrames)+1;
                    end
                    
                    nearest_approved_low = find(ParticleIsOn & p.FlaggingInfo.SpotZTraceApproved & (1:NumFrames > nearest_off_low) & (1:NumFrames < p_frame_index), 1, 'last');
                    nearest_approved_high = find(ParticleIsOn & p.FlaggingInfo.SpotZTraceApproved & (1:NumFrames < nearest_off_high) & (1:NumFrames > p_frame_index), 1);
                    if isempty(nearest_approved_low) & isempty(nearest_approved_high)
                        if nearest_off_high-nearest_off_low + 1 < 10
                            ParticleIsOn(p_frame_index) = false;
                        end
                    elseif isempty(nearest_approved_low) | isempty(nearest_approved_high)
                        if isempty(nearest_approved_low)
                            nearest_approved_low = nearest_off_low;
                        end
                        if isempty(nearest_approved_high)
                            nearest_approved_high = nearest_off_high;
                        end
                        
                        if nearest_approved_high-nearest_approved_low + 1 > 10
                            ParticleIsOn(p_frame_index) = false;
                        end
                        
                        
                    end
                    
                    
                    
                    
                end
            end
            
        end
        
        
        
        
        for p_frame_index = 2:NumFrames-1
            if ~ParticleIsOn(p_frame_index)
                if ParticleIsOn(p_frame_index -1) && ParticleIsOn(p_frame_index+1)
                    ParticleIsOn(p_frame_index) = true;
                else
                    nearest_approved_low = find(ParticleIsOn & (1:NumFrames < p_frame_index), 1, 'last');
                    nearest_approved_high = find(ParticleIsOn & (1:NumFrames > p_frame_index), 1);
                    if ~isempty(nearest_approved_low) & ~isempty(nearest_approved_high)
                        MiddleFramesNotApproved = sum(~FrameIsApproved & (1:NumFrames > nearest_approved_low)& (1:NumFrames < nearest_approved_high));
                        if MiddleFramesNotApproved <= 3 & (nearest_approved_high-nearest_approved_low == MiddleFramesNotApproved -1)
                            ParticleIsOn(p_frame_index) = true;
                        end
                    end
                end
            end
        end
        
        if ParticleIsOn(1) & ~ParticleIsOn(2)
            ParticleIsOn(1) = false;
        end
        if ParticleIsOn(end) & ~ParticleIsOn(end-1)
            ParticleIsOn(end) = false;
        end
        for p_frame_index = 2:NumFrames-1
            if ParticleIsOn(p_frame_index)
                if ~ParticleIsOn(p_frame_index -1) && ~ParticleIsOn(p_frame_index+1)
                    ParticleIsOn(p_frame_index) = false;
                end
            end
        end
        
        p.FlaggingInfo.SchnitzOn = ParticleIsOn;
        p.FlaggingInfo.UseTraceFluo  = p.FlaggingInfo.UseTraceFluo & p.FlaggingInfo.SchnitzOn &...
            ~isnan(p.FlaggingInfo.FluoVector) & p.FlaggingInfo.FluoVector >= 0;
        
        if p.Approved >= 1
            if UsePositionInfoForApprovedParticles
                if all(p.FlaggingInfo.SchnitzAwayFromBoundary(p.FlaggingInfo.AllSchnitzFrames >=EarliestTurnOnTimes(sc.cycle-8)) == 1)
                    p.FlaggingInfo.SpotStateCanBeCalled  = true;
                end
                if all(~p.FlaggingInfo.SchnitzOn(p.FlaggingInfo.SchnitzAwayFromBoundary == 1))
                    p.FlaggingInfo.SchnitzOff(p.FlaggingInfo.SchnitzAwayFromBoundary == 1) = true;
                    p.FlaggingInfo.SpotStateDefinitive(p.FlaggingInfo.SchnitzPresent == 1) = true;
                else
                    for frame_index = 1:NumFrames
                        if ~p.FlaggingInfo.SchnitzOn(frame_index) && p.FlaggingInfo.SchnitzAwayFromBoundary(frame_index)
                            p.FlaggingInfo.SpotStateDefinitive(frame_index) = true;
                            if sum(p.FlaggingInfo.SchnitzOn(frame_index:end) == 1) > 0
                                p.FlaggingInfo.SchnitzTemporarilyOff(frame_index) = true;
                            else
                                p.FlaggingInfo.SchnitzFinishedTranscribing(frame_index) = true;
                            end
                        elseif p.FlaggingInfo.SchnitzOn(frame_index)&& p.FlaggingInfo.SchnitzAwayFromBoundary(frame_index)
                            p.FlaggingInfo.SpotStateDefinitive(frame_index) = true;
                        end
                        
                    end
                    
                    
                end
                
            else
                if all(p.FlaggingInfo.SchnitzPresent(p.FlaggingInfo.AllSchnitzFrames >=EarliestTurnOnTimes(sc.cycle-8)) == 1)
                    p.FlaggingInfo.SpotStateCanBeCalled  = true;
                end
                if all(~p.FlaggingInfo.SchnitzOn(p.FlaggingInfo.SchnitzPresent == 1))
                    p.FlaggingInfo.SchnitzOff(p.FlaggingInfo.SchnitzPresent == 1) = true;
                    p.FlaggingInfo.SpotStateDefinitive(p.FlaggingInfo.SchnitzPresent == 1) = true;
                else
                    for frame_index = 1:NumFrames
                        if ~p.FlaggingInfo.SchnitzOn(frame_index) && p.FlaggingInfo.SchnitzPresent(frame_index)
                            p.FlaggingInfo.SpotStateDefinitive(frame_index) = true;
                            if sum(p.FlaggingInfo.SchnitzOn(frame_index:end) == 1) > 0
                                p.FlaggingInfo.SchnitzTemporarilyOff(frame_index) = true;
                            else
                                p.FlaggingInfo.SchnitzFinishedTranscribing(frame_index) = true;
                            end
                        elseif p.FlaggingInfo.SchnitzOn(frame_index)&& p.FlaggingInfo.SchnitzPresent(frame_index)
                            p.FlaggingInfo.SpotStateDefinitive(frame_index) = true;
                        end
                        
                    end
                    
                    
                end
                
            end
            
        elseif p.Approved < 1
            if UsePositionInfoForApprovedParticles
                if all(p.FlaggingInfo.SchnitzAwayFromBoundary(p.FlaggingInfo.AllSchnitzFrames > EarliestTurnOnTimes(sc.cycle-8) ) == 1)
                    if sum(p.FlaggingInfo.SchnitzOn)/sum(p.FlaggingInfo.SchnitzAwayFromBoundary) < 0.2 & ...
                            (sum(p.FlaggingInfo.FluoVector(~isnan(p.FlaggingInfo.FluoVector)) < 0.3*MaxParticleFluo)/sum(p.FlaggingInfo.SchnitzAwayFromBoundary) < 0.2) & ...
                            (sum(p.FlaggingInfo.FluoVector(~isnan(p.FlaggingInfo.FluoVector)) > 0.2*MaxParticleFluo)/sum(p.FlaggingInfo.SchnitzAwayFromBoundary) < 0.05) & ...
                            (sum(p.FlaggingInfo.FluoVector(~isnan(p.FlaggingInfo.FluoVector)) > 0.3*MaxParticleFluo) == 0)
                        p.FlaggingInfo.SchnitzOff(p.FlaggingInfo.SchnitzAwayFromBoundary) = true;
                        p.FlaggingInfo.SpotStateDefinitive(p.FlaggingInfo.SchnitzAwayFromBoundary) = true;
                        p.FlaggingInfo.SpotStateCanBeCalled= true;
                        p.FlaggingInfo.UseTraceFluo(:) = false;
                        
                    end
                end
            else
                if all(p.FlaggingInfo.SchnitzPresent(p.FlaggingInfo.AllSchnitzFrames > EarliestTurnOnTimes(sc.cycle-8) ) == 1)
                    if sum(p.FlaggingInfo.SchnitzOn)/sum(p.FlaggingInfo.SchnitzPresent) < 0.2 & ...
                            (sum(p.FlaggingInfo.FluoVector(~isnan(p.FlaggingInfo.FluoVector)) < 0.3*MaxParticleFluo)/sum(p.FlaggingInfo.SchnitzPresent) < 0.2) & ...
                            (sum(p.FlaggingInfo.FluoVector(~isnan(p.FlaggingInfo.FluoVector)) > 0.2*MaxParticleFluo)/sum(p.FlaggingInfo.SchnitzPresent) < 0.05)& ...
                            (sum(p.FlaggingInfo.FluoVector(~isnan(p.FlaggingInfo.FluoVector)) > 0.3*MaxParticleFluo) == 0)
                        p.FlaggingInfo.SchnitzOff(p.FlaggingInfo.SchnitzPresent) = true;
                        p.FlaggingInfo.SpotStateDefinitive(p.FlaggingInfo.SchnitzPresent) = true;
                        p.FlaggingInfo.SpotStateCanBeCalled= true;
                        p.FlaggingInfo.UseTraceFluo(:) = false;
                        
                    end
                end
            end
        end
        
        
        
        % Make Histone and spot Fluorescence traces and z positions
        NFrames = length(p.FlaggingInfo.AllSchnitzFrames);
        
        
        
        if UseHistoneInfo
            MaxHisFluoLevel = NaN(1, NFrames);
            MaxHisZPos= NaN(1, NFrames);
            HisTraceApproved = zeros(1, NFrames, 'logical');
            for j = 1:length(sc.frames)
                schnitz_index = find(p.FlaggingInfo.AllSchnitzFrames == sc.frames(j), 1);
                HisVector = sc.HistoneFluo(j,2:zDim+1);
                MaxF = max(HisVector);
                if ~isempty(MaxF)
                    MaxHisFluoLevel(schnitz_index) = MaxF;
                    MaxZ = find(HisVector(his_zpadding:zDim-his_zpadding) == MaxF, 1);
                    if ~isempty(MaxZ)
                        MaxHisZPos(schnitz_index) = MaxZ+his_zpadding;
                        HisTraceApproved(schnitz_index) = true;
                    else
                        MaxZ = find(HisVector == MaxF, 1);
                        if ~isempty(MaxZ)
                            MaxHisZPos(schnitz_index) = MaxZ;
                            HisTraceApproved(schnitz_index) = false;
                        end
                    end
                    
                end
                
            end
            p.FlaggingInfo.MaxHisFluoLevel = MaxHisFluoLevel;
            p.FlaggingInfo.MaxHisZPos = MaxHisZPos;
            p.FlaggingInfo.HisTraceApproved = HisTraceApproved;
            p.FlaggingInfo.UseTraceFluo =  p.FlaggingInfo.UseTraceFluo  &  p.FlaggingInfo.HisTraceApproved;
        else
            
            p.FlaggingInfo.MaxHisFluoLevel =ones(1, NFrames);
            p.FlaggingInfo.MaxHisZPos = ones(1, NFrames);
            p.FlaggingInfo.HisTraceApproved = ones(1, NFrames, 'logical');
        end
        
        
        CompiledParticles{ChN}(i) = p;
        
    end
    
    
    
    %%
    if UsePositionInfo
        error('GM 7/24/22: This function would need to be updated to be compatibile with current filtering.')
        CompiledParticles = AddBoundaryPositionInfo(CompiledParticles, ChN, liveExperiment);
    else
        for i=1:length(CompiledParticles{ChN})%[test_idx]
            p = CompiledParticles{ChN}(i);
            NFrames = length(p.FlaggingInfo.AllSchnitzFrames);
            p.FlaggingInfo.PositionApproved = ones(1, NFrames, 'logical');
            CompiledParticles{ChN}(i) = p;
        end
    end
    
    if UseFluoInfo
        error('GM 7/24/22: This function would need to be updated to be compatibile with current filtering.')
        CompiledParticles = AddFluoFlaggingInfo(CompiledParticles, ChN, FluoString, FrameInfo, liveExperiment);
    else
        for i=1:length(CompiledParticles{ChN})%[test_idx]
            p = CompiledParticles{ChN}(i);
            NFrames = length(p.FlaggingInfo.AllSchnitzFrames);
            p.FlaggingInfo.FluoApproved = ones(1, NFrames, 'logical');
            CompiledParticles{ChN}(i) = p;
        end
    end
    %%
    
    for i=1:length(CompiledParticles{ChN})%
        p = CompiledParticles{ChN}(i);
        sc = p.schnitzcell;
        FI = p.FlaggingInfo;
        CompiledParticles{ChN}(i).ManualApproved = CompiledParticles{ChN}(i).Approved;
        CompiledParticles{ChN}(i).ManualFrameApproved = CompiledParticles{ChN}(i).FrameApproved;
        if FI.SickNucleus | ~FI.ApprovedNucleus
            CompiledParticles{ChN}(i).Approved = -1;
        end
        CompiledParticles{ChN}(i).FrameApproved  = FI.UseTraceFluo(ismember(FI.AllSchnitzFrames, p.Frame));
    end
    

    
    
    
end
%% Also add QC Info to schnitz cells for cells that don't have an associated particle 
clear p
for i=1:length(schnitzcells)
    schnitzcells(i).FlaggingInfo = {};
end
for i=1:length(schnitzcells)
    sc = schnitzcells(i);
    
    if sc.Flag == 6
        sc.FlaggingInfo.SickNucleus = true;
    else
        sc.FlaggingInfo.SickNucleus = false;
    end

    if ~isfield(sc, 'Approved')
        sc.Approved = 0;
    end
    if isempty(sc.Approved)
        sc.Approved = 0;
    end
    
    if  sc.Approved == 1
        sc.FlaggingInfo.ManualApproved = true;
    else
        sc.FlaggingInfo.ManualApproved = false;
    end
    
    
    
    % Find True Nuclear start frame by inference or direct observation
    if ~isempty(sc.inferredAnaphaseFrame)
        if sc.inferredAnaphaseFrame == 0
            if ~isempty(sc.anaphaseFrame) & sc.anaphaseFrame ~= 0
                if sc.anaphaseFrame >= nc_info(sc.cycle-8)-2
                    sc.FlaggingInfo.FirstNucFrame = min(sc.anaphaseFrame);
                    sc.FlaggingInfo.FNFinferred = false;
                else
                    sc.FlaggingInfo.FirstNucFrame=  min([sc.frames.'   max([1, nc_info(sc.cycle-8)])]);
                    sc.FlaggingInfo.FNFinferred = true;
                end
            else
                sc.FlaggingInfo.FirstNucFrame =  min([sc.frames.' max([1, nc_info(sc.cycle-8)])]);
                sc.FlaggingInfo.FNFinferred = true;
            end
        elseif sc.anaphaseFrame >= nc_info(sc.cycle-8)-2 & (nc_info(sc.cycle-8) > 0)
            sc.FlaggingInfo.FirstNucFrame =  sc.anaphaseFrame;
            sc.FlaggingInfo.FNFinferred = true;
        else
            sc.FlaggingInfo.FirstNucFrame =   min([sc.frames.'  max([1, nc_info(sc.cycle-8)])]);
            sc.FlaggingInfo.FNFinferred = true;
        end
    elseif ~isempty(sc.anaphaseFrame)
        if sc.anaphaseFrame >= nc_info(sc.cycle-8)-2
            sc.FlaggingInfo.FirstNucFrame = min([max([1, sc.anaphaseFrame]) sc.frames.' max([1, nc_info(sc.cycle-8)])]);
            sc.FlaggingInfo.FNFinferred = true;
        else
            sc.FlaggingInfo.FirstNucFrame =  min([sc.frames.' max([1, nc_info(sc.cycle-8)])]);
            sc.FlaggingInfo.FNFinferred = true;
        end
    else
        sc.FlaggingInfo.FirstNucFrame =  min([sc.frames.'  max([1, nc_info(sc.cycle-8)])]);
        sc.FlaggingInfo.FNFinferred = true;
    end
    % Find True Nuclear end frame by inference or direct observation
    
    sc.FlaggingInfo.LastNucFrame = max([1  sc.frames.' max([1 nc_info(sc.cycle-7)])]);
    if sc.FlaggingInfo.LastNucFrame ~= max(sc.frames.')
        sc.FlaggingInfo.FNEinferred = true;
    else
        sc.FlaggingInfo.FNEinferred = false;
    end
    
    sc.FlaggingInfo.AllSchnitzFrames =...
        sc.FlaggingInfo.FirstNucFrame:sc.FlaggingInfo.LastNucFrame;
    
    sc.FlaggingInfo.SpotStateCanBeCalled = false;
    sc.FlaggingInfo.SpotStateDefinitive = zeros(1, length(sc.FlaggingInfo.AllSchnitzFrames), 'logical');

    sc.FlaggingInfo.SchnitzFrameApproved = -1*ones(1, length(sc.FlaggingInfo.AllSchnitzFrames));
    sc.FlaggingInfo.SchnitzPresent = zeros(1, length(sc.FlaggingInfo.AllSchnitzFrames), 'logical');
    sc.FlaggingInfo.SchnitzAwayFromBoundary = zeros(1, length(sc.FlaggingInfo.AllSchnitzFrames), 'logical');

    
    
    NFrames = length(sc.FlaggingInfo.AllSchnitzFrames);

    
    if all(ismember(sc.frames.', sc.FlaggingInfo.AllSchnitzFrames))
        sc.FlaggingInfo.SchnitzPresent( ismember(sc.FlaggingInfo.AllSchnitzFrames, sc.frames.')) = true;
        sc.FlaggingInfo.SchnitzFrameApproved(ismember(sc.FlaggingInfo.AllSchnitzFrames, sc.frames.')) = sc.FrameApproved;
        
    else
        sc.FlaggingInfo.SchnitzPresent(find( ismember(sc.FlaggingInfo.AllSchnitzFrames,  sc.frames.'))) = true;
        SchnitzFrameApproved = sc.FrameApproved;
        SubFrameApprovedVector = SchnitzFrameApproved(ismember(sc.frames.', sc.FlaggingInfo.AllSchnitzFrames));
        sc.FlaggingInfo.SchnitzFrameApproved(find( ismember(sc.FlaggingInfo.AllSchnitzFrames, sc.frames.'))) = SubFrameApprovedVector;
    end
    

    
    
    
    for frame_index = 1:length(sc.FlaggingInfo.AllSchnitzFrames)
        if sc.FlaggingInfo.SchnitzPresent(frame_index)
            sc_frame_index = find(sc.frames.' == sc.FlaggingInfo.AllSchnitzFrames(frame_index));
            xpos = sc.cenx(sc_frame_index);
            ypos = sc.ceny(sc_frame_index);
            if (xpos > MinSchnitzDistanceToBoundary) & (xpos < xDim- MinSchnitzDistanceToBoundary) & ...
                    (ypos > MinSchnitzDistanceToBoundary) & (ypos < yDim- MinSchnitzDistanceToBoundary)
                sc.FlaggingInfo.SchnitzAwayFromBoundary(frame_index) = true;
                sc.FlaggingInfo.SpotStateDefinitive(frame_index) = true;
            end
        end
    end
    
    if all(sc.FlaggingInfo.SchnitzAwayFromBoundary(sc.FlaggingInfo.AllSchnitzFrames > EarliestTurnOnTimes(sc.cycle-8) )) & ...
            ~ismember(i, [CompiledParticles{1}(:).schnitz])
        sc.FlaggingInfo.SpotStateCanBeCalled = true;
    end
    
   
    
    schnitzcells(i) = sc;
    
end




if whos(var2str(schnitzcells)).bytes < 2E9
    save([DropboxFolder, filesep, Prefix, filesep, Prefix, '_lin.mat'], 'schnitzcells', '-append')
else
    save([DropboxFolder, filesep, Prefix, filesep, Prefix, '_lin.mat'], 'schnitzcells', '-append', '-nocompression')
end

end
