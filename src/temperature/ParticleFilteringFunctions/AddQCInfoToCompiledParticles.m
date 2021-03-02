function [CompiledParticles] = AddQCInfoToCompiledParticles(Prefix, CompiledParticles, varargin)

UseZInfo = true;
UseHistoneInfo = true;
UsePositionInfo = false;
Use2SpotPositionInfo = true;
UseFluoInfo = true;
FluoString = 'Fluo';

x = 0;
while x < length(varargin)
    if strcmp(lower(varargin{x}), 'nozinfo')
        UseZInfo = false;
    elseif strcmp(lower(varargin{x}), 'nohistone')
        UseHistoneInfo = false;
    elseif strcmp(lower(varargin{x}), 'noposition')
        UsePositionInfo = false;
    elseif strcmp(lower(varargin{x}), 'no2spot')
        Use2SpotPositionInfo = false;
    elseif strcmp(lower(varargin{x}), 'nofluo')
        UseFluoInfo = false;
    end
    x = x+1;
    
end
%% Add histone fluorescence information for each schnitz cell z-stack
liveExperiment = LiveExperiment(Prefix);
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
PixelSize = liveExperiment.pixelSize_um;
SnippetSize = liveExperiment.snippetSize_px;
nucleusDiameters = zeros(1, 6);
for nc=9:14
    nucleusDiameters(nc-8) = getDefaultParameters(FrameInfo,[['d', num2str(nc)]])/PixelSize; % in pixels
end
xDim = liveExperiment.xDim;
yDim = liveExperiment.yDim;
zDim = liveExperiment.zDim;
schnitzcells = integrateHistoneFluo(Prefix, schnitzcells, FrameInfo);

%  Add spot fluorescence information for each particle z-stack in CompiledParticles


[CompiledParticles] = GetFluoZInfo(liveExperiment, CompiledParticles);
%%
his_zpadding = 1;
spot_zpadding = 2;

for ChN=1:length(liveExperiment.spotChannels)
    for i =1:length(CompiledParticles{ChN})
        CompiledParticles{ChN}(i).schnitzcell = schnitzcells(CompiledParticles{ChN}(i).schnitz);
        CompiledParticles{ChN}(i).cycle = CompiledParticles{ChN}(i).schnitzcell.cycle;
        CompiledParticles{ChN}(i).FlaggingInfo = {};
        if isfield(CompiledParticles{ChN}, 'ManualFrameApproved')
            CompiledParticles{ChN}(i).FrameApproved = CompiledParticles{ChN}(i).ManualFrameApproved;
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
        
        if p.Approved == 1
            p.FlaggingInfo.ApprovedNucleus = true;
        else
            p.FlaggingInfo.ApprovedNucleus = false;
        end
        % Find True Nuclear start frame by inference or direct observation
        if ~isempty(sc.inferredAnaphaseFrame)
            if sc.inferredAnaphaseFrame == 0
                if ~isempty(sc.anaphaseFrame)
                    if sc.anaphaseFrame >= nc_info(sc.cycle-8)-1
                        p.FlaggingInfo.TrueNucStart = min(max([sc.anaphaseFrame, 1]),sc.frames(1));
                        p.FlaggingInfo.TNSinferred = false;
                    else
                        p.FlaggingInfo.TrueNucStart =  min(max([nc_info(sc.cycle-8), 1]), sc.frames(1));
                        p.FlaggingInfo.TNSinferred = true;
                    end
                else
                    p.FlaggingInfo.TrueNucStart =  min(max([1, nc_info(sc.cycle-8)]), sc.frames(1));
                    p.FlaggingInfo.TNSinferred = true;
                end
            elseif sc.anaphaseFrame >= nc_info(sc.cycle-8)-1
                p.FlaggingInfo.TrueNucStart = min(max([1, sc.anaphaseFrame]), sc.frames(1));
                p.FlaggingInfo.TNSinferred = true;
            else
                p.FlaggingInfo.TrueNucStart =  min(max([1,nc_info(sc.cycle-8)]), sc.frames(1));
                p.FlaggingInfo.TNSinferred = true;
            end
        elseif ~isempty(sc.anaphaseFrame)
            if sc.anaphaseFrame >= nc_info(sc.cycle-8)-1
                p.FlaggingInfo.TrueNucStart = min(max([1, sc.anaphaseFrame]), sc.frames(1));
                p.FlaggingInfo.TNSinferred = false;
            else
                p.FlaggingInfo.TrueNucStart =  min(max([1, nc_info(sc.cycle-8)]), sc.frames(1));
                p.FlaggingInfo.TNSinferred = true;
            end
        else
            p.FlaggingInfo.TrueNucStart =  min(max([1, nc_info(sc.cycle-8)]), sc.frames(1));
            p.FlaggingInfo.TNSinferred = true;
        end
        % Find True Nuclear end frame by inference or direct observation
        if p.NucEnd >= nc_info(sc.cycle-7)-1
            p.FlaggingInfo.TrueNucEnd = max([1, p.NucEnd]);
            p.FlaggingInfo.TNEinferred =false;
        else
            p.FlaggingInfo.TrueNucEnd = max([1, nc_info(sc.cycle-7)-1]);
            p.FlaggingInfo.TNEinferred =true;
        end
        
        p.FlaggingInfo.TrueFrames =...
            p.FlaggingInfo.TrueNucStart:p.FlaggingInfo.TrueNucEnd;
        
        
        
        % Make Histone and spot Fluorescence traces and z positions
        NFrames = length(p.FlaggingInfo.TrueFrames);
        if UseZInfo
            MaxSpotFluoLevel = NaN(1, NFrames);
            MaxSpotZPos= NaN(1, NFrames);
            SpotZTraceApproved = zeros(1, NFrames);
            DetectedSpot = zeros(1, NFrames);
            

            for j = 1:length(p.Frame)
                particle_index = find(p.FlaggingInfo.TrueFrames == p.Frame(j), 1);
                if isempty(find(sc.frames == p.Frame(j), 1))
                    continue
                elseif p.FrameApproved(j) == 0
                    continue
                end
                DetectedSpot(particle_index) = 1;
                SpotZSliceVector = p.FluoZInfo(j,:);
                MaxF = max(SpotZSliceVector);
                if ~isempty(MaxF)
                    MaxSpotFluoLevel(particle_index) = MaxF;
                    MaxZ = find(SpotZSliceVector(spot_zpadding:zDim-spot_zpadding) == MaxF, 1);
                    if ~isempty(MaxZ)
                        MaxSpotZPos(particle_index) = MaxZ+spot_zpadding;
                        SpotZTraceApproved(particle_index) = 1;
                    else
                        MaxZ = find(SpotZSliceVector == MaxF, 1);
                        if ~isempty(MaxZ)
                            MaxSpotZPos(particle_index) = MaxZ;
                            SpotZTraceApproved(particle_index) = -1;
                        end
                    end

                end

            end
            p.FlaggingInfo.MaxSpotFluoLevel = MaxSpotFluoLevel;
            p.FlaggingInfo.MaxSpotZPos = MaxSpotZPos;
            p.FlaggingInfo.SpotZTraceApproved = SpotZTraceApproved;

         
            p.FlaggingInfo.DetectedSpot = DetectedSpot;
            % Add spot position deviation velocity flags
        else
            p.FlaggingInfo.MaxSpotFluoLevel = ones(1, NFrames);
            p.FlaggingInfo.MaxSpotZPos = ones(1, NFrames);
            p.FlaggingInfo.SpotZTraceApproved = ones(1, NFrames);

         
            p.FlaggingInfo.DetectedSpot = ones(1, NFrames);
        end
        
        
        if UseHistoneInfo
            MaxHisFluoLevel = NaN(1, NFrames);
            MaxHisZPos= NaN(1, NFrames);
            HisTraceApproved = zeros(1, NFrames);
            DetectedNucleus = zeros(1, NFrames);
            for j = 1:length(sc.frames)
                schnitz_index = find(p.FlaggingInfo.TrueFrames == sc.frames(j), 1);
                DetectedNucleus(schnitz_index) = 1;
                HisVector = sc.HistoneFluo(j,2:zDim+1);
                MaxF = max(HisVector);
                if ~isempty(MaxF)
                    MaxHisFluoLevel(schnitz_index) = MaxF;
                    MaxZ = find(HisVector(his_zpadding:zDim-his_zpadding) == MaxF, 1);
                    if ~isempty(MaxZ)
                        MaxHisZPos(schnitz_index) = MaxZ+his_zpadding;
                        HisTraceApproved(schnitz_index) = 1;
                    else
                        MaxZ = find(HisVector == MaxF, 1);
                        if ~isempty(MaxZ)
                            MaxHisZPos(schnitz_index) = MaxZ;
                            HisTraceApproved(schnitz_index) = -1;
                        end
                    end

                end

            end
            p.FlaggingInfo.MaxHisFluoLevel = MaxHisFluoLevel;
            p.FlaggingInfo.MaxHisZPos = MaxHisZPos;
            p.FlaggingInfo.HisTraceApproved = HisTraceApproved;
            p.FlaggingInfo.DetectedNucleus = DetectedNucleus;
        else

            p.FlaggingInfo.MaxHisFluoLevel =ones(1, NFrames);
            p.FlaggingInfo.MaxHisZPos = ones(1, NFrames);
            p.FlaggingInfo.HisTraceApproved = ones(1, NFrames);
            p.FlaggingInfo.DetectedNucleus = ones(1, NFrames);
        end
        CompiledParticles{ChN}(i) = p;
    end
    
    if Use2SpotPositionInfo
        CompiledParticles = Add2SpotFlaggingInfo(CompiledParticles, ChN, FluoString, FrameInfo, liveExperiment);
    else
        for i=1:length(CompiledParticles{ChN})%[test_idx]
            p = CompiledParticles{ChN}(i);
            NFrames = length(p.FlaggingInfo.TrueFrames);
            p.FlaggingInfo.TwoSpotApproved = ones(1, NFrames);
            CompiledParticles{ChN}(i) = p;
        end
    end
        
    %%
    if UsePositionInfo
        CompiledParticles = AddBoundaryPositionInfo(CompiledParticles, ChN, liveExperiment);
    else
        for i=1:length(CompiledParticles{ChN})%[test_idx]
            p = CompiledParticles{ChN}(i);
            NFrames = length(p.FlaggingInfo.TrueFrames);
            p.FlaggingInfo.PositionApproved = ones(1, NFrames);
            CompiledParticles{ChN}(i) = p;
        end
    end
    
    if UseFluoInfo 
       CompiledParticles = AddFluoFlaggingInfo(CompiledParticles, ChN, FluoString, FrameInfo, liveExperiment);
    else
        for i=1:length(CompiledParticles{ChN})%[test_idx]
            p = CompiledParticles{ChN}(i);
            NFrames = length(p.FlaggingInfo.TrueFrames);
            p.FlaggingInfo.FluoApproved = ones(1, NFrames);
            CompiledParticles{ChN}(i) = p;
        end
    end 
    %%
    
    for i=1:length(CompiledParticles{ChN})%
        p = CompiledParticles{ChN}(i);
        sc = p.schnitzcell;
        FI = p.FlaggingInfo;
        if FI.SickNucleus | ~FI.ApprovedNucleus
            CompiledParticles{ChN}(i).Approved = false;
            CompiledParticles{ChN}(i).ManualApproved = false;
            CompiledParticles{ChN}(i).ManualFrameApproved = CompiledParticles{ChN}(i).FrameApproved;
        else
            NFrames = length(FI.TrueFrames);
            FrameApprovedByFluoValue = FI.FluoApproved;
            FrameApprovedByZPos = FI.SpotZTraceApproved;
            FrameApprovedByHisZPos = FI.HisTraceApproved;
            FrameApprovedByBoundaryPosition = FI.PositionApproved;
            FrameApprovedBy2SpotPos = FI.TwoSpotApproved;
            
            FrameApprovedByAllMetrics = (FrameApprovedByFluoValue ~= -1) & ...
                (FrameApprovedByZPos ~= -1) & (FrameApprovedByHisZPos ~= -1) & ...
                (FrameApprovedByBoundaryPosition ~= -1) & (FrameApprovedBy2SpotPos ~= -1);
            
            CompiledParticles{ChN}(i).FlaggingInfo.FrameApprovedFinal = FrameApprovedByAllMetrics;
            FrameApproved = FrameApprovedByAllMetrics(ismember(FI.TrueFrames, p.Frame));
            CompiledParticles{ChN}(i).ManualFrameApproved = CompiledParticles{ChN}(i).FrameApproved;
            CompiledParticles{ChN}(i).FrameApproved = FrameApproved;
            CompiledParticles{ChN}(i).ManualApproved = CompiledParticles{ChN}(i).Approved;
            
        end
        
       
        

    end
end
