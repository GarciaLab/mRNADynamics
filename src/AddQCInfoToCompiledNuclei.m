function CompiledNuclei = AddQCInfoToCompiledNuclei(Prefix, varargin)

UseZInfo = true;
UseHistoneInfo = true;

x = 1;
while x <= length(varargin)
    if strcmp(lower(varargin{x}), 'nozinfo')
        UseZInfo = false;
    elseif strcmp(lower(varargin{x}), 'nohistone')
        UseHistoneInfo = false;
    end
    x = x+1;
    
end

%% Add histone fluorescence information for each schnitz cell z-stack
liveExperiment = LiveExperiment(Prefix);
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
load([liveExperiment.resultsFolder,filesep,'CompiledNuclei.mat']);
load([liveExperiment.resultsFolder,filesep,'HealthSummary.mat']);

FinishedNuclearTracking = HealthSummary.NuclearTrackingDone;

    

%  Add spot fluorescence information for each particle z-stack in CompiledParticles


for i =1:length(CompiledNuclei)
    CompiledNuclei(i).FlaggingInfo = {};
    CompiledNuclei(i).ManualFrameApproved = [];
    CompiledNuclei(i).Approved = [];
end

schnitzcells = getSchnitzcells(liveExperiment);
%%
NoHistone = true;
for ch = 1:length(liveExperiment.Channels)
    if contains(lower(liveExperiment.Channels{ch}), 'his') 
        NoHistone = false;
    end
end

if NoHistone
    UseHistoneInfo = false;
end



his_zpadding = 1;
input_zpadding = 2;
for i=1:length(CompiledNuclei)%[test_idx]
    cp = CompiledNuclei(i);
    sc = schnitzcells(cp.schnitz);
    if cp.Flag6 == 1
        cp.FlaggingInfo.SickNucleus = true;
    else
        cp.FlaggingInfo.SickNucleus = false;
    end

    cp.Approved = 1;

    % Find True Nuclear start frame by inference or direct observation
    if ~isempty(cp.inferredAnaphaseFrame)
        if cp.inferredAnaphaseFrame == 0
            if ~isempty(cp.anaphaseFrame)
                if cp.anaphaseFrame >= nc_info(cp.nc-8)-1
                    cp.FlaggingInfo.TrueNucStart = max([cp.anaphaseFrame, 1]);
                    cp.FlaggingInfo.TNSinferred = false;
                else
                    cp.FlaggingInfo.TrueNucStart =  max([nc_info(cp.nc-8), 1]);
                    cp.FlaggingInfo.TNSinferred = true;
                end
            else
                cp.FlaggingInfo.TrueNucStart =  max([1, nc_info(cp.nc-8)]);
                cp.FlaggingInfo.TNSinferred = true;
            end
        elseif cp.anaphaseFrame >= nc_info(cp.nc-8)-1
            cp.FlaggingInfo.TrueNucStart = max([1, cp.anaphaseFrame]);
            cp.FlaggingInfo.TNSinferred = true;
        else
            cp.FlaggingInfo.TrueNucStart =  max([1,nc_info(cp.nc-8)]);
            cp.FlaggingInfo.TNSinferred = true;
        end
    elseif ~isempty(cp.anaphaseFrame)
        if cp.anaphaseFrame >= nc_info(cp.nc-8)-1
            cp.FlaggingInfo.TrueNucStart = max([1, cp.anaphaseFrame]);
            cp.FlaggingInfo.TNSinferred = false;
        else
            cp.FlaggingInfo.TrueNucStart =  max([1, nc_info(cp.nc-8)]);
            cp.FlaggingInfo.TNSinferred = true;
        end
    else
        cp.FlaggingInfo.TrueNucStart =  max([1, nc_info(cp.nc-8)]);
        cp.FlaggingInfo.TNSinferred = true;
    end
    % Find True Nuclear end frame by inference or direct observation
    if sc.containsLastFrameOfCycle
        cp.FlaggingInfo.TrueNucEnd = max([1, max(sc.frames)]);
        cp.FlaggingInfo.TNEinferred =false;
    else
        cp.FlaggingInfo.TrueNucEnd = max([1, nc_info(cp.nc-7)-1]);
        cp.FlaggingInfo.TNEinferred =true;
    end
    
    cp.FlaggingInfo.TrueFrames =...
        cp.FlaggingInfo.TrueNucStart:cp.FlaggingInfo.TrueNucEnd;
    
    cp.FlaggingInfo.FrameApproved = ones(1, length(cp.FlaggingInfo.TrueFrames), 'logical');
    cp.ManualFrameApproved = cp.FrameApproved;
    % Make Histone and spot Fluorescence traces and z positions
    NFrames = length(cp.FlaggingInfo.TrueFrames);
    if UseZInfo
        MaxFluoLevel = NaN(1, NFrames);
        MaxZPos= NaN(1, NFrames);
        InputZTraceApproved = zeros(1, NFrames);
        DetectedNucleus = zeros(1, NFrames);
        if UseHistoneInfo
            MaxHisFluoLevel = NaN(1, NFrames);
            MaxHisZPos= NaN(1, NFrames);
            HisTraceApproved = zeros(1, NFrames);
            

            for j = 1:length(cp.Frames)
                schnitz_index = find(cp.FlaggingInfo.TrueFrames == cp.Frames(j), 1);
               
                HisVector = cp.HistoneFluo(j,2:zDim+1);
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
            cp.FlaggingInfo.MaxHisFluoLevel = MaxHisFluoLevel;
            cp.FlaggingInfo.MaxHisZPos = MaxHisZPos;
            cp.FlaggingInfo.HisTraceApproved = HisTraceApproved;
        end
        
        
        for j = 1:length(cp.Frames)
            cp_index = find(cp.FlaggingInfo.TrueFrames == cp.Frames(j), 1);
            if isempty(cp_index) | cp.FrameApproved(j) == 0
                continue
            end
            DetectedNucleus(cp_index) = 1;
            MaxFluoLevel(cp_index) = cp.FluoTimeTrace(j);
            MaxZPos(cp_index) = cp.FluoZInfo(j);
            if (MaxZPos(cp_index) > input_zpadding) & (MaxZPos(cp_index) <= zDim - input_zpadding)
                InputZTraceApproved(cp_index) = 1;
            else
                InputZTraceApproved(cp_index) = -1;
            end
           
            
        end
        cp.FlaggingInfo.MaxInputFluoLevel = MaxFluoLevel;
        cp.FlaggingInfo.MaxInputZPos = MaxZPos;
        cp.FlaggingInfo.InputZTraceApproved = InputZTraceApproved;
        
        cp.FlaggingInfo.DetectedNucleus = DetectedNucleus;

        % Add spot position deviation velocity flags

    end
    
    FrameApproved = ones(1, length(cp.ManualFrameApproved), 'logical');
    FrameApprovedFinal = zeros(1, length(cp.FlaggingInfo.TrueFrames), 'logical');
    for j = 1:length(cp.Frames)
       tf_index = find(cp.FlaggingInfo.TrueFrames == cp.Frames(j), 1);
       if ~isempty(tf_index)
           FrameApproved(j) = cp.ManualFrameApproved(j);
           
           if UseZInfo 
               FrameApproved(j) = FrameApproved(j) & ...
               (cp.FlaggingInfo.InputZTraceApproved(tf_index) == 1);
           end
           if UseHistoneInfo
              FrameApproved(j) = FrameApproved(j) & ...
                 (cp.FlaggingInfo.HisTraceApproved(tf_index) == 1);
           end
               
       else
           FrameApproved(j) = false;
       end
       FrameApprovedFinal(tf_index) =  FrameApproved(j);
    end
    cp.FrameApproved = FrameApproved;
    cp.FlaggingInfo.FrameApprovedFinal = FrameApprovedFinal;
    
    CompiledNuclei(i) = cp;
end

%%

try
    save([liveExperiment.resultsFolder,filesep,'CompiledNuclei.mat'],...
        'CompiledNuclei','-v6');
catch
    save([liveExperiment.resultsFolder,filesep,'CompiledNuclei.mat'],...
        'CompiledNuclei','-v7.3', '-nocompression');
end



