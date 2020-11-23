function [NCDivisionInfo,DivisionStdErrorInfo, TemperatureInfo, FitParams, ActivationEnergies ] = ...
    CalculateDivisionCycleDurations(this, varargin)

% Parse inputs

UseSingleNuclei = false;
x = 1;

while x <= length(varargin)
    if strcmp(lower(varargin{x}), 'usesinglenuclei')
        UseSingleNuclei = true;
    end
    x = x+1;
end

% Define Constants and get DropboxFoldr
R = 8.314*10^(-3); % kJ * K^(-1)*mol^(-1)
[~,~,DropboxFolder,~,~]=...
    DetermineLocalFolders;

% Calculate Division cycle times
if ~UseSingleNuclei
    numSets = length(this.ProcessedExperiments);
    
    APResolution = this.Experiments{this.ProcessedExperiments(1)}.APResolution;
    numAPbins = uint16(1/APResolution+1);
    
    Temp_obs = this.Temp_obs(this.ProcessedExperiments);
    
    NCDivisionInfo = cell(1, 4);
    DivisionStdErrorInfo  = cell(1, 4);
    TemperatureInfo = cell(1, 4);
    FitParams = cell(1,4);
    alphas = NaN(1, 4);
    for i=1:numSets
        EmbryoIndex = this.ProcessedExperiments(i);
        Prefix = this.ExperimentPrefixes{EmbryoIndex};
        load([DropboxFolder,filesep,Prefix,filesep,'APDivision.mat']);
        FrameInfo = this.Experiments{EmbryoIndex}.getFrameInfo;
        EmbryoNCInfo = cell(1, 4);
        [ObservedNCs, ObservedAPbins] = find(APDivision(1:13,:) > 0);
        EmbryoIncludedNCs = min(ObservedNCs):max(ObservedNCs);
        for nc_idx=1:length(EmbryoIncludedNCs)
            NC = EmbryoIncludedNCs(nc_idx);
            for j=1:numAPbins
                if APDivision(NC,j)  == 0
                    continue
                end
                if ~isempty(EmbryoNCInfo{NC-9})
                    EmbryoNCInfo{NC-9}(length(EmbryoNCInfo{NC-9}) + 1) = (FrameInfo(APDivision(NC+1, j)).Time-FrameInfo(APDivision(NC, j)).Time)/60;
                else
                    EmbryoNCInfo{NC-9}= [(FrameInfo(APDivision(NC+1, j)).Time-FrameInfo(APDivision(NC, j)).Time)/60];
                end
                
            end
            if ~isempty(EmbryoNCInfo{NC-9})
                if ~isempty(NCDivisionInfo{NC-9})
                    NCDivisionInfo{NC-9}(length(NCDivisionInfo{NC-9}) + 1) = mean(EmbryoNCInfo{NC-9});
                    DivisionStdErrorInfo{NC-9}(length(DivisionStdErrorInfo{NC-9}) + 1) = std(EmbryoNCInfo{NC-9})/sqrt(length(EmbryoNCInfo{NC-9}));
                    TemperatureInfo{NC-9}(length(TemperatureInfo{NC-9}) + 1) = this.Temp_obs(EmbryoIndex);
                else
                    NCDivisionInfo{NC-9} = [mean(EmbryoNCInfo{NC-9})];
                    DivisionStdErrorInfo{NC-9}= [std(EmbryoNCInfo{NC-9})/sqrt(length(EmbryoNCInfo{NC-9}))];
                    TemperatureInfo{NC-9} = [this.Temp_obs(EmbryoIndex)];
                end
            end
            
        end
        
        
    end
    
else
    EmbryosWithNuclearTrackingChecked = [];
    for i = this.ProcessedExperiments
        if this.ExperimentStatuses{i}.hasCheckedNuclearTracking
            EmbryosWithNuclearTrackingChecked = [EmbryosWithNuclearTrackingChecked, i];
        end
    end
    numSets = length(EmbryosWithNuclearTrackingChecked);
    
    APResolution = this.Experiments{EmbryosWithNuclearTrackingChecked(1)}.APResolution;
    numAPbins = uint16(1/APResolution+1);
    
    Temp_obs = this.Temp_obs(EmbryosWithNuclearTrackingChecked);
    
    NCDivisionInfo = cell(1, 4);
    DivisionStdErrorInfo  = cell(1, 4);
    TemperatureInfo = cell(1, 4);
    FitParams = cell(1,4);
    alphas = NaN(1, 4);
    for i=1:numSets
        EmbryoIndex = EmbryosWithNuclearTrackingChecked(i);
        Prefix = this.ExperimentPrefixes{EmbryoIndex};
        EmbryoNCInfo = cell(1, 4);
        load([DropboxFolder,filesep,Prefix,filesep,'CompiledNuclei.mat']);
        FrameInfo = this.Experiments{EmbryoIndex}.getFrameInfo;
        xDim = this.Experiments{EmbryoIndex}.xDim;
        yDim = this.Experiments{EmbryoIndex}.yDim;
        for k = 1:length(CompiledNuclei)
            if (CompiledNuclei(k).Flag1 == 0) & (CompiledNuclei(k).Flag2 <  0.5) & ...
                    (CompiledNuclei(k).Flag3 <  0.5) & (CompiledNuclei(k).Flag4 <  0.5) & ...
                    (CompiledNuclei(k).Flag5 <  0.75) & (CompiledNuclei(k).Flag7 <  0.2) & ...
                    (CompiledNuclei(k).nc <  14)
                FinalXPos = CompiledNuclei(k).xPos(end);
                FinalYPos = CompiledNuclei(k).yPos(end);
                if min([abs(FinalXPos-xDim), abs(FinalYPos-yDim), xDim, yDim]) > 20
                    NucleusMinFrame = min(CompiledNuclei(k).Frames);
                    NucleusMaxFrame = max(CompiledNuclei(k).Frames);
                    NucleusCycleDuration = (FrameInfo(NucleusMaxFrame).Time-FrameInfo(NucleusMinFrame).Time)/60;
                    if ~isempty(EmbryoNCInfo{CompiledNuclei(k).nc-9})
                        EmbryoNCInfo{CompiledNuclei(k).nc-9}(length(EmbryoNCInfo{CompiledNuclei(k).nc-9}) + 1) = ...
                            NucleusCycleDuration;
                    else
                        EmbryoNCInfo{CompiledNuclei(k).nc-9} = [NucleusCycleDuration];
                    end
                end
                
            end
        end
        for NC=10:13
            if ~isempty(EmbryoNCInfo{NC-9})
                if ~isempty(NCDivisionInfo{NC-9})
                    NCDivisionInfo{NC-9}(length(NCDivisionInfo{NC-9}) + 1) = mean(EmbryoNCInfo{NC-9});
                    DivisionStdErrorInfo{NC-9}(length(DivisionStdErrorInfo{NC-9}) + 1) = std(EmbryoNCInfo{NC-9})/sqrt(length(EmbryoNCInfo{NC-9}));
                    TemperatureInfo{NC-9}(length(TemperatureInfo{NC-9}) + 1) = this.Temp_obs(EmbryoIndex);
                else
                    NCDivisionInfo{NC-9} = [mean(EmbryoNCInfo{NC-9})];
                    DivisionStdErrorInfo{NC-9}= [std(EmbryoNCInfo{NC-9})/sqrt(length(EmbryoNCInfo{NC-9}))];
                    TemperatureInfo{NC-9} = [this.Temp_obs(EmbryoIndex)];
                end
            end
        end
        
        
    end
    
end
for NC=10:13
    if ~isempty(NCDivisionInfo{NC-9})
        TemperatureFittingVector = -1./(TemperatureInfo{NC-9}+273);
        FitParams{NC-9} = polyfit(TemperatureFittingVector,log(NCDivisionInfo{NC-9}),1);
        alphas(NC-9) = FitParams{NC-9}(1);
    end
end
ActivationEnergies = alphas*R;
end

