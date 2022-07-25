function AllCompiledEmbryos = NormalizeUnivScaledProfiles(AllCompiledEmbryos)


APbins = 0:0.025:1;
NumAPbins = length(APbins);
NChannels = 5;
AllSetInfo = GetFixedSetPrefixInfo;
StandardStainIndices = find(~AllSetInfo.Flipped);
for set_index = StandardStainIndices
    if ~isfield(AllCompiledEmbryos{set_index}, 'NormalizedUnivScaledProfiles')
        AllCompiledEmbryos{set_index}.NormalizedUnivScaledProfiles = {};
    end
    if ~isfield(AllCompiledEmbryos{set_index}, 'IndNormalizedUnivScaledProfiles')
        AllCompiledEmbryos{set_index}.NormalizedUnivScaledProfiles = {};
    end
    if ~isfield(AllCompiledEmbryos{set_index}, 'BgdSubNormalizedUnivScaledProfiles')
        AllCompiledEmbryos{set_index}.NormalizedUnivScaledProfiles = {};
    end
end
ProfileTypes = fieldnames(AllCompiledEmbryos{1}.BootstrappedUnivScaledProfiles);
NumMasterProfs = size(AllCompiledEmbryos{1}.BootstrappedUnivScaledProfiles.(ProfileTypes{1}).ControlScaling.TestSet.mean, 4);

for f_index = 1:length(ProfileTypes)
    for idx = 1:length(StandardStainIndices)
        set_index = StandardStainIndices(idx);
        AllCompiledEmbryos{set_index}.NormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet = {};
        AllCompiledEmbryos{set_index}.NormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.AllProfiles = NaN(size(AllCompiledEmbryos{set_index}.UnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling));
        AllCompiledEmbryos{set_index}.NormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.mean = NaN(size(AllCompiledEmbryos{set_index}.BootstrappedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.mean));
        AllCompiledEmbryos{set_index}.NormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.se = NaN(size(AllCompiledEmbryos{set_index}.BootstrappedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.mean));
        AllCompiledEmbryos{set_index}.NormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.counts  = AllCompiledEmbryos{set_index}.BootstrappedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.counts;
        AllCompiledEmbryos{set_index}.NormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.counts_above  = AllCompiledEmbryos{set_index}.BootstrappedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.counts_above;
        AllCompiledEmbryos{set_index}.NormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.counts_below  = AllCompiledEmbryos{set_index}.BootstrappedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.counts_below;
        AllCompiledEmbryos{set_index}.NC13.NormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling = {};
        AllCompiledEmbryos{set_index}.NC13.NormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.Test = NaN(size(AllCompiledEmbryos{set_index}.NC13.UnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.Test));
        
        AllCompiledEmbryos{set_index}.IndNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet = {};
        AllCompiledEmbryos{set_index}.IndNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.AllProfiles = NaN(size(AllCompiledEmbryos{set_index}.UnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling));
        AllCompiledEmbryos{set_index}.IndNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.mean = NaN(size(AllCompiledEmbryos{set_index}.BootstrappedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.mean));
        AllCompiledEmbryos{set_index}.IndNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.se = NaN(size(AllCompiledEmbryos{set_index}.BootstrappedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.mean));
        AllCompiledEmbryos{set_index}.IndNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.counts  = AllCompiledEmbryos{set_index}.BootstrappedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.counts;
        AllCompiledEmbryos{set_index}.IndNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.counts_above  = AllCompiledEmbryos{set_index}.BootstrappedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.counts_above;
        AllCompiledEmbryos{set_index}.IndNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.counts_below  = AllCompiledEmbryos{set_index}.BootstrappedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.counts_below;
        AllCompiledEmbryos{set_index}.NC13.IndNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling = {};
        AllCompiledEmbryos{set_index}.NC13.IndNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.Test = NaN(size(AllCompiledEmbryos{set_index}.NC13.UnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.Test));
        
        AllCompiledEmbryos{set_index}.BgdSubNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet = {};
        AllCompiledEmbryos{set_index}.BgdSubNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.AllProfiles = NaN(size(AllCompiledEmbryos{set_index}.UnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling));
        AllCompiledEmbryos{set_index}.BgdSubNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.mean = NaN(size(AllCompiledEmbryos{set_index}.BootstrappedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.mean));
        AllCompiledEmbryos{set_index}.BgdSubNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.se = NaN(size(AllCompiledEmbryos{set_index}.BootstrappedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.mean));
        AllCompiledEmbryos{set_index}.BgdSubNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.counts  = AllCompiledEmbryos{set_index}.BootstrappedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.counts;
        AllCompiledEmbryos{set_index}.BgdSubNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.counts_above  = AllCompiledEmbryos{set_index}.BootstrappedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.counts_above;
        AllCompiledEmbryos{set_index}.BgdSubNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.counts_below  = AllCompiledEmbryos{set_index}.BootstrappedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.counts_below;
        AllCompiledEmbryos{set_index}.NC13.BgdSubNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling = {};
        AllCompiledEmbryos{set_index}.NC13.BgdSubNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.Test = NaN(size(AllCompiledEmbryos{set_index}.NC13.UnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.Test));
    end
    for ch_index = 2:5
        for master_index = 1:NumMasterProfs
            ProfMaxes = zeros(1, length(StandardStainIndices));
            ProfMins = zeros(1, length(StandardStainIndices));
            for idx = 1:length(StandardStainIndices)
                set_index = StandardStainIndices(idx);
                SetChProf = AllCompiledEmbryos{set_index}.BootstrappedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.mean(:,:,ch_index,master_index);
                MovMeanProf = movmean(SetChProf, 5, 1, 'omitnan');
                ProfMaxes(idx) = max(MovMeanProf,[], 'all');
                ProfMins(idx) = min(MovMeanProf,[], 'all');
            end
            SharedMax = max(ProfMaxes);
            SharedMin = max(ProfMins);
            for idx = 1:length(StandardStainIndices)
                set_index = StandardStainIndices(idx);
                SetChProf = AllCompiledEmbryos{set_index}.BootstrappedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.mean(:,:,ch_index,master_index);
                AllCompiledEmbryos{set_index}.NormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.mean(:,:,ch_index, master_index) =...
                    (SetChProf-ProfMins(1))/(ProfMaxes(1)-ProfMins(1));
                SetChProfSE = AllCompiledEmbryos{set_index}.BootstrappedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.se(:,:,ch_index,master_index);
                AllCompiledEmbryos{set_index}.NormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.se(:,:,ch_index, master_index) =...
                    (SetChProfSE-ProfMins(1))/(ProfMaxes(1)-ProfMins(1));
                
                NC13Profs = AllCompiledEmbryos{set_index}.NC13.UnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.Test(:,:,ch_index, master_index);
                AllCompiledEmbryos{set_index}.NC13.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.Test(:,:,ch_index, master_index)= ...
                    (NC13Profs-ProfMins(1))/(ProfMaxes(1)-ProfMins(1));
                AllProfs = AllCompiledEmbryos{set_index}.UnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling(:,:,ch_index, master_index);
                AllCompiledEmbryos{set_index}.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(:,:,ch_index, master_index)= ...
                    (AllProfs-ProfMins(1))/(ProfMaxes(1)-ProfMins(1));
                
                AllCompiledEmbryos{set_index}.IndNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.mean(:,:,ch_index, master_index) =...
                    (SetChProf-ProfMins(idx))/(ProfMaxes(idx)-ProfMins(idx));
                SetChProfSE = AllCompiledEmbryos{set_index}.BootstrappedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.se(:,:,ch_index,master_index);
                AllCompiledEmbryos{set_index}.IndNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.se(:,:,ch_index, master_index) =...
                    (SetChProfSE-ProfMins(idx))/(ProfMaxes(idx)-ProfMins(idx));
                
                NC13Profs = AllCompiledEmbryos{set_index}.NC13.UnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.Test(:,:,ch_index, master_index);
                AllCompiledEmbryos{set_index}.NC13.IndNormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.Test(:,:,ch_index, master_index)= ...
                    (NC13Profs-ProfMins(idx))/(ProfMaxes(idx)-ProfMins(idx));
                AllProfs = AllCompiledEmbryos{set_index}.UnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling(:,:,ch_index, master_index);
                AllCompiledEmbryos{set_index}.IndNormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(:,:,ch_index, master_index)= ...
                    (AllProfs-ProfMins(idx))/(ProfMaxes(idx)-ProfMins(idx));
                
                
                AllCompiledEmbryos{set_index}.BgdSubNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.mean(:,:,ch_index, master_index) =...
                    (SetChProf-ProfMins(idx))/(ProfMaxes(1)-ProfMins(1));
                SetChProfSE = AllCompiledEmbryos{set_index}.BootstrappedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.se(:,:,ch_index,master_index);
                AllCompiledEmbryos{set_index}.BgdSubNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.se(:,:,ch_index, master_index) =...
                    (SetChProfSE-ProfMins(idx))/(ProfMaxes(1)-ProfMins(1));
                
                NC13Profs = AllCompiledEmbryos{set_index}.NC13.UnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.Test(:,:,ch_index, master_index);
                AllCompiledEmbryos{set_index}.NC13.BgdSubNormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.Test(:,:,ch_index, master_index)= ...
                    (NC13Profs-ProfMins(idx))/(ProfMaxes(1)-ProfMins(1));
                AllProfs = AllCompiledEmbryos{set_index}.UnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling(:,:,ch_index, master_index);
                AllCompiledEmbryos{set_index}.BgdSubNormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(:,:,ch_index, master_index)= ...
                    (AllProfs-ProfMins(idx))/(ProfMaxes(1)-ProfMins(1));
            end
        end
    end
end

%%
FlippedStainIndices = find(AllSetInfo.Flipped);
for set_index = FlippedStainIndices
    if ~isfield(AllCompiledEmbryos{set_index}, 'NormalizedUnivScaledProfiles')
        AllCompiledEmbryos{set_index}.NormalizedUnivScaledProfiles = {};
    end
    if ~isfield(AllCompiledEmbryos{set_index}, 'IndNormalizedUnivScaledProfiles')
        AllCompiledEmbryos{set_index}.NormalizedUnivScaledProfiles = {};
    end
    if ~isfield(AllCompiledEmbryos{set_index}, 'BgdSubNormalizedUnivScaledProfiles')
        AllCompiledEmbryos{set_index}.NormalizedUnivScaledProfiles = {};
    end
end
ProfileTypes = fieldnames(AllCompiledEmbryos{1}.BootstrappedUnivScaledProfiles);
NumMasterProfs = size(AllCompiledEmbryos{1}.BootstrappedUnivScaledProfiles.(ProfileTypes{1}).ControlScaling.TestSet.mean, 4);

for f_index = 1:length(ProfileTypes)
    for idx = 1:length(FlippedStainIndices)
        set_index = FlippedStainIndices(idx);
        AllCompiledEmbryos{set_index}.NormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet = {};
        AllCompiledEmbryos{set_index}.NormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.AllProfiles = NaN(size(AllCompiledEmbryos{set_index}.UnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling));
        AllCompiledEmbryos{set_index}.NormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.mean = NaN(size(AllCompiledEmbryos{set_index}.BootstrappedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.mean));
        AllCompiledEmbryos{set_index}.NormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.se = NaN(size(AllCompiledEmbryos{set_index}.BootstrappedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.mean));
        AllCompiledEmbryos{set_index}.NormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.counts  = AllCompiledEmbryos{set_index}.BootstrappedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.counts;
        AllCompiledEmbryos{set_index}.NormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.counts_above  = AllCompiledEmbryos{set_index}.BootstrappedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.counts_above;
        AllCompiledEmbryos{set_index}.NormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.counts_below  = AllCompiledEmbryos{set_index}.BootstrappedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.counts_below;
        AllCompiledEmbryos{set_index}.NC13.NormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling = {};
        AllCompiledEmbryos{set_index}.NC13.NormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.Test = NaN(size(AllCompiledEmbryos{set_index}.NC13.UnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.Test));
        
        AllCompiledEmbryos{set_index}.IndNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet = {};
        AllCompiledEmbryos{set_index}.IndNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.AllProfiles = NaN(size(AllCompiledEmbryos{set_index}.UnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling));
        AllCompiledEmbryos{set_index}.IndNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.mean = NaN(size(AllCompiledEmbryos{set_index}.BootstrappedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.mean));
        AllCompiledEmbryos{set_index}.IndNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.se = NaN(size(AllCompiledEmbryos{set_index}.BootstrappedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.mean));
        AllCompiledEmbryos{set_index}.IndNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.counts  = AllCompiledEmbryos{set_index}.BootstrappedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.counts;
        AllCompiledEmbryos{set_index}.IndNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.counts_above  = AllCompiledEmbryos{set_index}.BootstrappedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.counts_above;
        AllCompiledEmbryos{set_index}.IndNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.counts_below  = AllCompiledEmbryos{set_index}.BootstrappedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.counts_below;
        AllCompiledEmbryos{set_index}.NC13.IndNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling = {};
        AllCompiledEmbryos{set_index}.NC13.IndNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.Test = NaN(size(AllCompiledEmbryos{set_index}.NC13.UnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.Test));
        
        AllCompiledEmbryos{set_index}.BgdSubNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet = {};
        AllCompiledEmbryos{set_index}.BgdSubNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.AllProfiles = NaN(size(AllCompiledEmbryos{set_index}.UnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling));
        AllCompiledEmbryos{set_index}.BgdSubNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.mean = NaN(size(AllCompiledEmbryos{set_index}.BootstrappedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.mean));
        AllCompiledEmbryos{set_index}.BgdSubNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.se = NaN(size(AllCompiledEmbryos{set_index}.BootstrappedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.mean));
        AllCompiledEmbryos{set_index}.BgdSubNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.counts  = AllCompiledEmbryos{set_index}.BootstrappedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.counts;
        AllCompiledEmbryos{set_index}.BgdSubNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.counts_above  = AllCompiledEmbryos{set_index}.BootstrappedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.counts_above;
        AllCompiledEmbryos{set_index}.BgdSubNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.counts_below  = AllCompiledEmbryos{set_index}.BootstrappedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.counts_below;
        AllCompiledEmbryos{set_index}.NC13.BgdSubNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling = {};
        AllCompiledEmbryos{set_index}.NC13.BgdSubNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.Test = NaN(size(AllCompiledEmbryos{set_index}.NC13.UnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.Test));
    end
    for ch_index = 2:5
        for master_index = 1:NumMasterProfs
            ProfMaxes = zeros(1, length(FlippedStainIndices));
            ProfMins = zeros(1, length(FlippedStainIndices));
            for idx = 1:length(FlippedStainIndices)
                set_index = FlippedStainIndices(idx);
                SetChProf = AllCompiledEmbryos{set_index}.BootstrappedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.mean(:,:,ch_index,master_index);
                MovMeanProf = movmean(SetChProf, 5, 1, 'omitnan');
                ProfMaxes(idx) = max(MovMeanProf,[], 'all');
                ProfMins(idx) = min(MovMeanProf,[], 'all');
            end
            SharedMax = max(ProfMaxes);
            SharedMin = max(ProfMins);
            for idx = 1:length(FlippedStainIndices)
                set_index = FlippedStainIndices(idx);
                SetChProf = AllCompiledEmbryos{set_index}.BootstrappedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.mean(:,:,ch_index,master_index);
                AllCompiledEmbryos{set_index}.NormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.mean(:,:,ch_index, master_index) =...
                    (SetChProf-ProfMins(1))/(ProfMaxes(1)-ProfMins(1));
                SetChProfSE = AllCompiledEmbryos{set_index}.BootstrappedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.se(:,:,ch_index,master_index);
                AllCompiledEmbryos{set_index}.NormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.se(:,:,ch_index, master_index) =...
                    (SetChProfSE-ProfMins(1))/(ProfMaxes(1)-ProfMins(1));
                
                NC13Profs = AllCompiledEmbryos{set_index}.NC13.UnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.Test(:,:,ch_index, master_index);
                AllCompiledEmbryos{set_index}.NC13.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.Test(:,:,ch_index, master_index)= ...
                    (NC13Profs-ProfMins(1))/(ProfMaxes(1)-ProfMins(1));
                AllProfs = AllCompiledEmbryos{set_index}.UnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling(:,:,ch_index, master_index);
                AllCompiledEmbryos{set_index}.NormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(:,:,ch_index, master_index)= ...
                    (AllProfs-ProfMins(1))/(ProfMaxes(1)-ProfMins(1));
                
                AllCompiledEmbryos{set_index}.IndNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.mean(:,:,ch_index, master_index) =...
                    (SetChProf-ProfMins(idx))/(ProfMaxes(idx)-ProfMins(idx));
                SetChProfSE = AllCompiledEmbryos{set_index}.BootstrappedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.se(:,:,ch_index,master_index);
                AllCompiledEmbryos{set_index}.IndNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.se(:,:,ch_index, master_index) =...
                    (SetChProfSE-ProfMins(idx))/(ProfMaxes(idx)-ProfMins(idx));
                
                NC13Profs = AllCompiledEmbryos{set_index}.NC13.UnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.Test(:,:,ch_index, master_index);
                AllCompiledEmbryos{set_index}.NC13.IndNormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.Test(:,:,ch_index, master_index)= ...
                    (NC13Profs-ProfMins(idx))/(ProfMaxes(idx)-ProfMins(idx));
                AllProfs = AllCompiledEmbryos{set_index}.UnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling(:,:,ch_index, master_index);
                AllCompiledEmbryos{set_index}.IndNormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(:,:,ch_index, master_index)= ...
                    (AllProfs-ProfMins(idx))/(ProfMaxes(idx)-ProfMins(idx));
                
                
                AllCompiledEmbryos{set_index}.BgdSubNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.mean(:,:,ch_index, master_index) =...
                    (SetChProf-ProfMins(idx))/(ProfMaxes(1)-ProfMins(1));
                SetChProfSE = AllCompiledEmbryos{set_index}.BootstrappedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.se(:,:,ch_index,master_index);
                AllCompiledEmbryos{set_index}.BgdSubNormalizedUnivScaledProfiles.(ProfileTypes{f_index}).ControlScaling.TestSet.se(:,:,ch_index, master_index) =...
                    (SetChProfSE-ProfMins(idx))/(ProfMaxes(1)-ProfMins(1));
                
                NC13Profs = AllCompiledEmbryos{set_index}.NC13.UnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.Test(:,:,ch_index, master_index);
                AllCompiledEmbryos{set_index}.NC13.BgdSubNormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.Test(:,:,ch_index, master_index)= ...
                    (NC13Profs-ProfMins(idx))/(ProfMaxes(1)-ProfMins(1));
                AllProfs = AllCompiledEmbryos{set_index}.UnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling(:,:,ch_index, master_index);
                AllCompiledEmbryos{set_index}.BgdSubNormalizedUnivScaledProfiles.SlideRescaledDorsalAvgAPProfiles.ControlScaling.AllProfiles(:,:,ch_index, master_index)= ...
                    (AllProfs-ProfMins(idx))/(ProfMaxes(1)-ProfMins(1));
                
               
            end
        end
    end
end
