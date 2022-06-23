function CombineNuclearDorsalProfiles(SetPrefixes, SetLabel)
savedir = ['S:/Gabriella/Dropbox/FixedEmbryoProfiles/',filesep, SetLabel];
if ~isdir(savedir)
mkdir(savedir);
end
%%
for i = 1:length(SetPrefixes)
    Prefix = SetPrefixes{i};
    liveExperiment = LiveExperiment(Prefix);

    DorsalNuclearPath = [liveExperiment.resultsFolder, filesep, 'DorsalNuclearProfiles.mat'];
    load(DorsalNuclearPath, 'DorsalAvgAPProfiles', 'DorsalAvgNarrowAPProfiles', 'NarrowProfileNucleiFluoInfo',...
    'DorsalAPProfiles', 'DorsalNarrowAPProfiles','ProfileNucleiFluoInfo',...
    'ProfileNarrowNucleiFluoInfo','AllDorsalNucleiFluoInfo','DorsalStdAPProfiles',...
    'DorsalStdNarrowAPProfiles','DorsalCountAPProfiles','DorsalCountNarrowAPProfiles');
    if i == 1
        TempDorsalAvgAPProfiles = DorsalAvgAPProfiles;
        TempDorsalAvgNarrowAPProfiles = DorsalAvgNarrowAPProfiles;
        TempNarrowProfileNucleiFluoInfo = NarrowProfileNucleiFluoInfo;
        TempDorsalAPProfiles = DorsalAPProfiles;
        TempDorsalNarrowAPProfiles = DorsalNarrowAPProfiles;
        TempProfileNucleiFluoInfo = ProfileNucleiFluoInfo;
        TempProfileNarrowNucleiFluoInfo = ProfileNarrowNucleiFluoInfo;
        TempAllDorsalNucleiFluoInfo = AllDorsalNucleiFluoInfo;
        TempDorsalStdAPProfiles = DorsalStdAPProfiles;
        TempDorsalStdNarrowAPProfiles = DorsalStdNarrowAPProfiles;
        TempDorsalCountAPProfiles = DorsalCountAPProfiles;
        TempDorsalCountNarrowAPProfiles = DorsalCountNarrowAPProfiles;
        
        SetIDs = ones(1, size(DorsalAvgAPProfiles, 1));
    else
        TempDorsalAvgAPProfiles = [TempDorsalAvgAPProfiles ; DorsalAvgAPProfiles];
        TempDorsalAvgNarrowAPProfiles = [TempDorsalAvgNarrowAPProfiles ;DorsalAvgNarrowAPProfiles];
        TempNarrowProfileNucleiFluoInfo = [TempNarrowProfileNucleiFluoInfo ; NarrowProfileNucleiFluoInfo];
        TempDorsalAPProfiles = [TempDorsalAPProfiles; DorsalAPProfiles];
        TempDorsalNarrowAPProfiles = [TempDorsalNarrowAPProfiles; DorsalNarrowAPProfiles];
        TempProfileNucleiFluoInfo = [TempProfileNucleiFluoInfo; ProfileNucleiFluoInfo];
        TempProfileNarrowNucleiFluoInfo = [TempProfileNarrowNucleiFluoInfo  ProfileNarrowNucleiFluoInfo];
        TempAllDorsalNucleiFluoInfo = [TempAllDorsalNucleiFluoInfo; AllDorsalNucleiFluoInfo];
        TempDorsalStdAPProfiles = [TempDorsalStdAPProfiles; DorsalStdAPProfiles];
        TempDorsalStdNarrowAPProfiles = [TempDorsalStdNarrowAPProfiles; DorsalStdNarrowAPProfiles];
        TempDorsalCountAPProfiles = [TempDorsalCountAPProfiles; DorsalCountAPProfiles];
        TempDorsalCountNarrowAPProfiles = [TempDorsalCountNarrowAPProfiles; DorsalCountNarrowAPProfiles];
        SetIDs = [SetIDs  i*ones(1, size(DorsalAvgAPProfiles, 1))];
    end
end

%%
DorsalAvgAPProfiles = TempDorsalAvgAPProfiles;
DorsalAvgNarrowAPProfiles = TempDorsalAvgNarrowAPProfiles;
NarrowProfileNucleiFluoInfo = TempNarrowProfileNucleiFluoInfo;
DorsalAPProfiles = TempDorsalAPProfiles;
DorsalNarrowAPProfiles = TempDorsalNarrowAPProfiles;
ProfileNucleiFluoInfo = TempProfileNucleiFluoInfo;
ProfileNarrowNucleiFluoInfo = TempProfileNarrowNucleiFluoInfo;
AllDorsalNucleiFluoInfo = TempAllDorsalNucleiFluoInfo;
DorsalStdAPProfiles = TempDorsalStdAPProfiles;
DorsalStdNarrowAPProfiles = TempDorsalStdNarrowAPProfiles;
DorsalCountAPProfiles = TempDorsalCountAPProfiles;
DorsalCountNarrowAPProfiles = TempDorsalCountNarrowAPProfiles;
save([savedir,filesep, 'DorsalNuclearProfiles.mat'], 'DorsalAvgAPProfiles', 'DorsalAvgNarrowAPProfiles', 'NarrowProfileNucleiFluoInfo',...
    'DorsalAPProfiles', 'DorsalNarrowAPProfiles','ProfileNucleiFluoInfo',...
    'ProfileNarrowNucleiFluoInfo','AllDorsalNucleiFluoInfo','DorsalStdAPProfiles',...
    'DorsalStdNarrowAPProfiles','DorsalCountAPProfiles','DorsalCountNarrowAPProfiles', 'SetIDs');
