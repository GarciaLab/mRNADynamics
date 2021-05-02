function ltmp = AddAutocorrElongationTimes(ltmp, TimingParamPath)
 
ltmp.UniqueTemperatures = fliplr(unique(ltmp.Temp_sps));



if ~isempty(TimingParamPath)
    if isfile([TimingParamPath filesep 'autocorrElongationTimes.mat'])
        load([TimingParamPath filesep 'autocorrElongationTimes.mat']);
        ltmp.AutocorrElongationTimeInfo = ElongationTimeInfo;
        ElongTemps =  ElongationTimeInfo.Temperatures;
        ltmp.AutocorrElongationTimes = NaN(1, length(ltmp.UniqueTemperatures));
        for i = 1:length(ElongTemps)
            match_index = find(round(ltmp.UniqueTemperatures, 1) == round(ElongTemps(i), 1));
            if ~isempty(match_index)
                ltmp.AutocorrElongationTimes(match_index) = ElongationTimeInfo.ElongationTimes(i);
            end
        end 
    else
        ltmp.AutocorrElongationTimes= NaN(1, length(ltmp.UniqueTemperatures));
        warning('No autocorrElongationTimes.mat file found in TimingParamPath');
    end
 
else
    ltmp.AutocorrElongationTimes = NaN(1, length(ltmp.UniqueTemperatures));
end

