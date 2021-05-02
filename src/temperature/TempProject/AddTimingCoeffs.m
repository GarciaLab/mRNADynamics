function ltmp = AddTimingCoeffs(ltmp, TimingParamPath)
 
ltmp.UniqueTemperatures = fliplr(unique(ltmp.Temp_sps));



if ~isempty(TimingParamPath)
    if isfile([TimingParamPath filesep 'DevTimeCoefficients.mat'])
        load([TimingParamPath filesep 'DevTimeCoefficients.mat']);
        ltmp.TimeScalingInfo = TimingInfo;
        DevTimeTemps =  TimingInfo.Temperatures;
        ltmp.TimingCoeffs = NaN(1, length(ltmp.UniqueTemperatures));
        for i = 1:length(DevTimeTemps)
            match_index = find(round(ltmp.UniqueTemperatures, 1) == round(DevTimeTemps(i), 1));
            if ~isempty(match_index)
                ltmp.TimingCoeffs(match_index) = TimingInfo.TimeCoeffs(i);
            end
        end 
    else
        ltmp.TimingCoeffs= NaN(1, length(ltmp.UniqueTemperatures));
        warning('No DevTimeCoefficients.mat file found in TimingParamPath');
    end
 
else
    ltmp.TimingCoeffs = NaN(1, length(ltmp.UniqueTemperatures));
end

