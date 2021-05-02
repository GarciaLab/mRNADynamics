function ltmp = AddFluoCoeffs(ltmp, FluoParamPath)
 
ltmp.UniqueTemperatures = fliplr(unique(ltmp.Temp_sps));


if ~isempty(FluoParamPath)
    if isfile([FluoParamPath filesep 'FluoCoefficients.mat'])
        load([FluoParamPath filesep 'FluoCoefficients.mat']);
        ltmp.FluoScalingInfo = FluoInfo;
        FluoTemps =  FluoInfo.Temperatures;
        ltmp.FluoCoeffs = NaN(1, length(ltmp.UniqueTemperatures));
        for i = 1:length(FluoTemps)
            match_index = find(round(ltmp.UniqueTemperatures, 1) == round(FluoTemps(i), 1));
            if ~isempty(match_index)
                ltmp.FluoCoeffs(match_index) = FluoInfo.FluoCoeffs(i);
            end
        end  
    else
        ltmp.FluoCoeffs = NaN(1, length(ltmp.UniqueTemperatures));
        warning('No FluoCoefficients.mat file found in FluoParamPath');
    end
 
else
    ltmp.FluoCoeffs = NaN(1, length(ltmp.UniqueTemperatures));
    
end