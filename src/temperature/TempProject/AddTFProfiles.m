function ltp = AddTFProfiles(ltp)
    ltp.TFProfiles = cell(1, length(ltp.ExperimentPrefixes));
    disp('Adding Mean Profiles.')
    for i = 1:length(ltp.ExperimentPrefixes)
        disp(num2str(i))
        if ~ismember(i, ltp.ProcessedExperiments)
            continue
        end
        if ltp.ExperimentStatuses{i}.hasCompiledNuclearProtein
            ltp.TFProfiles{i} = CalculateNBAPProfiles(ltp.ExperimentPrefixes{i}, ltp.time_delta);
        end
    end
end