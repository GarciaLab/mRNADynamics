function ltmp = AddMeanProfiles(ltmp)
    deltaT = ltmp.time_delta;
    MeanProfiles = {};
    disp('Adding Mean Profiles.')
    for i = 1:length(ltmp.ExperimentPrefixes)
        disp([num2str(i),'/', num2str(length(ltmp.ExperimentPrefixes))]);
        if ltmp.ExperimentStatuses{i}.hasCompiledParticles
            MeanProfiles{i} = CalculateMS2APProfiles(ltmp.ExperimentPrefixes{i}, deltaT);
        else
            MeanProfiles{i} = [];
        end
    end
    ltmp.MeanProfiles = MeanProfiles;
end