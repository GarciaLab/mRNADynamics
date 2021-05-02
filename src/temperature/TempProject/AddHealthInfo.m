function ltmp = AddHealthInfo(ltmp)
NumSets = length(ltmp.Experiments);




ltmp.EmbryoStats.SchnitzCount = NaN(NumSets, 6);
ltmp.EmbryoStats.FractionSickNuclei = NaN(NumSets, 6);
ltmp.EmbryoStats.FractionRejectedNuclei = NaN(NumSets, 6);
ltmp.EmbryoStats.FractionCompleteNuclei = NaN(NumSets, 6);
ltmp.EmbryoStats.FractionFirstLastFrameNuclei = NaN(NumSets, 6);
ltmp.EmbryoStats.MeanTotalXDistanceTraveled = NaN(NumSets, 6);
ltmp.EmbryoStats.MeanTotalYDistanceTraveled = NaN(NumSets, 6);
ltmp.EmbryoStats.MeanTotalDistanceTraveled = NaN(NumSets, 6);
ltmp.EmbryoStats.MeanDistanceTraveledPerSecond = NaN(NumSets, 6);
ltmp.EmbryoStats.MeanTotalXDisplacement = NaN(NumSets, 6);
ltmp.EmbryoStats.MeanTotalYDisplacement = NaN(NumSets, 6);
ltmp.EmbryoStats.MeanTotalDisplacement = NaN(NumSets, 6);
ltmp.EmbryoStats.MeanDisplacementPerSecond = NaN(NumSets, 6);
ltmp.EmbryoStats.StdTotalXDistanceTraveled = NaN(NumSets, 6);
ltmp.EmbryoStats.StdTotalYDistanceTraveled = NaN(NumSets, 6);
ltmp.EmbryoStats.StdTotalDistanceTraveled = NaN(NumSets, 6);
ltmp.EmbryoStats.StdDistanceTraveledPerSecond = NaN(NumSets, 6);
ltmp.EmbryoStats.StdTotalXDisplacement = NaN(NumSets, 6);
ltmp.EmbryoStats.StdTotalYDisplacement = NaN(NumSets, 6);
ltmp.EmbryoStats.StdTotalDisplacement = NaN(NumSets, 6);
ltmp.EmbryoStats.StdDisplacementPerSecond = NaN(NumSets, 6);
ltmp.EmbryoStats.NCDivisionInfo = NaN(NumSets, 6);
ltmp.EmbryoStats.DivisionStdInfo = NaN(NumSets, 6);





for SetIndex = 1:NumSets
    UseHealthSummaryInfo = false;
    if ~ismember(SetIndex, ltmp.ProcessedExperiments)
        continue
    end
    
    liveExperiment = ltmp.Experiments{SetIndex};
    HealthSummaryPath = [liveExperiment.resultsFolder, 'HealthSummary.mat'];
    if isfile(HealthSummaryPath)
        load(HealthSummaryPath);
        if HealthSummary.NuclearTrackingDone
            UseHealthSummaryInfo = true;
        end
    end
    
    if UseHealthSummaryInfo
        
        for NC = ltmp.IncludedNCs
            nc_idx = NC-9;
            ltmp.EmbryoStats.SchnitzCount(SetIndex, NC-8) = HealthSummary.SchnitzCount(nc_idx);
            ltmp.EmbryoStats.FractionSickNuclei(SetIndex, NC-8)  =HealthSummary.FractionSickNuclei(nc_idx);
            ltmp.EmbryoStats.FractionRejectedNuclei(SetIndex, NC-8)  = HealthSummary.FractionRejectedNuclei(nc_idx);
            ltmp.EmbryoStats.FractionCompleteNuclei(SetIndex, NC-8)  = HealthSummary.FractionCompleteNuclei(nc_idx);
            ltmp.EmbryoStats.FractionFirstLastFrameNuclei(SetIndex, NC-8)  = HealthSummary.FractionFirstLastFrameNuclei(nc_idx);
            ltmp.EmbryoStats.MeanTotalXDistanceTraveled(SetIndex, NC-8)   = HealthSummary.MeanTotalXDistanceTraveled(nc_idx);
            ltmp.EmbryoStats.MeanTotalYDistanceTraveled(SetIndex, NC-8)   = HealthSummary.MeanTotalYDistanceTraveled(nc_idx);
            ltmp.EmbryoStats.MeanTotalDistanceTraveled(SetIndex, NC-8)   = HealthSummary.MeanTotalDistanceTraveled(nc_idx);
            ltmp.EmbryoStats.MeanDistanceTraveledPerSecond(SetIndex, NC-8)   = HealthSummary.MeanDistanceTraveledPerSecond(nc_idx);
            ltmp.EmbryoStats.MeanTotalXDisplacement(SetIndex, NC-8)   = HealthSummary.MeanTotalXDisplacement(nc_idx);
            ltmp.EmbryoStats.MeanTotalYDisplacement(SetIndex, NC-8)   = HealthSummary.MeanTotalYDisplacement(nc_idx);
            ltmp.EmbryoStats.MeanTotalDisplacement(SetIndex, NC-8) = HealthSummary.MeanTotalDisplacement(nc_idx);
            ltmp.EmbryoStats.MeanDisplacementPerSecond(SetIndex, NC-8) = HealthSummary.MeanDisplacementPerSecond(nc_idx);
            ltmp.EmbryoStats.StdTotalXDistanceTraveled(SetIndex, NC-8) = HealthSummary.StdTotalXDistanceTraveled(nc_idx);
            ltmp.EmbryoStats.StdTotalYDistanceTraveled(SetIndex, NC-8) =HealthSummary.StdTotalYDistanceTraveled(nc_idx);
            ltmp.EmbryoStats.StdTotalDistanceTraveled(SetIndex, NC-8) = HealthSummary.StdTotalDistanceTraveled(nc_idx);
            ltmp.EmbryoStats.StdDistanceTraveledPerSecond(SetIndex, NC-8) =HealthSummary.StdDistanceTraveledPerSecond(nc_idx);
            ltmp.EmbryoStats.StdTotalXDisplacement(SetIndex, NC-8) = HealthSummary.StdTotalXDisplacement(nc_idx);
            ltmp.EmbryoStats.StdTotalYDisplacement(SetIndex, NC-8) = HealthSummary.StdTotalYDisplacement(nc_idx);
            ltmp.EmbryoStats.StdTotalDisplacement(SetIndex, NC-8) = HealthSummary.StdTotalDisplacement(nc_idx);
            ltmp.EmbryoStats.StdDisplacementPerSecond(SetIndex, NC-8) = HealthSummary.StdDisplacementPerSecond(nc_idx);
            if NC < 14
                ltmp.EmbryoStats.NCDivisionInfo(SetIndex, NC-8) = HealthSummary.NCDivisionInfo(nc_idx);
                ltmp.EmbryoStats.DivisionStdInfo(SetIndex, NC-8) = HealthSummary.DivisionStdInfo(nc_idx);
            end
        end
    end
end

end
