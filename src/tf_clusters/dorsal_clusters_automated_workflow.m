close all
clear

dynamics_results_path = 'S:\Meghan\Dropbox\DorsalClustersResults';
raw_data_path = 'P:\DorsalClusters\Data\RawDynamicsData';
date = {'2021-11-02', '2021-11-08', '2021-11-10'};
[~, prefixes_all] = dir_to_moviedatabase_entries(raw_data_path, date);
data_type = {'2xDl-Ven-clusters_twiPE','2xDl-Ven-clusters_snaBAC','2xDl-Ven-clusters_hbBAC'};
% data_type = {'2xDl-Ven-clusters_twiPE', '2xDl-Ven-clusters_snaBAC', '2xDl-Ven-clusters_hbBAC'};

exp_filter_sick = find(~contains(prefixes_all,'sick') .* ~contains(prefixes_all,'noSpots'));
prefixes_filtered = prefixes_all(exp_filter_sick);

for i = 3:numel(prefixes_filtered) 
    disp(['processing prefix ' prefixes_filtered{i}])
% %     ExportDataForLivemRNA(prefixes_filtered{i},'nuclearGUI',false)
% %     TrackNuclei(prefixes_filtered{i}, 'noTrack')
% % 
% %     % (1) Correct for bleaching so you can use the same Weka classifier for the
% %     %     whole movie (instead of making a new classifier for each time frame)
% %     normalizeForBleaching(prefixes_filtered{i})
% %     % (2) Mask out the cytoplasm so Weka only segments clusters that are inside
% %     %     the nuclei
% %     channelsToMask = [1];
% %     maskCytoplasmForWeka(prefixes_filtered{i}, 'includeChannels', channelsToMask, 'maskNormalizedImages')
% % (7) Use the mRNADynamics pipeline script, segmentSpots.m, to segment
% %      both the clusters and the MS2 transcription spots
%     segmentSpots(prefixes_filtered{i}, [5000,10015],'skipChannel',1,'fit3D');
%     segmentSpots(prefixes_filtered{i}, [5000,10015],'skipChannel',2,'fit3D');
%     try
        TrackNuclei(prefixes_filtered{i})
        segmentSpots(prefixes_filtered{i}, [5500,5500]);
        TrackmRNADynamics(prefixes_filtered{i},'displayFigures');
        CompileParticles(prefixes_filtered{i});
        pairwiseDistMS2Clusters(prefixes_filtered{i});
%         disp(['success!:' prefixes_filtered{i}])
%     catch
%         disp(['errored out:' prefixes_filtered{i}])
%     end
end

% [comp_settings, raw_settings] = compareExperimentSettings(data_type);

% Because I have to run the two channels through segmentSpots separately
% for Dl cluster analysis, I have to combine them at some point into a
% normal Spots.mat that the code will recognize

% for i = [8 9 10 11 12 13] %1:numel(Prefixes)
%     disp(['processing Spots.mat for prefix ' Prefixes{i}])
% %     load(['S:\Meghan\Dropbox\DorsalClustersResults\' Prefixes{i} '\Spots.mat'])
%     load(['S:\Meghan\Dropbox\DorsalClustersResults\' prefixes_filtered{i} '\Spots_ch01.mat'])
%     Spots_ch1 = Spots;
%     clear Spots
%     load(['S:\Meghan\Dropbox\DorsalClustersResults\' prefixes_filtered{i} '\Spots_ch02.mat'])
%     Spots_ch2 = Spots;
%     clear Spots
%     Spots{1,1} = Spots_ch1{1,1};
%     Spots{1,2} = Spots_ch2{1,2};
%     save(['S:\Meghan\Dropbox\DorsalClustersResults\' prefixes_filtered{i} '\Spots.mat'],'Spots')
%     clear Spots Spots_ch1 Spots_ch2
% end