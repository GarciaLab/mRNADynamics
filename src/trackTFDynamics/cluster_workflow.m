% This is an example workflow for how you would use the scripts in the
% "trackTFDynamics" folder to segment Dorsal clusters and make a histogram
% of their pairwise distances to MS2 transcription spots

% Example dataset for Dorsal clusters
Prefix = '2019-11-26-2xDl_Venus_snaBAC_MCPmCherry_Leica_Zoom45_21uW14uW_01';

% (1) Correct for bleaching so you can use the same Weka classifier for the
%     whole movie (instead of making a new classifier for each time frame)
normalizeForBleaching(Prefix)

% (2) Mask out the cytoplasm so Weka only segments clusters that are inside
%     the nuclei
channelsToMask = [1];
maskCytoplasmForWeka(Prefix, 'includeChannels', channelsToMask, 'maskNormalizedImages')

% (3) Use Weka to segment clusters (and MS2 spots). You can find
%     instructions for this step in the LivemRNA user manual

% (4) Use the mRNADynamics pipeline script, segmentSpots.m, to segment
%      both the clusters and the MS2 transcription spots

% (5) Calculate the pairwise distances between the MS2 spot in a nucleus 
%     and all clusters detected in the same nucleus
ms2ClusterDistances = pairwiseDistMS2Clusters(Prefix);

% Repeat (1) - (5) as needed to analyze all target and control data

% (6) Generate plots of MS2-cluster distance distributions for target, 
%     control, and comparison datasets
%     Note: This is not a function, just a script. You must edit the
%           Prefixes variable at the beginning of the script.
plotMS2ClusterDistances
