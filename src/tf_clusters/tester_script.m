results_folder = 'S:\Meghan\Dropbox\DorsalClustersResults';
prefix = '2022-03-29-Dl-mNeonGreen-MCP-mCh_snaBAC_settingsTest06_embryo04';
lineage_file = '_lin.mat';
load([results_folder, filesep, prefix, filesep, prefix, lineage_file])

max_trace = numel(schnitzcells(1).frames);
min_trace_length = floor(max_trace ./ 2);

for i = 1:numel(schnitzcells)
    if numel(schnitzcells(i).frames) < min_trace_length
        schnitzcells(i).Approved = 0;
    end
    
end

save([results_folder, filesep, prefix, filesep, prefix, lineage_file],'schnitzcells')