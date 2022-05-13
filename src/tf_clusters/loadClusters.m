function Clusters = loadClusters(outputFolder)

%initialize Clusters. 
Clusters = cell(1,1);

% Replace if it already exists
if exist([outputFolder, filesep, 'Clusters.mat'], 'file')
    
    load([outputFolder, filesep, 'Clusters.mat'], 'Clusters')

%     % If there's only one channel, Particles, Spots and other structures 
%     % are not saved as cells. We turn them into a cell to ensure
%     % compatibility with downstream code.
%     if nCh==1 && ~iscell(Clusters)
%         Clusters = {Clusters};
%     end
end
