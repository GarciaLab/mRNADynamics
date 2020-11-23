% To be passed as arguments
ltp = LTP('HbJB3', 'input');
prefix_regex = '(?<year>\d+)-(?<month>\d+)-(?<day>\d+)-HbJB3-(?<temperature>[0-9_]+)C-(?<region>[A-Za-z]+)-Embryo(?<embryo_number>\d+)';
for i=1:length(ltp.ExperimentPrefixes)
    prefix_split = regexp(ltp.ExperimentPrefixes(i), prefix_regex, 'names');
    ltp.LegendLabels{i} = [prefix_split{1}.year,'-',prefix_split{1}.month,'-',...
        prefix_split{1}.day, ' ', strrep(prefix_split{1}.temperature, '_', '.'),'°C ',prefix_split{1}.region, ' Embryo ',...
        prefix_split{1}.embryo_number];
end

ExperimentsWithSegmentedCytoplasm = [];
for i=1:length(ltp.ExperimentPrefixes)
    if ltp.ExperimentStatuses{i}.hasSegmentedCytoplasm
        ExperimentsWithSegmentedCytoplasm = [ExperimentsWithSegmentedCytoplasm, i];
    end
end
[~,~,DropboxFolder,~,~]=...
    DetermineLocalFolders;
%% Load cytoplasm fluorescence information 
for i=2:length(ExperimentsWithSegmentedCytoplasm)
    EmbryoIndex = ExperimentsWithSegmentedCytoplasm(i);
    Prefix = ltp.ExperimentPrefixes{EmbryoIndex};
    segmentCytoplasm(Prefix);
    %lE = ltp.Experiments{EmbryoIndex};
    %load([DropboxFolder, filesep, Prefix, filesep, 'CytoplasmFluorescence.mat']);
end
%% 
UseImagesForCalculation = zeros(size(NuclearCount));
ncs = [lE.nc9, lE.nc10, lE.nc11, lE.nc12, lE.nc13, lE.nc14, size(NuclearCount, 1)];
for nc_idx = find(ncs(1:6) > 0);
NC = nc_idx + 8;
MaxCycleNuclearCount = max(max(NuclearCount(ncs(nc_idx):(ncs(nc_idx+1)-1),:)));

for frame=ncs(nc_idx):(ncs(nc_idx+1)-1)
    UseImagesForCalculation(frame,:) = NuclearCount(frame,:) >= .5*MaxCycleNuclearCount;
end
end









