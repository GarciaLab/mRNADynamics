function Temp_sp =  pullTempSpFromDataStatus(dataType, Prefix)

dataStatusFilename = 'DataStatus.*';    %This naming convention is now enforced inside findDataStatus.m

% Get all Dropbox/Results folders
[~, ~, ~, ~, ~, ~, ~, ~, ~, allDropboxFolders] =  DetermineLocalFolders;

% Find all DataStatus.xlsx files in all DropboxFolders
dataStatusFolders = findDataStatus(allDropboxFolders);

%Look in all DataStatus.xlsx files to find the tab specified by the dataType
%user input
dataStatusWithDataTypeFolder = findDataStatusTab(dataStatusFolders, dataType);

%Redefine the DropboxFolder according to the DataStatus.xlsx we'll use
dropboxFolder = dataStatusWithDataTypeFolder;

%Load the contents of the DataStatus.XLSX tab we just found
dataStatusDir = dir([dropboxFolder,filesep,dataStatusFilename]);
dataTypeTabContents = readcell([dropboxFolder,filesep,dataStatusDir(1).name], ...
    'Sheet', dataType);
dataStatusRownames = {};
for i = 1:size(dataTypeTabContents, 1)
    if ~ismissing(dataTypeTabContents{i, 1})
        dataStatusRownames{i} = dataTypeTabContents{i, 1};
    else
        dataStatusRownames{i} = '';
    end
end

allPrefixes = getProjectPrefixes(dataType);
ColIndex = find(ismember(allPrefixes, Prefix))+1;

TempSPIndex = find(contains(dataStatusRownames, 'Temp_sp'));
if isnumeric(dataTypeTabContents{TempSPIndex, ColIndex})
    Temp_sp = double(dataTypeTabContents{TempSPIndex, ColIndex});
end
end