function dataStatusWithDataType = findDataStatusTab(dataStatusFolders, dataType)

% dataStatusWithDataType = findDataStatusTab(dataStatusFolders, dataType)
%
% DESCRIPTION
% Finds the DataStatus.xlsx file that contains a tab which name matches the
% specified dataType variable
%
% INPUT
% dataStatusFolders: A cell array containing the folders of one or more
%                    DataStatus.xlsx files
% dataType: String specifying the name of the tab to find in DataStatus.xlsx
%
% OPTIONS
% N/A
%
% OUTPUT
% dataStatusWithDataType: Char array specifying the path to the folder
%                         that contains the DataStatus.xlsx file with the  
%                         specified dataType tab
%
%
% Author (contact): Meghan Turner (meghan_turner@berkeley.edu)
% Created: 5/15/2020
% Origin: Functionalized from code originally included in LoadMS2Sets.m,
%         written by Hernan Garcia (hggarcia@berkeley.edu)
% Last Updated: N/A

dataStatusWithDataType = '';
dataStatusToCheck = {};

for i = 1:length(dataStatusFolders)
    currDataStatusFolder = dataStatusFolders{i};
    currDataStatusDir = dir([currDataStatusFolder,filesep,'DataStatus.*']);
    
    dataStatusSheets = sheetnames([currDataStatusFolder,filesep,currDataStatusDir(1).name]);
    findSheets = strcmpi(dataStatusSheets,dataType);
    
    if sum(findSheets) == 1
        dataStatusToCheck{end + 1} = currDataStatusFolder;
    elseif sum(findSheets) > 1
        error(['More than one tab matching ', dataType, ' found. Please check ', currDataStatusFolder], '.')
    end
end

% Check that we found the tab name in one, and only one, DataStatus.xlsx
% file
if isempty(dataStatusToCheck)
    errormsg = ['No DataStatus.xlsx found with a tab named ', dataType,...
                '. \nDo you have the folder containing the DataStatus file listed as a ''Results'' or ''Dropbox'' folder in ComputerFolders.csv?'];
    error(sprintf(errormsg)) %need sprintf to get line break in error message
elseif length(dataStatusToCheck) > 1
    error(['More than one DataStatus.xlsx found with a tab named ', dataType, '.'])
else
    %Convert to char array during assignment for easier use of output
    dataStatusWithDataType = dataStatusToCheck{1};
end

