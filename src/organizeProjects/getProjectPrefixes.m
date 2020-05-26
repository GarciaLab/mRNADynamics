function prefixes = getProjectPrefixes(dataType,varargin)

% function allPrefixes = getProjectPrefixes(dataType)
%
% DESCRIPTION
% Returns the Prefixes for all experiments in a project tab of 
% DataStatus.xlsx. This is an even more barebones version of the
% 'justPrefixes' option in LoadMS2Set.
%
% INPUT
% dataType: This is a string that is identical to the name of the tab in
%           DataStatus.xlsx that you wish to analyze.
%
% OPTIONS
% 'onlyApproved': Limits the prefixes returned by this function to only
%                 those that have been marked 'ready' or 'ApproveAll' in
%                 the CompileParticles row of the DataStatus tab
%
% OUTPUT
% prefixes: n x 1 cell array containing the Prefixes for this project tab
%
% Author (contact): Meghan Turner (meghan_turner@berkeley.edu)
% Created: 5/17/2020
% Last Updated: N/A

onlyApproved = false;
prefixes = {};
dataStatusFilename = 'DataStatus.*';    %This naming convention is now enforced inside findDataStatus.m

% Determine if 'onlyApproved' option applies
for i= 1:length(varargin)
    if strcmpi(varargin{i}, 'onlyApproved')
        onlyApproved = true;
    end
end

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

%Get the Prefixes for all datasets
[allPrefixes,~] = getPrefixesFromDataStatusTab(dataTypeTabContents);

%% 'onlyApproved' Option
if onlyApproved
    %Find which datasets are "ready" or "approved". All others are "ignored".
    compileRow = find(strcmpi(dataTypeTabContents(:,1),'AnalyzeLiveData Compile Particles') |...
                      strcmpi(dataTypeTabContents(:,1),'CompileParticles') |...
                      strcmpi(dataTypeTabContents(:,1),'CompileNuclearProtein'));

    readySets = find(strcmpi(dataTypeTabContents(compileRow,:),'READY') | ...
                   strcmpi(dataTypeTabContents(compileRow,:),'ApproveAll'));
    readySets = readySets - 1;  %shift to excluded the label column

    if isempty(readySets)
        error('No ApproveAll or READY sets found.')
    end
    
    readyPrefixes = allPrefixes(readySets);
    
    prefixes = readyPrefixes;
else
    prefixes = allPrefixes;
end