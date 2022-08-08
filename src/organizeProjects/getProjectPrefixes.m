function prefixes = getProjectPrefixes(dataTypes,varargin)

% function prefixes = getProjectPrefixes(dataTypes,varargin)
%
% DESCRIPTION
% Returns the Prefixes for all experiments in the project tab(s) of 
% DataStatus.xlsx specified by 'dataTypes'
% This is an even more barebones version of the 'justPrefixes' option in
% LoadMS2Set.
%
% INPUT
% dataType: This is a cell array of strings that match the name(s) of the 
%           tab(s) in DataStatus.xlsx that you wish to analyze. 
%           (Note: They don't all need to be in the same DataStatus file!)
%
% OPTIONS
% 'onlyApproved': Limits the prefixes returned by this function to only
%                 those that have been marked 'ready' or 'ApproveAll' in
%                 the CompileParticles row of the DataStatus tab
% 'onlyUnapproved': Limits the prefixes returned by this function to only
%                   those that have NOT been marked 'ready' or 'ApproveAll'
%                   in the CompileParticles row of the DataStatus tab
% 'customApproved': Limits the prefixes returned by this function to only
%                   those that have been marked with a custom approval
%                   flag.
%                   Follow with the name of the row in DataStatus that 
%                   contains the approval flag (which must be a boolean,
%                   1 for approved and 0 for unapproved)
%                   e.g. getProjectPrefixes(dataType, 'customApproved',
%                        'ReadyForEnrichment')
%
% OUTPUT
% prefixes: n x 1 cell array containing the Prefixes for the specified
%           project tab(s)
%
% Author (contact): Meghan Turner (meghan_turner@berkeley.edu)
% Created: 5/17/2020
% Last Updated: MT 2022-06-14, added functionality to get prefixes for
%                              multiple DataStatus tabs at once

%% Getting everything set up 

onlyApproved = false; 
onlyUnapproved = false;
customApproved = false;
customApprovalFlag = '';
prefixes = {}; %#ok<NASGU>
dataStatusFilename = 'DataStatus.*';    %This naming convention is now enforced inside findDataStatus.m

% If user passed a single dataType as a character array, convert to cell
% array to ensure downstream compatibility
if ischar(dataTypes)
    projectList = {dataTypes};
else
    projectList = dataTypes;
end

% Determine which type of approval flag we're dealing with
for i= 1:length(varargin)
    if strcmpi(varargin{i}, 'onlyApproved')
        onlyApproved = true; 
    elseif strcmpi(varargin{i}, 'onlyUnapproved')
        onlyUnapproved = true;
    elseif strcmpi(varargin{i}, 'customApproved')
        customApproved = true;
        customApprovalFlag = varargin{i+1};
    end
end

% Get all Dropbox/Results folders
[~, ~, ~, ~, ~, ~, ~, ~, ~, allDropboxFolders] =  DetermineLocalFolders;

% Find all DataStatus.xlsx files in all DropboxFolders
dataStatusFolders = findDataStatus(allDropboxFolders);


%% Loop over all dataTypes input by user

allExperimentNames = {};

for i = 1:numel(projectList)
    currDataType = projectList{i};
    
    % Look in all DataStatus.xlsx files to find the tab specified by the 
    % current dataType
    dataStatusWithDataTypeFolder = findDataStatusTab(dataStatusFolders, currDataType);

    % Redefine the DropboxFolder according to the DataStatus.xlsx we'll use
    dropboxFolder = dataStatusWithDataTypeFolder;

    % Load the contents of the DataStatus.XLSX tab we just found
    dataStatusDir = dir([dropboxFolder,filesep,dataStatusFilename]);
    dataTypeTabContents = readcell([dropboxFolder,filesep,dataStatusDir(1).name], ...
                                   'Sheet', currDataType);

    %Get the Prefixes for all experiments in this dataType
    [currPrefixes,~] = getPrefixesFromDataStatusTab(dataTypeTabContents);


    % Handle 'Approved' dataset options
    if onlyApproved || onlyUnapproved
        %Find which datasets are "ready" or "approved". All others are "ignored".
        compileRow = find(strcmpi(dataTypeTabContents(:,1),'AnalyzeLiveData Compile Particles') |...
                          strcmpi(dataTypeTabContents(:,1),'CompileParticles') |...
                          strcmpi(dataTypeTabContents(:,1),'CompileNuclearProtein'));

        readyLogicalArray = strcmpi(dataTypeTabContents(compileRow, 2:end),'READY') | ...
                       strcmpi(dataTypeTabContents(compileRow, 2:end),'ApproveAll');

        if onlyApproved
            outSets = find(readyLogicalArray);

            if isempty(outSets)
                warning('No ApproveAll or READY sets found.')
            end

        elseif onlyUnapproved
            outSets = find(~readyLogicalArray);
        end

        outPrefixes = currPrefixes(outSets);
        selected_prefixes = outPrefixes;

    elseif customApproved
        % Find the row in the DataStatus tab that contains the custom approval flag 
        approvalFlagRow = find(strcmpi(dataTypeTabContents(:,1),customApprovalFlag));
        % Contents of the approval flag row much be int booleans (1s or 0s) and
        % CANNOT be empty
        % If you get an error on this line, fill in any empty columns with 0s
        approvalFlagLogicalArray = cell2mat(dataTypeTabContents(approvalFlagRow,2:end));

        approvedPrefixes = currPrefixes(find(approvalFlagLogicalArray));
        selected_prefixes = approvedPrefixes;
    else
        selected_prefixes = currPrefixes;
    end
    
    % Add the prefixes from this dataTypes into the main, compiled list
    nNewPrefixes = length(selected_prefixes);
    allExperimentNames(end+1:end+nNewPrefixes, 1) = selected_prefixes;
end

% Not quite ready to move fully away from the prefix naming scheme
prefixes = allExperimentNames;