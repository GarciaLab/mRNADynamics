function ProcessedEnrichmentFolder = getProcessedEnrichmentFolder()
        
% ProcessedEnrichmentFolder = getProcessedEnrichmentFolder()
%
% DESCRIPTION
% Extracts folder paths from the user's ComputerFolders.csv file.
%
% PARAMETERS
% None
%
% OPTIONS
% None
%
% OUTPUT
% ProcessedEnrichmentFolder 
%
%
% Author (contact): Gabriella Martini (martini@berkeley.edu)
% Created: 3/17/21


CONFIG_CSV_PATH =  'ComputerFolders.csv';

computerFoldersValues = csv2cell(CONFIG_CSV_PATH, 'fromfile');

defaultDropboxFolder = getConfigValue(computerFoldersValues, 'DropboxFolder');
% Find all the other results folders
ProcessedEnrichmentFolderRow = contains(computerFoldersValues(:,1),'ProcessedEnrichmentFolder');

ProcessedEnrichmentFolder = computerFoldersValues(ProcessedEnrichmentFolderRow,2);
ProcessedEnrichmentFolder = [ProcessedEnrichmentFolder{1}, filesep];
end
