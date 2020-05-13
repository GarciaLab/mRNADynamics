function dataStatusFilenames = findDataStatus(DropboxFolders)

% dataStatusFilenames = findDataStatus(DropboxFolder)
%
% DESCRIPTION
% Locates the DataStatus.csv file(s) inside the specified Dropbox (results)
% folder(s).
%
% PARAMETERS
% varagin: user input options
%
% OPTIONS
% N/A
%
% OUTPUT
% noCompiledNuclei
% justPrefixes
% inputOutputFits
% inputOutputModel
% localMovieDatabase
% dataStatusFolder
%
%
% Author (contact): Meghan Turner (meghan_turner@berkeley.edu)
% Created: 5/13/2020
% Origin: Functionalized from code originally included in LoadMS2Sets.m,
%         written by Hernan Garcia (hggarcia@berkeley.edu)
% Last Updated: N/A

dataStatusFilenames = {};

dataStatusDir = dir([DropboxFolders,filesep,'Data*Status.*']);   %May or may not have space between 'Data' and 'Status'

if length(dataStatusDir) == 1
    dataStatusFilenames = [DropboxFolders,filesep,dataStatusDir];
elseif length(dataStatusDir) > 1
    error(['More than one DataStatus.XLS found in folder ',DropboxFolders])
end

end