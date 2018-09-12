function [Prefix,folder] = getPrefixAndFolder
% getPrefixAndFolder
%
% DESCRIPTION:
% This function is to allow you to get the Prefix and Folder of the data of interest.
%
% ARGUMENTS
% None
%
% OPTIONS
% None
%
% Author (contact): Emma Luu (emma_luu@berkeley.edu)
% Last Updated: 1/25/2018
%
% Documented by: Emma Luu (emma_luu@berkeley.edu)

%Figure out the initial folders. We'll update the Drobpox one later on in the code.
[sourcePath,~,~,~, ~]=...
    DetermineLocalFolders;

folder=uigetdir(sourcePath,'Select folder with data');
slashPositions=strfind(folder,filesep);
Prefix=[folder((slashPositions(end-1)+1):(slashPositions(end)-1)),'-',...
    folder((slashPositions(end)+1):(end))];
end
