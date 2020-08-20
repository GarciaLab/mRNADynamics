function [Prefix,folder] = getPrefixAndFolder(varargin)
% getPrefixAndFolder
%
% DESCRIPTION
% This function is to allow you to get the Prefix and Folder of the data of interest.
%
% ARGUMENTS
% None
%
% OPTIONS
% path : Pass a string variable with the path that you want the function to
%        open to. The default folder it will open is your raw dynamics
%        folder. 
%        Ex. 
%           path = 'E:\';
%           prefix = getPrefixAndFolder(path);
%
% Author (contact): Emma Luu (emma_luu@berkeley.edu)
% Last Updated: 2/12/2019
%
% Documented by: Emma Luu (emma_luu@berkeley.edu)

%Figure out the initial folders. 
% The default is your raw dynamics folders
if ~isempty(varargin) % user input
    sourcePath = varargin{1};
else % Default
    [sourcePath,~,~,~, ~]=...
        DetermineLocalFolders;
end



folder=uigetdir(sourcePath,'Select folder with data');
slashPositions=strfind(folder,filesep);
Prefix=[folder((slashPositions(end-1)+1):(slashPositions(end)-1)),'-',...
    folder((slashPositions(end)+1):(end))];
end
