function [noCompiledNuclei, justPrefixes, inputOutputFits, ...
          inputOutputModel, localMovieDatabase] ...
          = determineLoadMS2SetsOptions(varargin)
      
% [noCompiledNuclei, justPrefixes, inputOutputFits, ...
%          inputOutputModel, localMovieDatabase,dataStatusFolder] ...
%          = determineLoadMS2SetsOptions(varargin)
%
% DESCRIPTION
% Processes the user input options for LoadMS2Sets.m.
% Functionalized from code originally included in LoadMS2Sets.m, written by 
% Hernan Garcia (hggarcia@berkeley.edu).
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
%
%
% Author (contact): Meghan Turner (meghan_turner@berkeley.edu)
% Created: 5/13/2020
% Origin: Functionalized from code originally included in LoadMS2Sets.m,
%         written by Hernan Garcia (hggarcia@berkeley.edu)
% Last Updated: N/A
%


% Default options
noCompiledNuclei = false;
justPrefixes = false;
inputOutputFits = false;
inputOutputModel = '';
localMovieDatabase = false;


for i= 1:length(varargin)
    if strcmpi(varargin{i}, 'noCompiledNuclei')
        noCompiledNuclei = true;
    elseif strcmpi(varargin{i}, 'justPrefixes')
        justPrefixes = true;
    elseif strcmpi(varargin{i}, 'inputOutputFits')
        inputOutputFits = true;
        inputOutputModel = varargin{i+1};
    elseif strcmpi(varargin{i}, 'localMovieDatabase')
        localMovieDatabase = true;
    end
end

end