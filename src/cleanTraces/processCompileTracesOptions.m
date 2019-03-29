function [fnames, ncs_used, DropboxFolder, Prefixes] = processCompileTracesOptions(varargin)
%PROCESSCOMPILETRACESOPTIONS Summary of this function goes here
%   Detailed explanation goes here
    
    varargin = varargin{1};
    
    Prefixes = varargin{1};
    ncs_used = varargin{2};
    
    % promps user to choose data set if none were handed in as arguments
    if isempty(Prefixes)
        FolderTemp=uigetdir(DefaultDropboxFolder,'Select folder with data to analyze');
        Dashes=strfind(FolderTemp,filesep);
        Prefixes= {FolderTemp((Dashes(end)+1):end)};
    end
    
    % extracts addresses for each Prefix
    fnames = cell(size(Prefixes));
    for pidx = 1:length(Prefixes)
        Prefix = Prefixes{pidx};
        [~,~,DropboxFolder,~, ~,...
            ~, ~, ~, ~, ~,~] = readMovieDatabase(Prefix);
        fnames{pidx} = [DropboxFolder,filesep,Prefix];
    end
    
end

