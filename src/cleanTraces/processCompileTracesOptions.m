function [fnames, ncs_used, dataPath, movieDatabasePath, Prefixes] = processCompileTracesOptions(varargin)
%PROCESSCOMPILETRACESOPTIONS Summary of this function goes here
%   Detailed explanation goes here
    
    varargin = varargin{1};
    
    Prefixes = varargin{1};
    if numel(varargin) >= 2
        ncs_used = varargin{2};
    else
        warning("desired nuclear cycle unspecified. Defaulting to nc14")
        ncs_used = 14;
    end
    
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
%         [~,~,DropboxFolder,~, ~,...
%             ~, ~, ~, ~, ~,~] = readMovieDatabase(Prefix);
        [~, ~, dataPath, ~, ~, ~, movieDatabasePath] = DetermineLocalFolders(Prefix);
        slashes = strfind(movieDatabasePath,'\');
        movieDatabasePath = movieDatabasePath(1:slashes(end));
        fnames{pidx} = [dataPath,filesep,Prefix];        
    end
    
end

