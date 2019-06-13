function [fnames, ncs_used, dataPath, movieDatabasePath, Prefixes] = processCompileTracesOptions(varargin)
%PROCESSCOMPILETRACESOPTIONS Summary of this function goes here
%   Detailed explanation goes here
    
    varargin = varargin{1};
    
    Prefixes = varargin{1};
    
    dataPathOpt = false;
    ncOpt = false;
    for  i = 2:numel(varargin)-1
        if ischar(varargin{i}) 
            eval([varargin{i} ' =  varargin{i+1};'])
        end
        dataPathOpt = strcmpi(varargin{i},'dataPath');
        ncOpt = strcmpi(varargin{i},'ncs_used');
    end
    if ~ncOpt 
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
        [~, ~, dataPathDefault, ~, ~, ~, movieDatabasePath] = DetermineLocalFolders(Prefix);
        slashes = strfind(movieDatabasePath,'\');
        movieDatabasePath = movieDatabasePath(1:slashes(end));
        fnames{pidx} = [dataPathDefault,filesep,Prefix];        
    end
    if ~dataPathOpt
        dataPath = dataPathDefault;
    end
end

