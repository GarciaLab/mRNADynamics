function [Prefix, SkipFrames, ProjectionType, PreferredFileNameForTest,...
    keepTifs, generateTifStacks, nuclearGUI,...
    skipExtraction, rootFolder, zslicesPadding,...
    dataType, skipNuclearProjection]...
    ...
    = exportDataForLivemRNA_processInputParameters(varargin)

%Look at parameters
shouldTrackNuclei = true;
SkipFrames = [];
Prefix = '';
%Default setting for z-projection is maximum-based.
%This may fail when high intensity reflections are present
ProjectionType = 'midsumprojection';
%Added new argument to specify a preferred file name and enable automatic testing
PreferredFileNameForTest = '';
keepTifs = false;
generateTifStacks = false;
nuclearGUI = true;
skipExtraction = false;
rootFolder = '';
zslicesPadding = false;
dataType = '';
exportMovieFiles = true;
ignoreCh3 = false;
skipNuclearProjection = false;

k=1;
while k<=length(varargin)
    if strcmpi(varargin{k},'skipframes')
        SkipFrames=varargin{k+1};
        k=k+1;
        warning('SkipFrame mode.')
    elseif strcmpi(varargin{k},'ProjectionType')
        ProjectionType = varargin{k+1};
    elseif isobject(varargin{k}) && isa(varargin{k}, 'PreferredFileForTest')
        PreferredFileForTest = varargin{k};
        PreferredFileNameForTest = PreferredFileForTest.fileName;
    elseif strcmpi(varargin{k}, 'keepTifs')
        keepTifs = true;
    elseif strcmpi(varargin{k}, 'generateTifs')
        generateTifStacks = true;
    elseif strcmpi(varargin{k}, 'dataType')
        dataType = varargin{k+1};
    elseif strcmpi(varargin{k}, 'nuclearGUI')
        nuclearGUI = varargin{k+1};
    elseif strcmpi(varargin{k}, 'ignoreCh3')
        ignoreCh3 = true;
    elseif strcmpi(varargin{k}, 'skipExtraction')
        skipExtraction = true;
    elseif strcmpi(varargin{k}, 'rootFolder')
        rootFolder = varargin{k+1};
    elseif strcmpi(varargin{k}, 'shouldTrackNuclei')
        shouldTrackNuclei = varargin{k+1};
    elseif strcmpi(varargin{k}, 'zslicesPadding')
        zslicesPadding = true;
    elseif strcmpi(varargin{k}, 'skipNuclearProjection')
        skipNuclearProjection = true;
    elseif strcmpi(varargin{k}, 'exportMovieFiles')
        exportMovieFiles = varargin{k+1};
    else
        %prefix can only go in first position
        if k == 1 && isempty(rootFolder)
            Prefix = varargin{k};
            disp(['ExportDataForLivemRNA using prefix: ', Prefix]);
        end
    end
    k=k+1;
end
end