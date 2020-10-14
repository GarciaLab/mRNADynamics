function StitchFullEmbryoImages(Prefix, varargin)
% author: Gabriella Martini
% date created: 12/30/19
% date last modified: 8/23/20

%% Parse Inputs 
if ~exist('Prefix')
    FolderTemp=uigetdir(DropboxFolder,'Choose folder with files to analyze');
    Dashes=strfind(FolderTemp,filesep);
    Prefix=FolderTemp((Dashes(end)+1):end);
end 

FullyAutomate = false;
StitchManually = false;
keepExistingStitching = false;
manualStitchOrder = false;
selectRegions = false;
useSurfStitching = false;
manualSeeding = false;
x = 1;
while x <= length(varargin)
    switch varargin{x}
        case{'keepExistingStitching'}
            keepExistingStitching = true;
        case {'FullyAutomate'}
            FullyAutomate = true;
            fprintf('Stitching fully automated.\n')
        case {'StitchManually'}
            StitchManually = true;
            fprintf('Stitching to be performed manually.\n')
        case {'MaxDeltaR'}
            MaxDeltaR = varargin{x+1};
            x = x+1;
            fprintf('Max change in row overlap to be used in stitching loop: %d\n', MaxDeltaR)
        case{'MaxDeltaC'}
            MaxDeltaC = varargin{x+1};
            x = x+1;
            fprintf('Max change in column overlap to be used in stitching loop: %d\n', MaxDeltaC)
        case{'manualStitchOrder'}
            manualStitchOrder = true;
        case{'selectStitchingRegions'}
            selectRegions=true;
        case{'useSurfStitchingInfo'}
            useSurfStitching=true;
        case{'manualSeeding'}
            manualSeeding=true;
        otherwise
            error(['Flag "', varargin{1}{x},'" not valid'])
    end
    x = x +1;
end

%% Stitch Surface and MidSaggital plane full embryo images


if ~keepExistingStitching
    varargin2 = {};
    if FullyAutomate
        varargin2{length(varargin2) + 1} = 'FullyAutomate';
    end
    if StitchManually 
        varargin2{length(varargin2) + 1} = 'StitchManually';
    end
    if manualStitchOrder
        varargin2{length(varargin2) + 1} = 'manualStitchOrder';
    end
    if selectRegions
        varargin2{length(varargin2) + 1} = 'selectStitchingRegions';
    end
    if exist('MaxDeltaR', 'var')
        varargin2{length(varargin2) + 1} = 'MaxDeltaR';
        varargin2{length(varargin2) + 1} = MaxDeltaR;
    end
    if exist('MaxDeltaC', 'var')
        varargin2{length(varargin2) + 1} = 'MaxDeltaC';
        varargin2{length(varargin2) + 1} = MaxDeltaC;
    end
    if useSurfStitching
        varargin2{length(varargin2) + 1} = 'useSurfStitchingInfo';
    end
    if manualSeeding
        varargin2{length(varargin2) + 1} = 'manualSeeding';
    end
    
    if length(varargin2) > 0
        EmbryoTileStitch(Prefix, 'Surf', varargin2);
        EmbryoTileStitch(Prefix, 'Mid', varargin2);
    else
        EmbryoTileStitch(Prefix, 'Surf');
        EmbryoTileStitch(Prefix, 'Mid');
    end
    
end
%% 
AlignFullEmbryoImages(Prefix)
