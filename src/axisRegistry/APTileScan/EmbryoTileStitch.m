function EmbryoTileStitch(Prefix, ID, varargin)
% author: Gabriella Martini
% date: 12/30/19
% last modified: 7/26/20
% Function takes full embryo images taken as tile scans on the Leica SP8
% and stitches them into tifs to be used for defining the AP axis
% and indentifying AP position. 
%% 

% Parse inputs

if ~exist('ID', 'var')
    prompt = ['Enter an "ID" for stitching.',...
    ' The standard inputs "Mid" and "Surf" will stitch the',...
    ' Midsagittal and Surface Full Embryo images using the',...
    ' "MidTile.lif" and "SurfTile.lif" files respectively.'];
    ID = input(prompt,'s');

end
ID = [upper(ID(1)), ID(2:end)];


FullyAutomate = false;
StitchManually = false;
selectRegions = false;
manualStitchOrder = false;
if ~isempty(varargin)
    x = 1;
    while x <= length(varargin{1})
        switch varargin{1}{x}
            case{'FullyAutomate'}
                FullyAutomate = true;
                fprintf('Stitching fully automated.\n')
            case {'StitchManually'}
                StitchManually = true;
                fprintf('Stitching to be performed manually.\n')
            case {'MaxDeltaR'}
                MaxDeltaR = varargin{1}{x+1};
                x = x+1;
                fprintf('Max change in row overlap to be used in stitching loop: %d\n', MaxDeltaR)
            case{'MaxDeltaC'}
                MaxDeltaC = varargin{1}{x+1};
                x = x+1;
                fprintf('Max change in column overlap to be used in stitching loop: %d\n', MaxDeltaC)
            case{'manualStitchOrder'}
                manualStitchOrder = true;
            case{'selectStitchingRegions'}
                selectRegions=true;
        end
        x = x +1;
    end
end



if StitchManually && FullyAutomate
    error(['EmbryoTileStitch cannot be run with both "StitchManually"',...
        'and "FullyAutomate" being true.']);
end

varargin2 = {};
if manualStitchOrder
    varargin2{length(varargin2) + 1} = 'manualStitchOrder';
end
if selectRegions
    varargin2{length(varargin2) + 1} = 'selectStitchingRegions';
end


if ~exist('ID', 'var')
    prompt = ['Enter an "ID" for stitching.',...
        ' The standard inputs "Mid" and "Surf" will stitch the',...
        ' Midsagittal and Surface Full Embryo images using the',...
        ' "MidTile.lif" and "SurfTile.lif" files respectively.'];
    ID = input(prompt,'s');
    if length(ID) > 1
        ID = [upper(ID(1)), ID(2:end)];
    else
        ID = upper(ID);
    end
end
if ~exist('MaxDeltaR', 'var') && ~StitchManually
    disp('Using default value of MaxDeltaR=100.')
    MaxDeltaR = 100;
end
if ~exist('MaxDeltaC', 'var') && ~StitchManually
    disp('Using default value of MaxDeltaC=100.')
    MaxDeltaC = 100;
end


%% Get the relevant folders for this data set

[RawDynamicsPath, ~, DefaultDropboxFolder, DropboxFolder, MS2CodePath,...
    PreProcPath,configValues,movieDatabasePath, movieDatabase] =...
    DetermineAllLocalFolders(Prefix);


% Parse Date and EmbryoName from Prefix
Dashes=findstr(Prefix,'-');
Date=Prefix(1:Dashes(3)-1);
EmbryoName=Prefix(Dashes(3)+1:end);


%% Generate Stitching Information

% Generate Initial Tile Array from Metadata files and raw data
% NewTileArrayFromMetadata(Prefix, ID);
% 
% % Get user input to improve initial seed for tile stitching
% if ~FullyAutomate
%     ManualStitchingCorrection(Prefix, ID);
% end

% Finalize Tile Array Stitching Positions using FindStitchingPositions
if ~StitchManually
    if length(varargin2) > 0
        tile_array = OptimizeStitching(Prefix, ID, MaxDeltaR, MaxDeltaC, varargin2);
    else
        tile_array = OptimizeStitching(Prefix, ID, MaxDeltaR, MaxDeltaC);
    end
    imagesc(imstitchTile(tile_array))       
end

end