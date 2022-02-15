function EmbryoTileStitch(Prefix, ID, varargin)
% author: Gabriella Martini
% date: 12/30/19
% last modified: 12/31/19
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

if ~isempty(varargin)
    x = 1;
    while x <= length(varargin{1})
        switch varargin{1}{x}
            case{'NIterations'}
                NIterations = varargin{1}{x+1};
                x = x + 1;
                fprintf('Number of Iterations: %d\n', NIterations);
            case{'FullyAutomate'}
                FullyAutomate = true;
                fprintf('Stiarginching fully automated.\n')
            case {'StitchManually'}
                StitchManually = true;
                fprintf('Stitching to be performed manually.\n')
            case {'MaxStep'}
                MaxStep = varargin{1}{x+1};
                x = x+1;
                fprintf('Max Step Size to be used in stitching loop: %d\n', MaxStep)
            case{'MaxOverlap'}
                MaxOverlap = varargin{1}{x+1};
                x = x+1;
                fprintf('Max overlap between adjacent tiles to be used in stitching loop: %d\n', MaxOverlap)
        end
        x = x +1;
    end
end

if StitchManually && FullyAutomate
    error(['EmbryoTileStitch cannot be run with both "StitchManually"',...
        'and "FullyAutomate" being true.']);
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
if ~exist('NIterations', 'var') && ~StitchManually
    disp('Using default value of NIterations=100.')
    NIterations = 100;
end
if ~exist('MaxStep', 'var') && ~StitchManually
    disp('Using default value of MaxStep=2.')
    MaxStep = 2;
end
if ~exist('MaxOverlap', 'var') && ~StitchManually
    disp('Using default value of MaxOverlap=200.')
    MaxOverlap = 500;
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
NewTileArrayFromMetadata(Prefix, ID);

% Get user input to improve initial seed for tile stitching
if ~FullyAutomate
    ManualStitchingCorrection(Prefix, ID);
end

% Finalize Tile Array Stitching Positions using FindStitchingPositions
if ~StitchManually

    tile_array = FindStitchingPositions(Prefix, ID, MaxStep, MaxOverlap, NIterations);
    close all
    imshow(imstitchTile(tile_array))
    if ~FullyAutomate
        prompt = ['Do you want to adjust the existing tiling and try stitching again (y/n)?'];
        keepFitting = input(prompt,'s');
        while keepFitting=='y'
            ManualStitchingCorrection(Prefix, ID);
            tile_array = FindStitchingPositions(Prefix, ID, MaxStep, MaxOverlap, NIterations);
            close all
            imshow(imstitchTile(tile_array))
            prompt = ['Do you want to adjust the existing tiling and try stitching again (y/n)?'];
            keepFitting = input(prompt,'s');
        end
    end
            
end

end