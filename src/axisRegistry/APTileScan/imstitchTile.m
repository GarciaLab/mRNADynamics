function imm2 = imstitchTile(tile_array, varargin)
% author: Gabriella Martini
% date created: 12/26/19
% date last modified: 12/31/19
% 
% Takes a tile_array input and uses its stitching configuration to position
% tiles in a stitched image. 
%
% output: imm2 is a stitched image. 
%% Parse User input 
NoFilter = false;
if ~isempty(varargin)
    x = 1;
    while x <= length(varargin)
        switch varargin{x}
            case{'NoFilter'}
                NoFilter = varargin{x+1};
                x = x+1;
        end
        x = x+1;
    end
end
%% Generate stitched gray-scale image
if NoFilter % Creates stitched image using max-projected raw data
    tiles = tile_array.tiles;
else % Creates stitched image using gaussian filtered max projections of the raw data
    tiles = tile_array.imgs;
end

% Get stitching information from tile_array
NTiles = length(tiles);
heights = []; % tile heights
widths = []; % tile widths
for n=1:NTiles
    heights(n) = size(tiles{n}, 1);
    widths(n) = size(tiles{n}, 2);
end

rows = [tile_array.rows{:}];
cols = [tile_array.cols{:}];
rbounds = rows + heights-1;
cbounds = cols + widths-1;
numrows = max(rbounds);
numcols = max(cbounds);
imm2 = zeros(numrows, numcols);
imcounts = zeros(numrows, numcols);
for n =1:NTiles
    r = tile_array.rows{n};
    c = tile_array.cols{n};
    [h, w] = size(tiles{n});
    imcounts(r:(r+h-1),c:(c+w-1)) = imcounts(r:(r+h-1),c:(c+w-1)) +1;
    imm2(r:(r+h-1), c:(c+w-1)) =...
        imm2(r:(r+h-1), c:(c+w-1))+tiles{n};
end
% Image pixel intensity is normalized by the number of tiles contributing
% to each pixel
imm2 = imm2./imcounts;
end
