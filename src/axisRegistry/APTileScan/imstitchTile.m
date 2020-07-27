function imm2 = imstitchTile(tile_array, varargin)
% author: Gabriella Martini
% date created: 12/26/19
% date last modified: 7/27/20
% 
% Takes a tile_array input and uses its stitching configuration to position
% tiles in a stitched image. 
%
% output: imm2 is a stitched image. 
%% Parse User input 
Filter = false;
ZStack = false;
if ~isempty(varargin)
    x = 1;
    while x <= length(varargin)
        switch varargin{x}
            case{'Filter'}
                Filter = varargin{x+1};
                x = x+1;
            case{'ZStack'}
                ZStack = varargin{x+1};
                x = x+1;
        end
        x = x+1;
    end
end
%% Generate stitched gray-scale image
if ZStack % Creates stitched image using max-projected raw data
    tiles = tile_array.zstacks;
    numz = size(tile_array.zstacks{1}, 3);
elseif ~Filter
    tiles = tile_array.tiles;
    numz = 1;
else % Creates stitched image using gaussian filtered max projections of the raw data
    tiles = tile_array.imgs;
    numz = 1;
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
imm2 = zeros(numrows, numcols, numz, 'uint16');
imcounts = zeros(numrows, numcols, numz, 'uint16');
for n =1:NTiles
    tn = tiles{n};
    r = rows(n);
    c = cols(n);
    h = heights(n);
    w = widths(n);
    imcounts(r:(r+h-1),c:(c+w-1), :) = imcounts(r:(r+h-1),c:(c+w-1), :) +1;
    imm2(r:(r+h-1),c:(c+w-1), :) = imm2(r:(r+h-1),c:(c+w-1), :) +tn;
end
imm2 = imm2./imcounts;
imcounts = zeros(numrows, numcols, 'uint16');
for n =1:NTiles
    tn = tiles{n};
    r = rows(n);
    c = cols(n);
    h = heights(n);
    w = widths(n);
    gr = tile_array.grid_positions{n}(1);
    gc = tile_array.grid_positions{n}(2);
    rmin = r;
    cmin = c;
    rmax = r + h -1;
    cmax = c + w -1;
    rminprime = rmin;
    cminprime = cmin;
    rmaxprime = rmax;
    cmaxprime = cmax;
    for m=1:NTiles
        r2min = rows(m);
        r2max = rows(m) + heights(m) -1;
        c2min = cols(m);
        c2max = c2min + widths(m) -1;
        g2r = tile_array.grid_positions{m}(1);
        g2c = tile_array.grid_positions{m}(2);
        if (gr == g2r) & (gc == g2c +1)
            ovlp = c2max-cmin;
            cminprime = cmin + ceil(ovlp/2);
        elseif (gr == g2r) & (gc == g2c -1)
            ovlp = cmax-c2min;
            cmaxprime = cmax - ceil(ovlp/2);
        elseif (gr == g2r + 1) & (gc == g2c)
            ovlp = r2max-rmin;
            rminprime = rmin+ceil(ovlp/2);
        elseif (gr == g2r-1) & (gc == g2c)
            ovlp = rmax-r2min;
            rmaxprime = rmax -ceil(ovlp/2);
        end
    end
    amin = rminprime - rmin + 1;
    amax = rmaxprime - rmin + 1;
    bmin = cminprime - cmin + 1;
    bmax = cmaxprime - cmin + 1;
    tsub = tn(amin:amax, bmin:bmax,:);
    imcounts(rminprime:rmaxprime,cminprime:cmaxprime) = n;%imcounts(r:(r+h-1),c:(c+w-1)) +1;
    imm2(rminprime:rmaxprime,cminprime:cmaxprime,:) = tsub;
end
% Image pixel intensity is normalized by the number of tiles contributing
% to each pixel
%imm2 = imm2./imcounts;
end
