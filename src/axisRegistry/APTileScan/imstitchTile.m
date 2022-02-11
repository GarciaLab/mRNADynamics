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
try
    good_idx = [];
    for i = 1:NTiles
        if tile_array.use_tiles{i}
            good_idx(length(good_idx) + 1) = i;
        end 
    end
catch
    disp(['Old tile_array object does not have option to ignore tiles']);
    good_idx = 1:NTiles;
end
rows = [tile_array.rows{good_idx}];
cols = [tile_array.cols{good_idx}];
heights = heights(good_idx);
widths = widths(good_idx);
rbounds = rows + heights-1;
cbounds = cols + widths-1;

numrows = max(rbounds);
numcols = max(cbounds);
imm2 = zeros(numrows, numcols, numz, 'uint16');
imcounts = zeros(numrows, numcols, numz, 'uint16');
tiles_idx  = 1:NTiles;
for n =tiles_idx(good_idx)
    tn = tiles{n};
    r = tile_array.rows{n};
    c = tile_array.cols{n};
    h = size(tn, 1);
    w = size(tn, 2);
    imcounts(r:(r+h-1),c:(c+w-1), :) = imcounts(r:(r+h-1),c:(c+w-1), :) +1;
    imm2(r:(r+h-1),c:(c+w-1), :) = imm2(r:(r+h-1),c:(c+w-1), :) +tn;
end
imm2 = imm2./imcounts;
imcounts = zeros(numrows, numcols, 'uint16');
for n =tiles_idx(good_idx)
    try
        if ~tile_array.use_tiles{n}
            continue
        end
    catch
        disp(['Old tile_array object does not have option to ignore tiles']);
    end
    tn = tiles{n};
    r = tile_array.rows{n};
    c = tile_array.cols{n};
    h = size(tn, 1);
    w = size(tn, 2);
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
    for m=tiles_idx(good_idx)
        r2min = tile_array.rows{m};
        r2max = tile_array.rows{m} + size(tiles{m}, 1) -1;
        c2min = tile_array.cols{m};
        c2max = c2min + size(tiles{m}, 2) -1;
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
    amin = max([rminprime - rmin + 1, 1]);
    amax = min([rmaxprime - rmin + 1,h]);
    bmin = max([cminprime - cmin + 1, 1]);
    bmax = min([cmaxprime - cmin + 1, w]);
    tsub = tn(amin:amax, bmin:bmax,:);
    imcounts(rminprime:rmaxprime,cminprime:cmaxprime) = n;%imcounts(r:(r+h-1),c:(c+w-1)) +1;
    if (rmaxprime-rminprime+1 == size(tsub, 1)) & (cmaxprime-cminprime+1 == size(tsub, 2))
        imm2(rminprime:rmaxprime,cminprime:cmaxprime,:) = tsub;
    else
        rmaxprime = rminprime + size(tsub, 1)-1;
        cmaxprime = cminprime + size(tsub, 2)-1;
        imm2(rminprime:rmaxprime,cminprime:cmaxprime,:) = tsub;
    end
end
% Image pixel intensity is normalized by the number of tiles contributing
% to each pixel
%imm2 = imm2./imcounts;
end
