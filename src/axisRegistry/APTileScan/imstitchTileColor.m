
function imm2 = imstitchTileColor(tile_array)
    numTiles = length(tile_array.imgs);
    hs = [];
    ws = [];
    for n=1:numTiles
        hs(n) = size(tile_array.imgs{n}, 1);
        ws(n) = size(tile_array.imgs{n}, 2);
    end
    try
        good_idx = [];
        for i = 1:numTiles
            if tile_array.use_tiles{i}
                good_idx(length(good_idx) + 1) = i;
            end 
        end
    catch
        disp(['Old tile_array object does not have option to ignore tiles']);
        good_idx = 1:numTiles;
    end
    rs = [tile_array.rows{good_idx}];
    cs = [tile_array.cols{good_idx}];
    hs = hs(good_idx);
    ws = ws(good_idx);
    rbounds = rs + hs-1;
    cbounds = cs + ws-1;
    numrows = max(rbounds);
    numcols = max(cbounds);
    imm2 = zeros(numrows, numcols, 3, 'uint16');
    tiles_idx  = 1:numTiles;
    for n =tiles_idx(good_idx)
        r = tile_array.rows{n};
        c = tile_array.cols{n};
        [h, w] = size(tile_array.imgs{n});
        gr = tile_array.grid_positions{n}(1);
        gc = tile_array.grid_positions{n}(2);
        tile = tile_array.imgs{n};
        
        if (mod(gr, 2) == 1) & (mod(gc, 2) == 1)
            scale_factor = 6*10^4/max(max(tile));
            RGB = cat(3, tile, zeros(size(tile), 'uint16'), tile);
            RGB2 = RGB*scale_factor;
            imm2(r:(r+h-1), c:(c+w-1),:) = imm2(r:(r+h-1), c:(c+w-1),:) + RGB2;
        elseif  (mod(gr, 2) == 1) & (mod(gc, 2) == 0)
            scale_factor = 6*10^4/max(max(tile));
            RGB = cat(3, tile, tile, zeros(size(tile), 'uint16'));
            RGB2 = RGB*scale_factor;
            imm2(r:(r+h-1), c:(c+w-1),:) = imm2(r:(r+h-1), c:(c+w-1),:) + RGB2;
        elseif (mod(gr, 2) == 0) & (mod(gc, 2) == 1)
            scale_factor = 6*10^4/max(max(tile));
            RGB = cat(3, tile, tile, zeros(size(tile), 'uint16'));
            RGB2 = RGB*scale_factor;
            imm2(r:(r+h-1), c:(c+w-1),:) = imm2(r:(r+h-1), c:(c+w-1),:) + RGB2;
        else
            scale_factor = 6*10^4/max(max(tile));
            RGB = cat(3, zeros(size(tile), 'uint16'), tile, tile);
            RGB2 = RGB*scale_factor;
            imm2(r:(r+h-1), c:(c+w-1),:) = imm2(r:(r+h-1), c:(c+w-1),:) + RGB2;
        
        end
    end
end
