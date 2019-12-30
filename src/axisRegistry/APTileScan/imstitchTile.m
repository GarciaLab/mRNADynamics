function imm2 = imstitchTile(tile_array)
    numTiles = length(tile_array.imgs);
    hs = [];
    ws = [];
    for n=1:numTiles
        hs(n) = size(tile_array.imgs{n}, 1);
        ws(n) = size(tile_array.imgs{n}, 2);
    end

    rs = [tile_array.rows{:}];
    cs = [tile_array.cols{:}];
    rbounds = rs + hs-1;
    cbounds = cs + ws-1;
    numrows = max(rbounds);
    numcols = max(cbounds);
    imm2 = zeros(numrows, numcols);
    imcounts = zeros(numrows, numcols);
    for n =1:numTiles
        r = tile_array.rows{n};
        c = tile_array.cols{n};
        [h, w] = size(tile_array.imgs{n});
        imcounts(r:(r+h-1),c:(c+w-1)) = imcounts(r:(r+h-1),c:(c+w-1)) +1;
        imm2(r:(r+h-1), c:(c+w-1)) =...
            imm2(r:(r+h-1), c:(c+w-1))+tile_array.imgs{n};
    end
    imm2 = imm2./imcounts;
end
