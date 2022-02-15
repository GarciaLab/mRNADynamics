
function imm2 = imstitchTileColor(tile_array)
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
    imm2 = zeros(numrows, numcols, 3);
    for n =1:numTiles
        r = tile_array.rows{n};
        c = tile_array.cols{n};
        [h, w] = size(tile_array.imgs{n});
        gr = tile_array.grid_positions{n}(1);
        gc = tile_array.grid_positions{n}(2);
        tile = tile_array.imgs{n};
        if (gr == 1) && (gc == 1)
            %RGB1 = cat(3, imm2, imm2, imm2);  % information stored in intensity
            RGB2 = cat(3, tile, tile, tile)/2;

            RGB2(:, :, 2:3) = 0;  % All information in red channel
            imm2Index = uint8(floor(tile * 255));
            imm2(r:(r+h-1), c:(c+w-1),:) = imm2(r:(r+h-1), c:(c+w-1),:) + RGB2;
        elseif (gr == 1) && (gc == 2)
            %RGB1 = cat(3, imm2, imm2, imm2);  % information stored in intensity
            RGB2 = cat(3, tile, tile, tile)/2;
            RGB2(:, :, 1) = 0;  % All information in green channel
            RGB2(:, :, 3) = 0;  % All information in green channel
            imm2Index = uint8(floor(tile * 255));
            imm2(r:(r+h-1), c:(c+w-1),:) = imm2(r:(r+h-1), c:(c+w-1),:) + RGB2;
        elseif (gr == 1) && (gc == 3)
            %RGB1 = cat(3, imm2, imm2, imm2);  % information stored in intensity
            RGB2 = cat(3, tile, tile, tile)/2;
            RGB2(:, :, 2:3) = 0;  % All information in red channel
            imm2Index = uint8(floor(tile * 255));
            imm2(r:(r+h-1), c:(c+w-1),:) = imm2(r:(r+h-1), c:(c+w-1),:) + RGB2;
        elseif (gr == 2) && (gc == 1)
            %RGB1 = cat(3, imm2, imm2, imm2);  % information stored in intensity

            RGB2 = cat(3, tile, tile, tile)/2;
            RGB2(:, :, 1:2) = 0;  % All information in blue channel
            imm2(r:(r+h-1), c:(c+w-1),3) = RGB2(:,:,3);
        elseif (gr == 2) && (gc == 2)
            %RGB1 = cat(3, imm2, imm2, imm2);  % information stored in intensity
            %imm2Index = uint8(floor(tile * 255));
            %Map       = copper(255);
            %RGB2      = ind2rgb(imm2Index, Map);
            %RGB2 = cat(3, tile, tile, tile);
            RGB2 = cat(3, tile, tile, tile)/2;

            %RGB2(:, :, 2:3) = 0;  % All information in blue channel
            imm2(r:(r+h-1), c:(c+w-1),:) = imm2(r:(r+h-1), c:(c+w-1),:) + RGB2;
        elseif (gr == 2) && (gc == 3)
            %RGB1 = cat(3, imm2, imm2, imm2);  % information stored in intensity

            RGB2 = cat(3, tile, tile, tile)/2;

            RGB2(:, :, 1:2) = 0;  % All information in blue channel
            imm2(r:(r+h-1), c:(c+w-1),3) = RGB2(:,:,3);
        elseif (gr == 3) && (gc == 1)
            %RGB1 = cat(3, imm2, imm2, imm2);  % information stored in intensity

            RGB2 = cat(3, tile, tile, tile)/2;

            RGB2(:, :, 2:3) = 0;  % All information in red channel
            imm2(r:(r+h-1), c:(c+w-1),:) = imm2(r:(r+h-1), c:(c+w-1),:) + RGB2;
        elseif (gr == 3) && (gc == 2)
            %RGB1 = cat(3, imm2, imm2, imm2);  % information stored in intensity

            RGB2 = cat(3, tile, tile, tile)/2;

            RGB2(:, :, 1) = 0;  % All information in green channel
            RGB2(:, :, 3) = 0;  % All information in green channel
            imm2Index = uint8(floor(tile * 255));
            imm2(r:(r+h-1), c:(c+w-1),:) = imm2(r:(r+h-1), c:(c+w-1),:) + RGB2;
        elseif (gr == 3) && (gc == 3)
            %RGB1 = cat(3, imm2, imm2, imm2);  % information stored in intensity

            RGB2 = cat(3, tile, tile, tile)/2;
            RGB2(:, :, 2:3) = 0;  % All information in red channel
            imm2Index = uint8(floor(tile * 255));
            imm2(r:(r+h-1), c:(c+w-1),:) = imm2(r:(r+h-1), c:(c+w-1),:) + RGB2;
        else
            RGB2 = cat(3, tile, tile, tile)/2;
            imm2(r:(r+h-1), c:(c+w-1),:) = imm2(r:(r+h-1), c:(c+w-1),:) + RGB2;
        end
    c = c;
    end
end
