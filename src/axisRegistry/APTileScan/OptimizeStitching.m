function tile_array = OptimizeStitching(Prefix, ID, MaxDeltaR, MaxDeltaC, varargin)
% OptimizeStitching.m 
% Gabriella Martini
% 1/21/2020
% Last Modifed: 8/15/2020


%% Parse input arguments

manualStitchOrder = false;
selectRegions = false;
useSurfStitching = false;
manualSeeding = false;
if ~isempty(varargin)
    x = 1;
    while x <= length(varargin{1})
        switch varargin{1}{x}
            case{'manualStitchOrder'}
                manualStitchOrder = true;
            case{'selectStitchingRegions'}
                selectRegions=true;
            case{'useSurfStitchingInfo'}
                useSurfStitching=true;
            case{'manualSeeding'}
                manualSeeding=true;
            otherwise
                error('Flag not valid')
        end
        x = x +1;
    end
end


 %% Load existing tile array information and relevant folder info 


%Get relevant folder information 
[SourcePath,FISHPath,DropboxFolder,MS2CodePath]=...
    DetermineLocalFolders(Prefix);

stitchingDataFolder = [DropboxFolder,filesep,Prefix,filesep,'FullEmbryoStitching'];

% load TileArray with stitching information 

    
if exist([stitchingDataFolder, filesep, ID, 'TileArray.mat'], 'file') & ~useSurfStitching
    load([stitchingDataFolder, filesep, ID, 'TileArray.mat']);
elseif (useSurfStitching) & (lower(ID) == 'mid')
    if exist([stitchingDataFolder, filesep, 'SurfTileArray.mat'], 'file')
        load([stitchingDataFolder, filesep, 'Surf', 'TileArray.mat']);
        surf_tile_array = tile_array;
        NewTileArrayFromMetadata(Prefix, ID);
        load([stitchingDataFolder, filesep, ID, 'TileArray.mat']);
        tile_array.rows = surf_tile_array.rows;
        tile_array.cols = surf_tile_array.cols; 
    else
        disp('No MidTileArray data stored. Seed a new TileArray using "NewTileArrayFromMetadata".')
        NewTileArrayFromMetadata(Prefix, ID);
        load([stitchingDataFolder, filesep, ID, 'TileArray.mat']);
    end
else
     disp('No TileArray data stored. Seed a new TileArray using "NewTileArrayFromMetadata".') 
     NewTileArrayFromMetadata(Prefix, ID);
     load([stitchingDataFolder, filesep, ID, 'TileArray.mat']);
end

% if manualSeeding
%     ManualStitchingCorrection(Prefix, ID)
%     load([stitchingDataFolder, filesep, ID, 'TileArray.mat']);
% end
%% 
outputFolder = [DropboxFolder,filesep,Prefix,filesep,'FullEmbryoStitching'];
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end
saveVars = {};
saveVars = [saveVars, 'tile_array'];
outputDatafile = [upper(ID(1)), ID(2:end), 'TileArray.mat'];

% ERROR: NEED TO IMPLEMENT SOMETHING TO REQUIRE THAT INITIAL SEED
% DOESN'T VIOLATE MAX OVERLAP CONDITION

%% 
imm2 = imstitchTile(tile_array);
imshow(imm2)
%% 

NTiles = length(tile_array.imgs);
if ~manualStitchOrder
    stitchOrder = getStitchingOrder(tile_array);
else 
    stitchOrder = ManualStitchOrder(Prefix, ID);
end

for m=2:length(stitchOrder)
    fprintf('Stitching Tile %d of %d\n', [m, NTiles]);
    tA_ind = stitchOrder(m);gr = tile_array.grid_positions{tA_ind}(1);
    gc = tile_array.grid_positions{tA_ind}(2);
    tileA = tile_array.imgs{tA_ind};
    [hA, wA] = size(tileA);
    tile_array.prevrows{size(tile_array.prevrows, 2)+1} = tile_array.rows;
    tile_array.prevcols{size(tile_array.prevcols, 2)+1} = tile_array.cols;
    tAr = tile_array.rows{tA_ind}; 
    tAc = tile_array.cols{tA_ind};
    [tAr_min, tAr_max, tAc_min, tAc_max ] = ...
        GetRowColLimits(tile_array, tA_ind, MaxDeltaR, MaxDeltaC, stitchOrder, m);
    
    rmins = [tile_array.rows{:}];
    rmins(length(rmins)+1) = tAr_min;
    top_limit = min(rmins);
    if top_limit < 1
        tAr_min = tAr_min + (1-top_limit);
        tAr_max = tAr_max + (1-top_limit); 
        for r =1:length(tile_array.rows)
            tile_array.rows{r} = tile_array.rows{r} + (1-top_limit);
        end
    end

    cmins = [tile_array.cols{:}];
    cmins(length(cmins)+1) = tAc_min;
    left_limit = min(cmins);
    if left_limit < 1
        tAc_min = tAc_min + (1-left_limit);
        tAc_max = tAc_max + (1-left_limit); 
        for c =1:length(tile_array.cols)
            tile_array.cols{c} = tile_array.cols{c} + (1-left_limit);
        end
    end
    tArRange = tAr_min:1:tAr_max;
    tAcRange = tAc_min:1:tAc_max;
    tileBs = [];
    for t=1:(m-1)
        if isequal(tile_array.grid_positions{stitchOrder(t)}, [gr+1, gc])
            tileBs(length(tileBs)+1) = stitchOrder(t);
        elseif isequal(tile_array.grid_positions{stitchOrder(t)}, [gr-1, gc])
            tileBs(length(tileBs)+1) = stitchOrder(t);
        elseif isequal(tile_array.grid_positions{stitchOrder(t)}, [gr, gc+1])
            tileBs(length(tileBs)+1) = stitchOrder(t);
        elseif isequal(tile_array.grid_positions{stitchOrder(t)}, [gr, gc-1])
            tileBs(length(tileBs)+1) = stitchOrder(t);
        elseif isequal(tile_array.grid_positions{stitchOrder(t)}, [gr+1, gc+1])
            tileBs(length(tileBs)+1) = stitchOrder(t);
        elseif isequal(tile_array.grid_positions{stitchOrder(t)}, [gr+1, gc-1])
            tileBs(length(tileBs)+1) = stitchOrder(t);
        elseif isequal(tile_array.grid_positions{stitchOrder(t)}, [gr-1, gc+1])
            tileBs(length(tileBs)+1) = stitchOrder(t);
        elseif isequal(tile_array.grid_positions{stitchOrder(t)}, [gr-1, gc-1])
            tileBs(length(tileBs)+1) = stitchOrder(t);
        end
    end
    scores = zeros(length(tArRange), length(tAcRange), length(tileBs));
    areas = zeros(length(tArRange), length(tAcRange), length(tileBs));
    counter = 0;
    
    if ~selectRegions
        for tB_ind =tileBs
            Apos = tile_array.grid_positions{tA_ind};
            Bpos = tile_array.grid_positions{tB_ind};
            counter = counter + 1;
            tileB = tile_array.imgs{tB_ind};
            [hB, wB] = size(tileB);
            tBr = tile_array.rows{tB_ind};
            tBc = tile_array.cols{tB_ind};
            rminB = tBr; 
            rmaxB = tBr+hB-1;
            cminB = tBc;
            cmaxB = tBc+wB-1;
            for i=1:length(tArRange)
                r = tArRange(i);
                rminA = r;
                rmaxA = r+hA-1;
                rminAB = max([rminA, rminB]);
                rmaxAB = min([rmaxA, rmaxB]);
                for j = 1:length(tAcRange)
                    c = tAcRange(j);
                    cminA = c;
                    cmaxA = c+wA-1;
                    cminAB = max([cminA, cminB]);
                    cmaxAB = min([cmaxA, cmaxB]);
                    windowA = tileA(rminAB-r +1:rmaxAB-r, cminAB-c + 1:cmaxAB-c);
                    windowB = tileB(rminAB-tBr+1:rmaxAB-tBr, cminAB-tBc+1:cmaxAB-tBc);
                    diff = abs(windowA-windowB);
                    area = (rmaxAB-rminAB)*(cmaxAB-cminAB);
                    scores(i,j, counter) = sum(diff(:));
                    areas(i,j, counter) = area;
                end
            end
        end
        summed_scores = sum(scores, 3);
        summed_areas = sum(areas, 3);
        normed_scores = summed_scores./summed_areas;


        [newr_ind, newc_ind] = find(normed_scores == min(min(normed_scores)));
        newr = tArRange(newr_ind);
        newc = tAcRange(newc_ind);

        tile_array.rows{tA_ind} = newr;
        tile_array.cols{tA_ind} = newc;
        rmins = [tile_array.rows{:}];
        top_limit = min(rmins);
        for rr =1:length(tile_array.rows)
            tile_array.rows{rr} = tile_array.rows{rr} + (1-top_limit);
        end

        cmins = [tile_array.cols{:}]; 
        left_limit = min(cmins);
        for cc =1:length(tile_array.cols)
            tile_array.cols{cc} = tile_array.cols{cc} + (1-left_limit);
        end
    else
        [tB_inds, skip_tile] = SelectStitchingPartners(tile_array, tA_ind, tileBs);
         
        
        if skip_tile
            rmins = [tile_array.rows{:}];
            top_limit = min(rmins);
            for rr =1:length(tile_array.rows)
                tile_array.rows{rr} = tile_array.rows{rr} + (1-top_limit);
            end

            cmins = [tile_array.cols{:}]; 
            left_limit = min(cmins);
            for cc =1:length(tile_array.cols)
                tile_array.cols{cc} = tile_array.cols{cc} + (1-left_limit);
            end
        else
            [tile_array, TempMaxDeltaR, TempMaxDeltaC] = ManualSeedCorrection(tile_array, tA_ind, tB_inds);
            if isempty(TempMaxDeltaR)
                TempMaxDeltaR = MaxDeltaR;
            end
            if isempty(TempMaxDeltaC)
                TempMaxDeltaC = MaxDeltaC;
            end
            tile_array = ManualTileStitch(tile_array, tA_ind, tB_inds, stitchOrder, TempMaxDeltaR, TempMaxDeltaC);
        end
        %disp('Checking')

    end
    save([outputFolder, filesep, outputDatafile],saveVars{:});
end


%% 

good_rows = [tile_array.rows{stitchOrder}];
min_good_row = min(good_rows);
max_good_row = max(good_rows);
good_columns = [tile_array.cols{stitchOrder}];
min_good_column = min(good_columns);
max_good_column = max(good_columns);
for i=1:NTiles
    if ~ismember(i, stitchOrder)
        tile_array.use_tiles{i} = false;
        if tile_array.rows{i} < min_good_row
            tile_array.rows{i} = min_good_row;
        elseif tile_array.rows{i} > max_good_row
             tile_array.rows{i} = max_good_row;
        end
        
        if tile_array.cols{i} < min_good_column
            tile_array.cols{i} = min_good_column;
        elseif tile_array.cols{i} > max_good_column
             tile_array.cols{i} = max_good_column;
        end
        
    end
    tile_array.rows{i}= tile_array.rows{i} + (1-min_good_row);
    tile_array.cols{i}= tile_array.cols{i} + (1-min_good_column);
end



save([outputFolder, filesep, outputDatafile],saveVars{:});
GenerateStitchedData(Prefix, ID);





end









%% 

function [stitchOrder] = getStitchingOrder(tile_array)
    NTiles = length(tile_array.imgs);
    tileRows = zeros(NTiles, 1);
    tileCols = zeros(NTiles, 1);
    for i=1:NTiles
        tileRows(i) = tile_array.grid_positions{i}(1);
        tileCols(i) = tile_array.grid_positions{i}(2);
    end
    maxRow = max(tileRows);
    maxCol = max(tileCols);

    availableTiles = 1:NTiles;
    stitchOrder = [];
    initRow = ceil(maxRow/2);
    initCol = ceil(maxCol/2);
    for i=1:length(availableTiles)
        n = availableTiles(i);
        if isequal(tile_array.grid_positions{n}, [initRow, initCol])
            stitchOrder(length(stitchOrder)+1) = n;
            availableTiles(i) = [];
            break
        end
    end

    if length(availableTiles) > 0
        cr = initRow+1;
        while cr <= maxRow
            for i=1:length(availableTiles)
                n = availableTiles(i);
                if isequal(tile_array.grid_positions{n}, [cr, initCol])
                    stitchOrder(length(stitchOrder)+1) = n;
                    availableTiles(i) = [];
                    cr = cr+1;
                    break
                end
            end
        end
    end
    if length(availableTiles) > 0
        cr = initRow-1;
        while cr >= 1
            for i=1:length(availableTiles)
                n = availableTiles(i);
                if isequal(tile_array.grid_positions{n}, [cr, initCol])
                    stitchOrder(length(stitchOrder)+1) = n;
                    availableTiles(i) = [];
                    cr = cr-1;
                    break
                end
            end
        end
    end
    if length(availableTiles) > 0
        cc = initCol+1;
        while cc <= maxCol
            for i=1:length(availableTiles)
                n = availableTiles(i);
                if isequal(tile_array.grid_positions{n}, [initRow, cc])
                    stitchOrder(length(stitchOrder)+1) = n;
                    availableTiles(i) = [];
                    cc = cc+1;
                    break
                end
            end
        end
    end
    if length(availableTiles) > 0
        cc = initCol-1;
        while cc >= 1
            for i=1:length(availableTiles)
                n = availableTiles(i);
                if isequal(tile_array.grid_positions{n}, [initRow, cc])
                    stitchOrder(length(stitchOrder)+1) = n;
                    availableTiles(i) = [];
                    cc = cc-1;
                    break
                end
            end
        end
    end
    if length(availableTiles) > 0
        cr = initRow+1;
        cc = initCol-1;
        while cr <= maxRow
            while cc >= 1
                for i=1:length(availableTiles)
                    n = availableTiles(i);
                    if isequal(tile_array.grid_positions{n}, [cr, cc])
                        stitchOrder(length(stitchOrder)+1) = n;
                        availableTiles(i) = [];
                        cc = cc-1;
                        break
                    end
                end
            end
            cr = cr+1;
        end
    end
    if length(availableTiles) > 0
        cr = initRow-1;
        cc = initCol-1;
        while cr >= 1
            while cc >= 1
                for i=1:length(availableTiles)
                    n = availableTiles(i);
                    if isequal(tile_array.grid_positions{n}, [cr, cc])
                        stitchOrder(length(stitchOrder)+1) = n;
                        availableTiles(i) = [];
                        cc = cc-1;
                        break
                    end
                end
            end
            cr = cr-1;
        end
    end
    if length(availableTiles) > 0
        cr = initRow+1;
        cc = initCol+1;
        while cr <= maxRow
            while cc <= maxCol
                for i=1:length(availableTiles)
                    n = availableTiles(i);
                    if isequal(tile_array.grid_positions{n}, [cr, cc])
                        stitchOrder(length(stitchOrder)+1) = n;
                        availableTiles(i) = [];
                        cc = cc+1;
                        break
                    end
                end
            end
            cr = cr+1;
        end
    end
    if length(availableTiles) > 0
        cr = initRow-1;
        cc = initCol+1;
        while cr >= 1
            while cc <= maxCol
                for i=1:length(availableTiles)
                    n = availableTiles(i);
                    if isequal(tile_array.grid_positions{n}, [cr, cc])
                        stitchOrder(length(stitchOrder)+1) = n;
                        availableTiles(i) = [];
                        cc = cc+1;
                        break
                    end
                end
            end
            cr = cr-1;
        end
    end
end







