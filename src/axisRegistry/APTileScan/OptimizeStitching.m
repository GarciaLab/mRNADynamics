function tile_array = OptimizeStitching(Prefix, ID, MaxDeltaR, MaxDeltaC, varargin)
% OptimizeStitching.m 
% Gabriella Martini
% 1/21/2020
% Last Modifed: 7/25/2020

 %% Load existing tile array information and relevant folder info 

NewTileArrayFromMetadata(Prefix, ID);
%Get relevant folder information 
[SourcePath,FISHPath,DropboxFolder,MS2CodePath]=...
    DetermineLocalFolders(Prefix);

stitchingDataFolder = [DropboxFolder,filesep,Prefix,filesep,'FullEmbryoStitching'];

% load TileArray with stitching information 
if exist([stitchingDataFolder, filesep, ID, 'TileArray.mat'], 'file')
     load([stitchingDataFolder, filesep, ID, 'TileArray.mat']);
else
     error('No TileArray data stored. Seed a new TileArray using "NewTileArrayFromMetadata".')  
end


manualStitchOrder = false;
x = 1;
while x <= length(varargin)
    switch varargin{x}
        case{'manualStitchOrder'}
            manualStitchOrder = true;
    end
    x = x +1;
end

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

for m=2:NTiles
    fprintf('Stitching Tile %d of %d\n', [m, NTiles]);
    tA_ind = stitchOrder(m);
    gr = tile_array.grid_positions{tA_ind}(1);
    gc = tile_array.grid_positions{tA_ind}(2);
    tileA = tile_array.imgs{tA_ind};
    [hA, wA] = size(tileA);
    tile_array.prevrows{size(tile_array.prevrows, 2)+1} = tile_array.rows;
    tile_array.prevcols{size(tile_array.prevcols, 2)+1} = tile_array.cols;
    tAr = tile_array.rows{tA_ind}; 
    tAc = tile_array.cols{tA_ind};
    [tAr_min, tAr_max, tAc_min, tAc_max ] = ...
        getRowColLimits(tile_array, tA_ind, MaxDeltaR, MaxDeltaC, stitchOrder, m);
    
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
        end
    end
    scores = zeros(length(tArRange), length(tAcRange), length(tileBs));
    areas = zeros(length(tArRange), length(tAcRange), length(tileBs));
    counter = 0;
    
    for tB_ind =tileBs
        Apos = tile_array.grid_positions{tA_ind};
        Bpos = tile_array.grid_positions{tB_ind};
        counter = counter + 1;
        tileB = tile_array.imgs{tB_ind};
        [hB, wB] = size(tileB);
        tBr = tile_array.rows{tB_ind};
        tBc = tile_array.cols{tB_ind};
        %pB = polyshape([tBr, tBr, tBr+ hB, tBr+hB], [tBc, tBc+wB, tBc+wB, tBc]);
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
                %pA = polyshape([r, r, r+ hA, r+hA], [c, c+wA, c+wA, c]);
                %overlap = intersect(pA, pB);
                %if ~isempty(overlap.Vertices) 
                    %rabsmin = min(overlap.Vertices(:,1));
                    %rabsmax = max(overlap.Vertices(:,1));
                    %cabsmin = min(overlap.Vertices(:,2));
                    %cabsmax = max(overlap.Vertices(:,2));
                    windowA = tileA(rminAB-r +1:rmaxAB-r, cminAB-c + 1:cmaxAB-c);
                    windowB = tileB(rminAB-tBr+1:rmaxAB-tBr, cminAB-tBc+1:cmaxAB-tBc);
                    diff = abs(windowA-windowB);
                    %area = (rabsmax-rabsmin)*(cabsmax-cabsmin);
                    area = (rmaxAB-rminAB)*(cmaxAB-cminAB);
                    scores(i,j, counter) = sum(diff(:));
                    areas(i,j, counter) = area;
                %end
            end
%             if mod(i, 50) == 0
%                 fprintf('%d/%d\n', [i, length(tArRange)]);
%             elseif mod(i, 25) == 0
%                 fprintf('%d/%d', [i, length(tArRange)]);
%             else
%                 fprintf('.');
%             end
        end
    end
    summed_scores = sum(scores, 3);
    summed_areas = sum(areas, 3);
    normed_scores = summed_scores./summed_areas;
    [newr_ind, newc_ind] = find(normed_scores == min(min(normed_scores)));
    newr = tArRange(newr_ind); 
    newc= tAcRange(newc_ind);
    close all
    figure(1)
    imagesc(normed_scores)
    hold on 
    scatter(newc_ind, newr_ind, 100, 'r.') 
    colorbar
    title(['Moving tile: ', num2str(tA_ind)])
    hold off
    figure(2) 
    imagesc(tileA)
    title(['Tile A: ', num2str(tA_ind)])
    figure(3) 
    imagesc(tileB)
    title(['Tile B: ', num2str(tB_ind)])
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
    rmins = rmins;
    
end


%% 
outputFolder = [DropboxFolder,filesep,Prefix,filesep,'FullEmbryoStitching'];
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end


saveVars = {};
saveVars = [saveVars, 'tile_array'];
outputDatafile = [upper(ID(1)), ID(2:end), 'TileArray.mat'];
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


function [r_min, r_max, c_min, c_max] = getRowColLimits(tile_array, tA_ind,...
    MaxDeltaR, MaxDeltaC, stitchOrder, m)
    r_min = tile_array.rows{tA_ind} - MaxDeltaR; 
    r_max = tile_array.rows{tA_ind} + MaxDeltaR; 
    c_min = tile_array.cols{tA_ind} - MaxDeltaC; 
    c_max = tile_array.cols{tA_ind} + MaxDeltaC; 
    gr = tile_array.grid_positions{tA_ind}(1);
    gc = tile_array.grid_positions{tA_ind}(2);
    NTiles = length(tile_array.imgs);
    for n=stitchOrder(1:(m-1))
        if isequal(tile_array.grid_positions{n}, [gr, gc+1])
            [c_min, c_max] = ...
                     correctColLimitsRightTile(tile_array, tA_ind,...
                     n, c_min, c_max);
        elseif isequal(tile_array.grid_positions{n}, [gr,gc-1])
            [c_min, c_max] = ...
                 correctColLimitsLeftTile(tile_array, tA_ind,...
                 n, c_min, c_max);
        elseif isequal(tile_array.grid_positions{n}, [gr+1,gc])
            [r_min, r_max] = ...
                correctRowLimitsLowerTile(tile_array, tA_ind, n,...
                r_min, r_max);
        elseif isequal(tile_array.grid_positions{n}, [gr-1,gc])
            [r_min, r_max] = ...
                correctRowLimitsUpperTile(tile_array, tA_ind, n,...
                r_min, r_max);
        end  
    end
end

function [r_min, r_max] = correctRowLimitsLowerTile(tile_array, tA_ind, n, r_min, r_max)
    hA = size(tile_array.imgs{tA_ind}, 1);
    r_min = max([r_min, (tile_array.rows{n} - hA +1)]);
end

function [r_min, r_max] = correctRowLimitsUpperTile(tile_array, tA_ind, n, r_min, r_max)
    h = size(tile_array.imgs{n}, 1);
    r_max = min([r_max, (tile_array.rows{n} + h-1)]);
end

function [c_min, c_max] = correctColLimitsRightTile(tile_array, tA_ind, n, c_min, c_max)
    wA = size(tile_array.imgs{tA_ind}, 2);
    c_min = max([c_min, (tile_array.cols{n} - wA +1)]);
end

function [c_min, c_max] = correctColLimitsLeftTile(tile_array, tA_ind, n, c_min, c_max)
    w = size(tile_array.imgs{n}, 2);
    c_max = min([c_max, (tile_array.cols{n} + w-1)]);
end







