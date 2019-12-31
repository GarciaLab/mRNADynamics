function tile_array = FindStitchingPositions(Prefix, ID, MaxStep, MaxOverlap, NIterations)
% author: Gabriella Martini
% date created: 12/30/19
% date last modified: 12/30/19

%% Parse inputs

if ~exist('Prefix')
     FolderTemp=uigetdir(DropboxFolder,'Choose folder with files to analyze');
     Dashes=strfind(FolderTemp,filesep);
     Prefix=FolderTemp((Dashes(end)+1):end);
end
if ~exist('ID', 'var')
    prompt = ['Enter an "ID" for stitching.',...
    'The standard inputs "Mid" and "Surf" will stitch the',...
    'Midsagittal and Surface Full Embryo images using the',...
    '"MidTile.lif" and "SurfTile.lif" files respectively.'];
    ID = input(prompt,'s');

end
ID = [upper(ID(1)), ID(2:end)];
if ~exist('MaxStep', 'var')
    prompt = 'Choose the maximum step size to use for a single loop iteration: ';
    MaxStep = input(prompt);
end
if ~exist('MaxOverlap', 'var')
    prompt = 'Enter the MaxOverlap between adjacent tiles: ';
    MaxOverlap = input(prompt);
end
if ~exist('NIterations', 'var')
    prompt = 'Enter the Number of Iterations ("NIterations") to be used in loop: ';
    NIterations = input(prompt);
end

%% Load existing tile array information and relevant folder info 

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

%% 

imm2 = imstitchTile(tile_array);
imshow(imm2)
    % Start by moving tile1 relative to tiles 2, 3, and 4
    % allow x and y positions that are within 100 of the current value, without
    % failing to overlap neighbors and sstaying within prescribed limits:
    NTiles = length(tile_array.imgs);
    for iter=1:NIterations
        for tA_ind=1:NTiles
            tileA = tile_array.imgs{tA_ind};
            [hA, wA] = size(tileA);
            tile_array.prevrows{size(tile_array.prevrows, 2)+1} = tile_array.rows;
            tile_array.prevcols{size(tile_array.prevcols, 2)+1} = tile_array.cols;
            tAr = tile_array.rows{tA_ind}; tAc = tile_array.cols{tA_ind};
            [tAr_min, tAr_max, tAc_min, tAc_max ] = ...
                getRowColLimits(tile_array, tA_ind, MaxOverlap, MaxStep);

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
            scores = zeros(length(tArRange), length(tAcRange), length(tile_array.imgs)-1);
            areas = zeros(length(tArRange), length(tAcRange), length(tile_array.imgs)-1);
            counter = 0;
            %figure(2)
            for tB_ind =1:length(tile_array.imgs)
                if tA_ind == tB_ind
                    continue
                end
                Apos = tile_array.grid_positions{tA_ind};
                Bpos = tile_array.grid_positions{tB_ind};
                if ~(all(abs(Apos-Bpos)<=1))
                    continue
                end
                counter = counter + 1;
                tileB = tile_array.imgs{tB_ind};
                [hB, wB] = size(tileB);
                tBr = tile_array.rows{tB_ind};
                tBc = tile_array.cols{tB_ind};

                pB = polyshape([tBr, tBr, tBr+ hB, tBr+hB], [tBc, tBc+wB, tBc+wB, tBc]);
                for i=1:length(tArRange)
                    r = tArRange(i);
                    for j = 1:length(tAcRange)
                        c = tAcRange(j);
                        pA = polyshape([r, r, r+ hA, r+hA], [c, c+wA, c+wA, c]);
                        overlap = intersect(pA, pB);
                        rabsmin = min(overlap.Vertices(:,1));
                        rabsmax = max(overlap.Vertices(:,1));
                        cabsmin = min(overlap.Vertices(:,2));
                        cabsmax = max(overlap.Vertices(:,2));
                        windowA = tileA(rabsmin-r +1:rabsmax-r, cabsmin-c + 1:cabsmax-c);
                        windowB = tileB(rabsmin-tBr+1:rabsmax-tBr, cabsmin-tBc+1:cabsmax-tBc);
                        diff = abs(windowA-windowB);
                        area = (rabsmax-rabsmin)*(cabsmax-cabsmin);
                        scores(i,j, counter) = sum(diff(:));
                        areas(i,j, counter) = area;
                    end    
                end
            end
            summed_scores = sum(scores, 3);
            summed_areas = sum(areas, 3);
            normed_scores = summed_scores./summed_areas;
            [newr_ind, newc_ind] = find(normed_scores == min(min(normed_scores)));
            newr = tArRange(newr_ind); newc= tAcRange(newc_ind);
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

        end
        numIters = length(tile_array.prevrows);
        oldrows = [tile_array.prevrows{numIters-(NTiles-1)}{:}];
        oldcols = [tile_array.prevcols{numIters-(NTiles-1)}{:}];
        newrows = [tile_array.rows{:}];
        newcols = [tile_array.cols{:}];
        if mod(iter, 10) == 0
            fprintf('%d/%d', [iter, NIterations]);
        else
            fprintf('.');
        end
        if mod(iter, 50) == 0
            fprintf('\n');
        end
        
        if (isequal(oldrows, newrows) && isequal(oldcols, newcols))
            fprintf('Done!\n');
            imm = imstitchTile(tile_array);
            figure            
            imshow(imm)

            break
        end
        
    end
    if ~(isequal(oldrows, newrows) || ~isequal(oldcols, newcols))
        disp('Tile Positions failed to stabilize. Saving current stitching configuration.');
        imm = imstitchTile(tile_array);
        figure            
        imshow(imm)
        %Create the output folder if it doesn't exist 
        
    end
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


function [r_min, r_max, c_min, c_max] = getRowColLimits(tile_array, tA_ind, MaxOverlap, MaxStep)
    r_min = tile_array.rows{tA_ind} - MaxStep; 
    r_max = tile_array.rows{tA_ind} + MaxStep; 
    c_min = tile_array.cols{tA_ind} - MaxStep; 
    c_max = tile_array.cols{tA_ind} + MaxStep; 
    NTiles = length(tile_array.imgs);
    if isequal(tile_array.grid_positions{tA_ind}, [1,1])
        for n=1:NTiles
            if isequal(tile_array.grid_positions{n}, [1,2])
                 [c_min, c_max] = ...
                     correctColLimitsRightTile(tile_array, tA_ind,...
                     n, c_min, c_max, MaxOverlap);
            elseif isequal(tile_array.grid_positions{n}, [2,1])
                [r_min, r_max] = ...
                    correctRowLimitsLowerTile(tile_array, tA_ind, n,...
                    r_min, r_max, MaxOverlap);
            end    
        end
    elseif isequal(tile_array.grid_positions{tA_ind}, [1,2])
        for n=1:NTiles
            if isequal(tile_array.grid_positions{n}, [1,1])
                [c_min, c_max] = ...
                     correctColLimitsLeftTile(tile_array, tA_ind,...
                     n, c_min, c_max, MaxOverlap);
            elseif isequal(tile_array.grid_positions{n}, [1,3])
                [c_min, c_max] = ...
                     correctColLimitsRightTile(tile_array, tA_ind,...
                     n, c_min, c_max, MaxOverlap);
            elseif isequal(tile_array.grid_positions{n}, [2,2])
                [r_min, r_max] = ...
                    correctRowLimitsLowerTile(tile_array, tA_ind, n,...
                    r_min, r_max, MaxOverlap);
            end    
        end   
    elseif isequal(tile_array.grid_positions{tA_ind}, [1,3])
        for n=1:NTiles
            if isequal(tile_array.grid_positions{n}, [1,2])
                [c_min, c_max] = ...
                     correctColLimitsLeftTile(tile_array, tA_ind,...
                     n, c_min, c_max, MaxOverlap);
            elseif isequal(tile_array.grid_positions{n}, [2,3])
                [r_min, r_max] = ...
                    correctRowLimitsLowerTile(tile_array, tA_ind, n,...
                    r_min, r_max, MaxOverlap);
            end    
        end
    elseif isequal(tile_array.grid_positions{tA_ind}, [2,1])
        for n=1:NTiles
            if isequal(tile_array.grid_positions{n}, [2,2])
                [c_min, c_max] = ...
                     correctColLimitsRightTile(tile_array, tA_ind,...
                     n, c_min, c_max, MaxOverlap);
            elseif isequal(tile_array.grid_positions{n}, [3,1])
                [r_min, r_max] = ...
                    correctRowLimitsLowerTile(tile_array, tA_ind, n,...
                    r_min, r_max, MaxOverlap);
            elseif isequal(tile_array.grid_positions{n}, [1,1])
                [r_min, r_max] = ...
                    correctRowLimitsUpperTile(tile_array, tA_ind, n,...
                    r_min, r_max, MaxOverlap);
            end    
        end
    elseif isequal(tile_array.grid_positions{tA_ind}, [2,2])
        for n=1:NTiles
            if isequal(tile_array.grid_positions{n}, [2,3])
                [c_min, c_max] = ...
                     correctColLimitsRightTile(tile_array, tA_ind,...
                     n, c_min, c_max, MaxOverlap);
            elseif isequal(tile_array.grid_positions{n}, [2,1])
                [c_min, c_max] = ...
                     correctColLimitsLeftTile(tile_array, tA_ind,...
                     n, c_min, c_max, MaxOverlap);
            elseif isequal(tile_array.grid_positions{n}, [3,2])
                [r_min, r_max] = ...
                    correctRowLimitsLowerTile(tile_array, tA_ind, n,...
                    r_min, r_max, MaxOverlap);
            elseif isequal(tile_array.grid_positions{n}, [1,2])
                [r_min, r_max] = ...
                    correctRowLimitsUpperTile(tile_array, tA_ind, n,...
                    r_min, r_max, MaxOverlap);
            end    
        end
    elseif isequal(tile_array.grid_positions{tA_ind}, [2,3])
        for n=1:NTiles
            if isequal(tile_array.grid_positions{n}, [2,2])
                [c_min, c_max] = ...
                     correctColLimitsLeftTile(tile_array, tA_ind,...
                     n, c_min, c_max, MaxOverlap);
            elseif isequal(tile_array.grid_positions{n}, [3,3])
                [r_min, r_max] = ...
                    correctRowLimitsLowerTile(tile_array, tA_ind, n,...
                    r_min, r_max, MaxOverlap);
            elseif isequal(tile_array.grid_positions{n}, [1,3])
                [r_min, r_max] = ...
                    correctRowLimitsUpperTile(tile_array, tA_ind, n,...
                    r_min, r_max, MaxOverlap);
            end    
        end
    elseif isequal(tile_array.grid_positions{tA_ind}, [3,1])
        for n=1:NTiles
            if isequal(tile_array.grid_positions{n}, [3,2])
                [c_min, c_max] = ...
                     correctColLimitsRightTile(tile_array, tA_ind,...
                     n, c_min, c_max, MaxOverlap);
            elseif isequal(tile_array.grid_positions{n}, [2,1])
                [r_min, r_max] = ...
                    correctRowLimitsUpperTile(tile_array, tA_ind, n,...
                    r_min, r_max, MaxOverlap);
            end    
        end
    elseif isequal(tile_array.grid_positions{tA_ind}, [3,2])
        for n=1:NTiles
            if isequal(tile_array.grid_positions{n}, [3,3])
                [c_min, c_max] = ...
                     correctColLimitsRightTile(tile_array, tA_ind,...
                     n, c_min, c_max, MaxOverlap);
            elseif isequal(tile_array.grid_positions{n}, [3,1])
                [c_min, c_max] = ...
                     correctColLimitsLeftTile(tile_array, tA_ind,...
                     n, c_min, c_max, MaxOverlap);
            elseif isequal(tile_array.grid_positions{n}, [2,2])
                [r_min, r_max] = ...
                    correctRowLimitsUpperTile(tile_array, tA_ind, n,...
                    r_min, r_max, MaxOverlap);
            end    
        end
    elseif isequal(tile_array.grid_positions{tA_ind}, [3,3])
        for n=1:NTiles
            if isequal(tile_array.grid_positions{n}, [3,2])
                [c_min, c_max] = ...
                     correctColLimitsLeftTile(tile_array, tA_ind,...
                     n, c_min, c_max, MaxOverlap);
            elseif isequal(tile_array.grid_positions{n}, [2,3])
                [r_min, r_max] = ...
                    correctRowLimitsUpperTile(tile_array, tA_ind, n,...
                    r_min, r_max, MaxOverlap);
            end    
        end        
    else
        error('need to write other cases')
    end

end

function [r_min, r_max] = correctRowLimitsLowerTile(tile_array, tA_ind, n, r_min, r_max, MaxOverlap)
    hA = size(tile_array.imgs{tA_ind}, 1);
    r_min = max([r_min, (tile_array.rows{n} - hA +1)]);
    r_max = min([r_max, (tile_array.rows{n} - hA+MaxOverlap)]);
end

function [r_min, r_max] = correctRowLimitsUpperTile(tile_array, tA_ind, n, r_min, r_max, MaxOverlap)
    h = size(tile_array.imgs{n}, 1);
    r_min = max([r_min, (tile_array.rows{n} + h - MaxOverlap)]);
    r_max = min([r_max, (tile_array.rows{n} + h-1)]);
end

function [c_min, c_max] = correctColLimitsRightTile(tile_array, tA_ind, n, c_min, c_max, MaxOverlap)
    wA = size(tile_array.imgs{tA_ind}, 2);
    c_min = max([c_min, (tile_array.cols{n} - wA +1)]);
    c_max = min([c_max, (tile_array.cols{n} - wA+MaxOverlap)]);
end

function [c_min, c_max] = correctColLimitsLeftTile(tile_array, tA_ind, n, c_min, c_max, MaxOverlap)
    w = size(tile_array.imgs{n}, 2);
    c_min = max([c_min, (tile_array.cols{n} + w - MaxOverlap)]);
    c_max = min([c_max, (tile_array.cols{n} + w-1)]);
end
