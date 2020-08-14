% ManualTileStitch.m
% author: Gabriella Martini
% date created: 8/13/20
% date last modified: 8/13/20
% Prefix='2020-08-09-HbJB3-27_5C-Anterior-Embryo1';
% ID = 'Surf';
% tA_idx = 2;
% tB_idxs = [1];
% stitchOrder = [1, 2, 3, 4, 6, 7, 8, 9, 5, 10];
% MaxDeltaR = 150;
% MaxDeltaC = 150;
% 
% %% 
% 
% NewTileArrayFromMetadata(Prefix, ID);
% %Get relevant folder information 
% [SourcePath,FISHPath,DropboxFolder,MS2CodePath]=...
%     DetermineLocalFolders(Prefix);
% 
% stitchingDataFolder = [DropboxFolder,filesep,Prefix,filesep,'FullEmbryoStitching'];
% 
% % load TileArray with stitching information 
% if exist([stitchingDataFolder, filesep, ID, 'TileArray.mat'], 'file')
%      load([stitchingDataFolder, filesep, ID, 'TileArray.mat']);
% else
%      error('No TileArray data stored. Seed a new TileArray using "NewTileArrayFromMetadata".')  
% end


%% 

function [tile_array] = ManualTileStitch(tile_array, tA_idx, tB_idxs, stitchOrder, MaxDeltaR, MaxDeltaC)
close all
NTiles = length(tile_array.rows);
m = find(stitchOrder == tA_idx);
gr = tile_array.grid_positions{tA_idx}(1);
gc = tile_array.grid_positions{tA_idx}(2);
tileA = tile_array.imgs{tA_idx};
[hA, wA] = size(tileA);
tile_array.prevrows{size(tile_array.prevrows, 2)+1} = tile_array.rows;
tile_array.prevcols{size(tile_array.prevcols, 2)+1} = tile_array.cols;
tAr = tile_array.rows{tA_idx}; 
tAc = tile_array.cols{tA_idx};
[tAr_min, tAr_max, tAc_min, tAc_max ] = ...
    GetRowColLimits(tile_array, tA_idx, MaxDeltaR, MaxDeltaC, stitchOrder, m);

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

%% 


hBs = arrayfun(@(x) size(tile_array.tiles{x},1), tB_idxs);
wBs = arrayfun(@(x) size(tile_array.tiles{x},2), tB_idxs);
tB_rows = [tile_array.rows{tB_idxs}];
tB_columns = [tile_array.cols{tB_idxs}];
tileBs = {};
Bpos = {};
for j=1:length(tB_idxs)
    tileBs{j} = tile_array.imgs{tB_idxs(j)};
    Bpos = tile_array.grid_positions{tB_idxs(j)};
end


scores = zeros(length(tArRange), length(tAcRange), length(tB_idxs));
areas = zeros(length(tArRange), length(tAcRange), length(tB_idxs));

Apos = tile_array.grid_positions{tA_idx};  
    
    %% 

for i=1:length(tArRange)
    
    r = tArRange(i);
    rminA = r;
    rmaxA = r+hA-1;
    
    for j = 1:length(tAcRange)
        c = tAcRange(j);
        cminA = c;
        cmaxA = c+wA-1;
        
        
           
        for k =1:length(tB_idxs)
            rminB = tB_rows(k); 
            rmaxB = tB_rows(k)+hBs(k)-1;
            cminB = tB_columns(k);
            cmaxB = tB_columns(k)+hBs(k)-1;
            rminAB = max([rminA, rminB]);
            rmaxAB = min([rmaxA, rmaxB]);
            cminAB = max([cminA, cminB]);
            cmaxAB = min([cmaxA, cmaxB]);
            windowA = tileA(rminAB-r +1:rmaxAB-r, cminAB-c + 1:cmaxAB-c);
            windowB = tileBs{k}(rminAB- tB_rows(k)+1:rmaxAB- tB_rows(k), cminAB- tB_columns(k)+1:cmaxAB-tB_columns(k));
            diff = abs(windowA-windowB);
            area = (rmaxAB-rminAB)*(cmaxAB-cminAB);
            scores(i,j, k) = sum(diff(:));
            areas(i,j, k) = area;

        end
    end
end

summed_scores = sum(scores, 3);
summed_areas = sum(areas, 3);
normed_scores = summed_scores./summed_areas;
%% 
StitchedImage = figure(1);
imAx = axes(StitchedImage);
MinColumn = 1;
MaxColumn = size(normed_scores, 2);
MinRow = 1;
MaxRow = size(normed_scores, 1);

CumScoresImage= figure(2);
scAx = axes(CumScoresImage);

if length(tB_idxs) > 1
    SubScoresImage = figure(3);
    subplots(length(tB_idxs), 1)
end
%% 
if exist('manual_minr', 'var') 
    clear manual_minr
end
if exist('manual_minc', 'var') 
    clear manual_minc
end
cc=1;
cc_value = double(cc);
% imagesc(tileA)

% figure(2) 
% imagesc(tileB)
useScoreMatrix = true;
total_row_shift = 0;
total_column_shift = 0;
temp_tile_array = tile_array;
scores_copy = normed_scores;
scores_mask = zeros(size(normed_scores));
scores_mask(MinRow:MaxRow, MinColumn:MaxColumn) = 1;
scores_copy(scores_mask == 0) = max(max(scores_copy));
[minr, minc] = find(scores_copy == min(min(scores_copy)));

newr = tArRange(minr);
newc = tAcRange(minc);
temp_tile_array.rows{tA_idx} = newr;
temp_tile_array.cols{tA_idx} = newc;
rmins = [temp_tile_array.rows{:}];
top_limit = min(rmins);
for rr =1:length(temp_tile_array.rows)
    temp_tile_array.rows{rr} = temp_tile_array.rows{rr} + (1-top_limit);
end

cmins = [temp_tile_array.cols{:}]; 
left_limit = min(cmins);
for cc =1:length(temp_tile_array.cols)
    temp_tile_array.cols{cc} = temp_tile_array.cols{cc} + (1-left_limit);
end

while (cc~='x')
    subStitchedImage = StitchSubset(temp_tile_array, tA_idx, tB_idxs);
    figure(StitchedImage)
    image(subStitchedImage)
    if ~((total_row_shift == 0) & (total_column_shift == 0))
        title({['Total Manual RowShift: ', num2str(total_row_shift)], ['Total Manual Column Shift: ', num2str(total_column_shift)]})
    end
    figure(CumScoresImage)
    imagesc(normed_scores)
    colorbar
    hold on 
    scatter(minc, minr, 100, 'r.') 
    colorbar

    xline(MinColumn,'r')

    xline(MaxColumn,'r')
    yline(MinRow,'r')

    yline(MaxRow,'r')
    if exist('manual_minr', 'var') & exist('manual_minc', 'var')
       scatter(manual_minc, manual_minr, 100, 'r*')  
    end
    hold off
    if useScoreMatrix
        

        figure(CumScoresImage)
        ct=waitforbuttonpress;
        cc=get(CumScoresImage,'currentcharacter');
        cc_value = double(cc);
        cm=get(scAx,'CurrentPoint');


        if (ct~=0)&(cc=='c')        %Clear all AP information
            MinColumn = 1;
            MaxColumn = size(normed_scores, 2);
            MinRow = 1;
            MaxRow = size(normed_scores, 1);
        elseif (ct~=0)&(cc=='l')	%Select anterior end
            [MinColumn,coordLy]=ginputc(1,'Color',[1,1,1]);
        elseif (ct~=0)&(cc=='r')    %Select posterior end
            [MaxColumn,coordRy]=ginputc(1,'Color',[1,1,1]);
        elseif (ct~=0)&(cc=='t')    %Select posterior end
            [coordTx,MinRow]=ginputc(1,'Color',[1,1,1]);
        elseif (ct~=0)&(cc=='b')    %Select posterior end
            [coordBx,MaxRow]=ginputc(1,'Color',[1,1,1]);
        elseif (ct~=0) & (cc == 's')
            useScoreMatrix = false;
            cc = 1;
        end

        MaxRow = uint16(round(MaxRow));
        MinRow = uint16(round(MinRow));
        MaxColumn = uint16(round(MaxColumn));
        MinColumn = uint16(round(MinColumn));
        
        disp('Using correlation matrix')
        scores_copy = normed_scores;
        scores_mask = zeros(size(normed_scores));
        scores_mask(MinRow:MaxRow, MinColumn:MaxColumn) = 1;
        scores_copy(scores_mask == 0) = max(max(scores_copy));
        [minr, minc] = find(scores_copy == min(min(scores_copy)));

        newr = tArRange(minr);
        newc = tAcRange(minc);
        temp_tile_array.rows{tA_idx} = newr;
        temp_tile_array.cols{tA_idx} = newc;
        rmins = [temp_tile_array.rows{:}];
        top_limit = min(rmins);
        for rr =1:length(temp_tile_array.rows)
            temp_tile_array.rows{rr} = temp_tile_array.rows{rr} + (1-top_limit);
        end

        cmins = [temp_tile_array.cols{:}]; 
        left_limit = min(cmins);
        for cx =1:length(temp_tile_array.cols)
            temp_tile_array.cols{cx} = temp_tile_array.cols{cx} + (1-left_limit);
        end
    else
       disp('Using stitched image')
       temprows = [temp_tile_array.rows{:}];
       tempcols = [temp_tile_array.cols{:}];
       tempheights = arrayfun(@(x) size(temp_tile_array.tiles{x},1), 1:NTiles);
       tempwidths = arrayfun(@(x) size(temp_tile_array.tiles{x},2), 1:NTiles);
       temprmaxs = temprows+tempheights-1;
       tempcmaxs = tempcols + tempwidths - 1;
       figure(StitchedImage) 
       
       ct=waitforbuttonpress;
       cc=get(StitchedImage,'currentcharacter');
       cc_value = double(cc);
       cm=get(imAx,'CurrentPoint');
       
       if (ct~=0) && (cc_value == 30)
            %disp('left arrow');
            temprows(tA_idx) = temprows(tA_idx)-1;
            temprmaxs(tA_idx) = temprmaxs(tA_idx)-1;
            rshift = -1;
            cshift = 0;
        elseif (ct~=0) && (cc_value == 31)
            %disp('right arrow')
            temprows(tA_idx) = temprows(tA_idx)+1;
            temprmaxs(tA_idx) = temprmaxs(tA_idx)+1;
            rshift = 1;
            cshift = 0;
        elseif (ct~=0) && (cc_value == 28)
            %disp('up arrow')
            tempcols(tA_idx) = tempcols(tA_idx)-1;
            tempcmaxs(tA_idx) = tempcmaxs(tA_idx)-1;
            rshift = 0;
            cshift = -1;
        elseif (ct~=0) && (cc_value == 29)
            %disp('down arrow')
            tempcols(tA_idx) = tempcols(tA_idx)+1;
            tempcmaxs(tA_idx) = tempcmaxs(tA_idx)+1;
            rshift = 0;
            cshift=1;
        elseif (ct~=0) && (cc_value == 43)
            %disp('-')
            temprows(tA_idx) = temprows(tA_idx)-10;
            temprmaxs(tA_idx) = temprmaxs(tA_idx)-10;
            rshift = -10;
            cshift = 0;
        elseif (ct~=0) && (cc_value == 45)
            %disp('+')

            temprows(tA_idx) = temprows(tA_idx)+10;
            temprmaxs(tA_idx) = temprmaxs(tA_idx)+10;
            rshift = 10;
            cshift = 0;
        elseif (ct~=0) && (cc_value == 62)
            %disp('>')
            tempcols(tA_idx) = tempcols(tA_idx)+10;
            tempcmaxs(tA_idx) = tempcmaxs(tA_idx)+10;
            rshift = 0;
            cshift= 10;

        elseif (ct~=0) && (cc_value == 60)
            %disp('<')
            tempcols(tA_idx) = tempcols(tA_idx)-10;
            tempcmaxs(tA_idx) = tempcmaxs(tA_idx)-10;
            rshift = 0;
            cshift = -10;

        elseif (ct~=0) && (cc=='s')
            useScoreMatrix = true;
            cc = 1;
            if exist('manual_minr', 'var') 
                clear manual_minr
            end
            if exist('manual_minc', 'var') 
                clear manual_minc
            end
        end
        rmin = min(temprows);
        cmin = min(tempcols);
        temprows = temprows-rmin + 1;
        tempcols = tempcols-cmin+1;
        temprmaxs = temprmaxs-rmin+1;
        tempcmaxs = tempcmaxs-cmin+1;
        for t=1:NTiles
            temp_tile_array.rows{t} = temprows(t);
            temp_tile_array.cols{t} = tempcols(t);

        end
        total_row_shift = total_row_shift + rshift;
        total_column_shift = total_column_shift + cshift;
        manual_minr = minr + total_row_shift;
        manual_minc = minc + total_column_shift;
    end

end
close all

%% 
tile_array = temp_tile_array;

