% ManualSeedCorrection.m
% author: Gabriella Martini
% date created: 8/18/20
% date last modified: 8/18/20

% Prefix = '2020-08-10-HbJB3-27_5C-Anterior-Embryo1';
% ID = 'Surf';
% 
% [SourcePath,FISHPath,DropboxFolder,MS2CodePath]=...
%     DetermineLocalFolders(Prefix);
% 
% stitchingDataFolder = [DropboxFolder,filesep,Prefix,filesep,'FullEmbryoStitching'];
% load([stitchingDataFolder, filesep, ID, 'TileArray.mat']);
% 
% tA_idx = 6;
% tB_idxs = [2];

%% 
function [temp_tile_array, TempMaxDeltaR, TempMaxDeltaC] = ManualSeedCorrection(tile_array, tA_idx, tB_idxs)
rmins = [tile_array.rows{:}];
top_limit = min(rmins);

cmins = [tile_array.cols{:}];
left_limit = min(cmins);


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


temp_tile_array = tile_array;
total_row_shift = 0;
total_column_shift = 0;
NTiles = length(temp_tile_array.tiles);
StitchedImage = figure;
imAx = axes(StitchedImage);
title('Shift Selected Tile (shown in magenta)')
rshift = 0;
cshift = 0;
cc = 1;
while (cc~='x')
    temprows = [temp_tile_array.rows{:}];
    tempcols = [temp_tile_array.cols{:}];
    tempheights = arrayfun(@(x) size(temp_tile_array.tiles{x},1), 1:NTiles);
    tempwidths = arrayfun(@(x) size(temp_tile_array.tiles{x},2), 1:NTiles);
    temprmaxs = temprows+tempheights-1;
    tempcmaxs = tempcols + tempwidths - 1;
    figure(StitchedImage)
    subStitchedImage = StitchSubset(temp_tile_array, tA_idx, tB_idxs);
    image(subStitchedImage)
    axis off
    title({'Shift Selected Tile (shown in magenta)',...
        ['total row shift: ', num2str(total_row_shift), ', total column shift: ', num2str(total_column_shift)]})
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
    
end

close all

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



row_prompt = ['Enter MaxDeltaR (press return to use default value):'];
TempMaxDeltaR = input(row_prompt);
column_prompt = ['Enter MaxDeltaC (press return to use default value):'];
TempMaxDeltaC = input(column_prompt);
  

