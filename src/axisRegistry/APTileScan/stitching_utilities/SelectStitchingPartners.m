% SelectStitchingPartners.m
% author: Gabriella Martini
% date created: 8/13/20
% date last modified: 8/18/20



function [tB_idx_return, skip_tile] = SelectStitchingPartners(tile_array, tA_idx, tB_idxs, varargin)
%% 
 
useFilteredImages = true;
if ~isempty(varargin) > 0
    x = 1;
    while x <= length(varargin)
        switch varargin{x}
            case{'useRawImages'}
                useFilteredImages = false;
            otherwise
                error('Flag not recognized')
        end
        x = x +1;
    end
end
%% 

skip_tile = false;

% Load positional information for potential tiles involved in stitching
NTiles = length(tile_array.tiles);
all_indx = [[tA_idx] tB_idxs];
tile_rows = [tile_array.rows{all_indx}];
tile_columns = [tile_array.cols{all_indx}];
tile_heights = arrayfun(@(x) size(tile_array.tiles{x},1), all_indx);
tile_widths = arrayfun(@(x) size(tile_array.tiles{x},2), all_indx);
tile_grid_rows = arrayfun(@(x) tile_array.grid_positions{x}(1), all_indx);
tile_grid_columns = arrayfun(@(x) tile_array.grid_positions{x}(2), all_indx);
row_min = min(tile_rows);
column_min = min(tile_columns);


%% Make image of all tiles

plot_rows = tile_rows - row_min + 1;
plot_columns = tile_columns - column_min + 1;
plot_row_bounds = plot_rows + tile_heights - 1;
plot_column_bounds = plot_columns + tile_widths -1;
im_rows = max(plot_row_bounds);
im_columns = max(plot_column_bounds);



%Now, select which tiles to use. 
cc=1;

neighbor_image = figure(1);
imAx = axes(neighbor_image);
tileB_incolor = true(1, length(tB_idxs));
%iters = 0;
while (cc~='x')
    %iters = iters+1;
    %disp([num2str(iters)])
    imm2 = zeros(im_rows, im_columns, 3, 'uint16');
    imcounts = zeros(im_rows, im_columns, 'uint16');
    im_tileIDs = zeros(im_rows, im_columns, 'uint16');

    %% First plot tile to be moved
    if useFilteredImages
        tileA = tile_array.imgs{tA_idx};
    else
        tileA = tile_array.tiles{tA_idx};
    end
    % Plot this in magenta
    scale_factor = 6*10^4/max(max(tileA));
    RGB_tileA = cat(3, tileA, zeros(size(tileA), 'uint16'), tileA);
    scaled_RGB_tileA = RGB_tileA*scale_factor;
    imm2(plot_rows(1):plot_row_bounds(1), plot_columns(1):plot_column_bounds(1),:) = ...
        imm2(plot_rows(1):plot_row_bounds(1), plot_columns(1):plot_column_bounds(1),:) + scaled_RGB_tileA;

    imcounts(plot_rows(1):plot_row_bounds(1), plot_columns(1):plot_column_bounds(1)) = ...
        imcounts(plot_rows(1):plot_row_bounds(1), plot_columns(1):plot_column_bounds(1)) + 1;
    im_tileIDs(plot_rows(1):plot_row_bounds(1), plot_columns(1):plot_column_bounds(1)) = 1;
    %% Next plot neighbor tiles
    tA_grid_row = tile_grid_rows(1);
    tA_grid_column = tile_grid_columns(1);

    for j=2:length(all_indx)
        if useFilteredImages 
            tileB = tile_array.imgs{all_indx(j)};
        else
            tileB = tile_array.tiles{all_indx(j)};
        end
        scale_factor = 6*10^4/max(max(tileB));
        if tileB_incolor(j-1)
            if ~((abs(tA_grid_row - tile_grid_rows(j)) == 1)&(abs(tA_grid_column - tile_grid_columns(j)) == 1))
                RGB_tileB = cat(3, tileB, tileB, zeros(size(tileB), 'uint16'));
                scaled_RGB_tileB = RGB_tileB*scale_factor;
            else
                RGB_tileB = cat(3, zeros(size(tileB), 'uint16'), tileB, tileB);
                scaled_RGB_tileB = RGB_tileB*scale_factor;  
            end
        else
           RGB_tileB = cat(3, tileB, tileB, tileB);
           scaled_RGB_tileB = RGB_tileB*scale_factor; 
        end
        imm2(plot_rows(j):plot_row_bounds(j), plot_columns(j):plot_column_bounds(j),:) = ...
            imm2(plot_rows(j):plot_row_bounds(j), plot_columns(j):plot_column_bounds(j),:) + scaled_RGB_tileB;
        imcounts(plot_rows(j):plot_row_bounds(j), plot_columns(j):plot_column_bounds(j)) = ...
            imcounts(plot_rows(j):plot_row_bounds(j), plot_columns(j):plot_column_bounds(j)) + 1;
        im_tileIDs(plot_rows(j):plot_row_bounds(j), plot_columns(j):plot_column_bounds(j)) = j;
    end
    %disp('Checkpoint 1')
    imm2(imm2 > 65535) = 65535;
    figure(neighbor_image)
    image(imAx, imm2)
    title({'Moving tile shown in magenta.', 'Selected neighboring tiles shown in cyan and yellow. Unselected tiles in gray.'})


    figure(neighbor_image)
    ct=waitforbuttonpress;
    cc=get(neighbor_image,'currentcharacter');
    cm=get(imAx,'CurrentPoint');
    
    
    if (ct~=0) & (cc == 'm')
        [select_column, select_row]=ginputc(1, 'Color',[1,1,1]);
        select_row = uint16(round(select_row));
        select_column = uint16(round(select_column));
        if imcounts(select_row, select_column) == 1
            select_idx = im_tileIDs(select_row, select_column);
            if select_idx > 1
                if tileB_incolor(select_idx-1)
                    tileB_incolor(select_idx-1) = false;
                else
                    tileB_incolor(select_idx-1) = true;
                end
                    
            end
        end
    elseif (cc == 's')
        if skip_tile
            skip_tile = false;
        else
            skip_tile = true;
        end
    elseif (cc == 'x')
        if (length(tB_idxs) == 0) & ~skip_tile
            cc = 1;
            disp('Invalid choice of tiles. Must select at least one neighboring tile for stitching.')
        end
    end
    

end

close all
tB_idx_return = tB_idxs(tileB_incolor);
end




