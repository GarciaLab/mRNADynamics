% StitchSubset.m
% author: Gabriella Martini
% date created: 8/13/20
% date last modified: 8/13/20


function [imm2] = StitchSubset(tile_array, tA_idx, tB_idxs, varargin)


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


%neighbor_image = figure;
%imAx = axes(neighbor_image);
imm2 = zeros(im_rows, im_columns, 3, 'uint16');
imcounts = zeros(im_rows, im_columns, 'uint16');
im_tileIDs = zeros(im_rows, im_columns, 'uint16');

%% First plot tile to be moved
tileA = tile_array.imgs{tA_idx};

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
    tileB = tile_array.imgs{all_indx(j)};
    scale_factor = 6*10^4/max(max(tileB));

    if abs(tA_grid_row - tile_grid_rows(j)) == 1
        RGB_tileB = cat(3, tileB, tileB, zeros(size(tileB), 'uint16'));
        scaled_RGB_tileB = RGB_tileB*scale_factor;
    else
        RGB_tileB = cat(3, zeros(size(tileB), 'uint16'), tileB, tileB);
        scaled_RGB_tileB = RGB_tileB*scale_factor;  
    end
    imm2(plot_rows(j):plot_row_bounds(j), plot_columns(j):plot_column_bounds(j),:) = ...
        imm2(plot_rows(j):plot_row_bounds(j), plot_columns(j):plot_column_bounds(j),:) + scaled_RGB_tileB;
end
imm2(imm2 > 65535) = 65535;

%image(imm2)
%title({'Moving tile shown in magenta.', 'Neighboring tiles shown in cyan and yellow.'})





hold off

   








