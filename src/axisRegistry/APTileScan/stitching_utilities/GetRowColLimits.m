% GetRowColLimits.m
% author: Gabriella Martini
% date created: 8/13/20
% date last modified: 8/13/20


function [r_min, r_max, c_min, c_max] = GetRowColLimits(tile_array, tA_ind,...
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

