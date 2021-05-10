function [schnitzcells, Ellipses] = readimariscsv()
%%
T = readtable("P:\Liz\MATLAB\Data\RawDynamicsData\2020-11-26\h2bmScar-MCPmed\ImarisResult\h2bmScar-MCPmed_Statistics\h2bmScar-MCPmed_Position.csv");

%drop needless features
T = removevars(T, {'Unit', 'Category', 'Collection', 'Var10'});

%groupby nucleus ID and aggregate into a list. 
T_groupedByID = groupsummary(T, 'TrackID', @(x) {x});

%I think these are all stray nuclei that weren't tracked through multiple
%frames.
T_groupedByID = rmmissing(T_groupedByID);

%unfortunately, the times are not sorted. let's fix that 
for row = 1:height(T_groupedByID)
    timevec = T_groupedByID(1, 'fun1_Time').fun1_Time{1};
    
    for col = 1:width(T_groupedByID)
        val = T_groupedByID.(col)(row);
        if iscell(val)
            g = [timevec, cell2mat(T_groupedByID.(col)(row))];
            newg = sortrows(g, 1);
            T_groupedByID.(col)(row) = {newg(:, 2)};
        end
    end
    
end

%modulo column names, this table is basically schnitzcells. 
%Also, the positions still need to be converted from um to voxels. 
schnitzcells = table2struct(T_groupedByID);

%Now, we'll make Ellipses by grouping by the Time column 
T_groupedByTime = groupsummary(T, 'Time', @(x) {x});

%We need to place the resulting data into the exact format of Ellipses
Ellipses = cell(length(timevec), 1);
for frame = 1:height(T_groupedByTime)
    val = T_groupedByTime.(3)(frame);
    vec = val{1};
    Ellipses{frame} = nan(numel(vec), 1);
    for col = 3:width(T_groupedByTime)
        val = T_groupedByTime.(col)(frame);
        vec = val{1};
        Ellipses{frame}(:, col-2) = vec;
    end
   
end

%Note that the orders the columns within each Ellipses cell may be wrong 
%and they should be reordered. Plus some appending needs to be done to make
%the width of each matrix 8 or 9. Also, ditto about conversion from um to
%pixels.
