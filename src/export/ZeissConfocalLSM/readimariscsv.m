function [schnitzcells, Ellipses] = readimariscsv(positionFile, pixelXSize, pixelYSize, pixelZSize)

    arguments
        positionFile char
        pixelXSize double
        pixelYSize double
        pixelZSize double
    end

    pixelSizes.cenx = pixelXSize;
    pixelSizes.ceny = pixelYSize;
    pixelSizes.cenz = pixelZSize;

    % first 3 rows in imaris CSV are blank and/or file header
    % we want to start parsing for column names at row 3 (zero-based, so it's row 4 of the file)
    opts = detectImportOptions(positionFile);
    opts.VariableNamesLine = 3;

    T = readtable(positionFile, opts);

    % drop needless features
    T = removevars(T, {'Unit', 'Category', 'Collection'});

    % groupby nucleus ID and aggregate into a list.
    % @(x) {x} function handle has the effect of returning the same element,
    % so the effect is that a list is created. Other optionw when grouping by are functions
    % like 'sum', 'mean', 'average', etc. But in this case, we wan the actual aggregation of values.
    T_groupedByID = groupsummary(T, 'TrackID', @(x) {x});

    % colums created by group by with a custom function handle will be named fun1_etc
    % let's rename them to match column names used in schnitzcells mat.
    T_groupedByID = renamevars(T_groupedByID,...
        ["fun1_PositionX", "fun1_PositionY", "fun1_PositionZ", "fun1_Time", "fun1_ID"],...
        ["cenx", "ceny", "cenz", "frames", "ImarisID"]);

    % This removes any nuclei with a missing "TrackID" value on the original file.
    % We think these are all stray nuclei that weren't tracked through multiple frames,
    % and it's ok to drop them.
    T_groupedByID = rmmissing(T_groupedByID);

    % sort all columns of type cell (the list results from the groupby) by frame number
    % just in case frames are not in proper order
    for row = 1:height(T_groupedByID)
        % take the frames array from the current row (for example [2,1])
        timevec = T_groupedByID(row, 'frames').frames{1};
        
        % iterate across columns of current row, looking for values of type cell
        for col = 1:width(T_groupedByID)

            val = T_groupedByID.(col)(row);
            if iscell(val)
                % create a temporary matrix with frame numbers as first column
                tempMatrix = [timevec, cell2mat(T_groupedByID.(col)(row))];
                
                % sort by first column (frame numbers)
                sortedMatrix = sortrows(tempMatrix, 1);
                
                % rewrite current col/row value with values sorted by frame number
                % (but without the actual frame numbers)
                T_groupedByID.(col)(row) = {sortedMatrix(:, 2)};

                columnName = T_groupedByID.Properties.VariableNames(col);
                columnName = columnName{1};

                if ismember(columnName, ['cenx', 'ceny', 'cenz'])
                
                    try
                        ratio = extractfield(pixelSizes, columnName);
                    catch
                        % if no pixel resolution info available for this column name,
                        % we leave the value unchanged
                        ratio = 1;
                    end

                    T_groupedByID.(col)(row) = convertSizeToPixels(T_groupedByID.(col)(row), ratio);
                end
            end
        end
        
    end

    % convert all positions from um to voxels, according to voxel size metadata


    % convert table to struct, this should basically match schnitzcells as if they were produced by TrackNuclei. 
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
end
