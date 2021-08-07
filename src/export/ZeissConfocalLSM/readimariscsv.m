function [schnitzcells, Ellipses, lastFrameNumber] = readimariscsv(frameNumberOffset, imarisStatisticsFolder, positionFile, pixelXSize, pixelYSize, pixelZSize)

    arguments
        frameNumberOffset double
        imarisStatisticsFolder char
        positionFile char
        pixelXSize double
        pixelYSize double
        pixelZSize double
    end

    fprintf('About to parse Imaris CSV data with frame offset %d\n', frameNumberOffset);

    pixelSizes.cenx = pixelXSize;
    pixelSizes.ceny = pixelYSize;
    pixelSizes.cenz = pixelZSize;


    % first 3 rows in imaris CSV are blank and/or file header
    % we want to start parsing for column names at row 3 (zero-based, so it's row 4 of the file)
    opts = detectImportOptions(positionFile);
    opts.VariableNamesLine = 3;

    T = readtable(positionFile, opts);
    
    % offset to use overall frame numbering.
    % See comment in the code that calls this function.
    T.Time = T.Time + frameNumberOffset;
    % after offset, save what is our current last frame number
    % this will be the offset for the next file
    lastFrameNumber = max(T.Time);

    % drop needless features
    T = removevars(T, {'Unit', 'Category', 'Collection'});

    % groupby nucleus ID and aggregate into a list.
    % @(x) {x} function handle has the effect of returning the same element,
    % so the effect is that a list is created. Other options when grouping-by are functions
    % like 'sum', 'mean', 'average', etc. But in this case, we want the actual aggregation of values.
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

                    % convert all positions from um to voxels, according to voxel size metadata
                    T_groupedByID.(col)(row) = convertSizeToPixels(T_groupedByID.(col)(row), ratio);
                end
            end
        end
        
    end

    
    % add legacy columns for backwards compatibility
    % also add non critical values as zeros or empty arrays
    tableHeight = height(T_groupedByID);
    emptyArrays = cell(tableHeight,1);
    zeroValues = zeros(tableHeight, 1);
    
    T_groupedByID.P = emptyArrays;
    T_groupedByID.E = emptyArrays;
    T_groupedByID.D = emptyArrays;
    T_groupedByID.len = emptyArrays;
    T_groupedByID.cellno = emptyArrays;
    
    T_groupedByID.AlreadyUsed = zeroValues;
    T_groupedByID.ExtendedIntoFutureAlready = zeroValues;
    T_groupedByID.StitchedTo = emptyArrays;
    T_groupedByID.StitchedFrom = emptyArrays;
    T_groupedByID.cycle = zeroValues;
    T_groupedByID.timeSinceAnaphase = emptyArrays;
    T_groupedByID.catchedAnaphase = zeroValues;
    
    % let's put the columns in the same order TrackNuclei puts them,
    % just in case.
    T_groupedByID = movevars(T_groupedByID,{'P', 'E', 'D', 'frames', 'cenx', 'ceny', 'len', 'cellno', 'AlreadyUsed', 'ExtendedIntoFutureAlready', 'StitchedTo', 'StitchedFrom', 'cycle', 'timeSinceAnaphase', 'catchedAnaphase'},'Before',1);
    
    % Now, we'll make Ellipses.
    % We start by reading all ellipsoid raw info from imaris CSV files
    [A, B, C] = readImarisEllipsoidTables(imarisStatisticsFolder);

    % combine (join) all three tables of ellipsoid axis lenghts A, B andC into one table for convenience
    TJoin = join(A, B);
    TJoin = join(TJoin, C);
    
    % we need to offset time data in these tables as well
    % so it matches the offset we applied before to the Positions table
    TJoin.Time = TJoin.Time + frameNumberOffset;

    % remove rows with missing column values
    T_withoutMissing = rmmissing(T);
    
    % now join the ellipsoid info to the positions table, to gather
    % all info into one big table.
    % we join all data by time, trackID, and ID columns
    TJoin_all = join(T_withoutMissing, TJoin, 'LeftKeys', [4,5,6], 'RightKeys', [2,3,4]);
    
    % now we group all this info by the Time column 
    T_groupedByTime = groupsummary(TJoin_all, 'Time', @(x) {x});

    T_groupedByTime = renamevars(T_groupedByTime, ["fun1_PositionX", "fun1_PositionY", "fun1_PositionZ", "fun1_TrackID", "fun1_ID", "fun1_EllipsoidAxisLengthA", "fun1_EllipsoidAxisLengthB", "fun1_EllipsoidAxisLengthC"],["cenx", "ceny", "cenz", "TrackID", "ID", "axisA", "axisB", "axisC"]);


    % We need to place the resulting data into the exact format of Ellipses
    Ellipses = cell(height(T_groupedByTime), 1);
    for frame = 1:height(T_groupedByTime)
        val = T_groupedByTime.cenx(frame);
        vec = val{1};
        Ellipses{frame} = zeros(numel(vec), 1);

        val = T_groupedByTime.cenx(frame);
        val = convertSizeToPixels(val, pixelSizes.cenx);
        vec = val{1};
        Ellipses{frame}(:, 1) = vec;

        val = T_groupedByTime.ceny(frame);
        val = convertSizeToPixels(val, pixelSizes.ceny);
        vec = val{1};
        Ellipses{frame}(:, 2) = vec;

        val = T_groupedByTime.axisA(frame);
        val = convertSizeToPixels(val, pixelSizes.cenx);
        vec = val{1};
        Ellipses{frame}(:, 3) = vec;


        val = T_groupedByTime.axisB(frame);
        val = convertSizeToPixels(val, pixelSizes.ceny);
        vec = val{1};
        Ellipses{frame}(:, 4) = vec;

        % column number 5 is the ellipses orientation angle, used in notEllipseCPT.m
        % to calculate angle "ang" by doing "pi - column5"
        % TO-DO grab orientation information from Imaris, if available.
        Ellipses{frame}(:, 5) = zeros(numel(vec), 1);
        Ellipses{frame}(:, 6) = zeros(numel(vec), 1);
        Ellipses{frame}(:, 7) = zeros(numel(vec), 1);
        Ellipses{frame}(:, 8) = zeros(numel(vec), 1);

        % we register here in Ellipses the indexes corresponding to where the nuclei are in schnitzcells structure
        TrackID = T_groupedByTime.TrackID(frame);
        [~,IndexInSchnitzcells] = ismember(TrackID{1}, T_groupedByID.TrackID);
        vec = IndexInSchnitzcells;
        Ellipses{frame}(:, 9) = vec;

        % and here we do the oposite, we register in schnitzcells where the nuclei are in ellipses structure
        TrackIDFromSchnitzcells = T_groupedByID.TrackID;
        [~,IndexInEllipses] = ismember(TrackIDFromSchnitzcells, TrackID{1});
        vec = IndexInEllipses;
        %C = arrayfun(@(a, b) [a, b], T_groupedByID.cellno, vec, 'UniformOutput', false);
        C = cellfun(@concatArrayElement, T_groupedByID.cellno, num2cell(vec), 'UniformOutput', false);
        T_groupedByID.cellno = C;

        
        % we add z as a non-standard 10th column to Ellipses for now
        % JP: still have to check backwards compatibility of doing this
        val = T_groupedByTime.cenz(frame);
        val = convertSizeToPixels(val, pixelSizes.cenz);
        vec = val{1};
        Ellipses{frame}(:, 10) = vec;
    end

    
    % convert table to struct, this should basically match schnitzcells as if they were produced by TrackNuclei. 
    % we already used TrackID column for Ellipses generation, we can now
    % remove it and other unwanted columns before generating the struct
    T_groupedByID = removevars(T_groupedByID, {'TrackID', 'GroupCount', 'ImarisID'});
    schnitzcells = table2struct(T_groupedByID);
end
