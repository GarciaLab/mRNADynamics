
%Hardcoded variables for initial implementation 
clear all, close all
SourcePath = 'E:/Gabriella/LivemRNA\Data\RawDynamicsData';
Prefix = '2019-12-09-4xVasamEGFPVK22Homo-HbNbHomo-Anterior-27_5C-NC13NC14';
Date = '2019-12-09';
EmbryoName = '4xVasamEGFPVK22Homo-HbNbHomo-Anterior-27_5C-NC13NC14';


%% 

r_maxoffset = 200;
c_maxoffset = 200;
r_range = 2;
c_range = 2;
ID = 'Mid';
NIterations = 100;
tile_array = newTileArrayFromMetadata(Prefix, ID);
%% 

tile_array  = ManualStitchingCorrection(Prefix, tile_array, ID);
%% 
tile_array = improveTileStitching(Prefix,...
    tile_array, r_maxoffset, c_maxoffset, r_range, c_range, NIterations, ID);
%% 

imm2 = imstitchTileColor(tile_array);
imshow(imm2)

%% 
function tile_array = newTileArrayFromMetadata(Prefix, ID)
    if ~exist('Prefix')
        FolderTemp=uigetdir(DropboxFolder,'Choose folder with files to analyze');
        Dashes=strfind(FolderTemp,filesep);
        Prefix=FolderTemp((Dashes(end)+1):end);
    end    
    [SourcePath,FISHPath,DropboxFolder,MS2CodePath]=...
        DetermineLocalFolders(Prefix);
    %Find out the date it was taken
    Dashes=findstr(Prefix,'-');
    Date=Prefix(1:Dashes(3)-1);
    EmbryoName=Prefix(Dashes(3)+1:end);
    if ~isempty(strfind(lower(ID), 'mid'))
        filename = 'MidTile';
    elseif ~isempty(strfind(lower(ID), 'surf'))
        filename = 'SurfTile';
    else
        filename = ID;
    end
    LIFPath = [SourcePath, filesep, Date, filesep, EmbryoName,filesep,...
    'FullEmbryo\', filename,'.lif'];

    r = bfGetReader();
    % Decorate the reader with the Memoizer wrapper
    r = loci.formats.Memoizer(r);
    r.setId(LIFPath);
    LIFImages = bfopen(LIFPath);
    LIFMeta = LIFImages{:, 4};
    r.close();
    
    [Prefix, SkipFrames, ProjectionType, PreferredFileNameForTest, keepTifs,...
    generateTifStacks, nuclearGUI, skipExtraction, rootFolder, zslicesPadding,...
    lowbit] = exportDataForLivemRNA_processInputParameters(Prefix);
    [rawDataPath, ~, DropboxFolder, ~, PreProcPath, rawDataFolder, Prefix, ExperimentType, Channel1, Channel2, ~,...
        Channel3] = readMovieDatabase(Prefix,'rootFolder', rootFolder);
    [NSeries, NFrames, NSlices, NPlanes, NChannels, Frame_Times] = getFrames(LIFMeta);
    if sum(NFrames)~=0
        [Frame_Times, First_Time] = obtainFrameTimes(XMLFolder, seriesPropertiesXML, NSeries, NFrames, NSlices, NChannels);
        [InitialStackTime, zPosition] = getFirstSliceTimestamp(NSlices, NSeries, NPlanes, NChannels, Frame_Times, XMLFolder, seriesXML);
    else
        InitialStackTime = [];
        zPosition = [];
    end
    framesIndex  = 1;
    NTiles = NSeries;
    FrameInfo = recordFrameInfo(NFrames, NSlices, InitialStackTime, LIFMeta, zPosition);

    [coatChannel, histoneChannel, fiducialChannel, inputProteinChannel, FrameInfo] =...
    LIFExportMode_interpretChannels(ExperimentType, Channel1, Channel2, Channel3, FrameInfo);

    %Create the output folder
    mkdir([DropboxFolder,filesep,Prefix,filesep,'StitchedEmbryoImages'])



    framesIndex = 1;
    tiles = {};
    for n=1:NTiles
        ImageSlices = generateHisSlicesTileScan(LIFImages, NSlices(n),...
            NChannels, fiducialChannel, framesIndex, n);
        temp = max(ImageSlices, [], 3);
        tiles{n} = temp;
    end


    PixelSize = double(LIFMeta.getPixelsPhysicalSizeX(1).value);% units: microns
    PixelSize_m = double(PixelSize)*10^(-6);
    ypos = [];
    xpos = [];
    for i=0:(NTiles-1)
        xpos(length(xpos)+1) = -1*double(LIFMeta.getPlanePositionX(i,0).value);
        ypos(length(ypos)+1) = double(LIFMeta.getPlanePositionY(i,0).value);
    end
    uxpos = sort(unique(xpos), 'ascend');
    uypos = sort(unique(ypos), 'ascend');
    xdim = length(uxpos); ydim = length(uypos);
    dx = round(abs((uxpos(xdim)-uxpos(xdim-1))/PixelSize_m), 0);
    dy = round(abs((uypos(ydim)-uypos(ydim-1))/PixelSize_m), 0);
    tile_array.rows = {};
    tile_array.cols = {};
    tile_array.grid_positions = {};
    tile_array.imgs = {};
    sigma = .6/PixelSize;
    for i=1:NTiles
        ri = round((xpos(i)-uxpos(1))/PixelSize_m, 0)+1;
        ci = round((ypos(i)-uypos(1))/PixelSize_m, 0)+1;
        gri = find(xpos(i) == uxpos);
        gci = find(ypos(i) == uypos);
        tile_array.rows{i} = ri;
        tile_array.cols{i} = ci;
        tile_array.grid_positions{i} = [gri, gci];
        tile_array.imgs{i} = imgaussfilt(tiles{i}, sigma);
    end

    tile_array.prevrows = {};
    tile_array.prevcols = {};
    
       %Save the information
    saveVars = {};
    saveVars = [saveVars, 'tile_array'];
    save([DropboxFolder,filesep,Prefix,filesep,ID,'TileStitch.mat'],saveVars{:});
    imm2 = imstitchTile(tile_array);
    mkdir([DropboxFolder,filesep,Prefix,filesep,'APDetection'])
    imwrite(imm2,[DropboxFolder,filesep,Prefix,filesep,'StitchedEmbryoImages',filesep,'FullEmbryo',ID,'.tif'],'compression','none');
    imwrite(imm2,[DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'FullEmbryo',ID,'.tif'],'compression','none');
end


function tile_array = improveTileStitching(Prefix,...
    tile_array, r_maxoffset, c_maxoffset, r_range, c_range, NIterations, ID)
    [SourcePath,FISHPath,DropboxFolder,MS2CodePath]=...
        DetermineLocalFolders(Prefix);
    imm2 = imstitchTile(tile_array);
    imshow(imm2)
    % Start by moving tile1 relative to tiles 2, 3, and 4
    % allow x and y positions that are within 100 of the current value, without
    % failing to overlap neighbors and sstaying within prescribed limits:
    numTiles = length(tile_array.imgs);
    for iter=1:NIterations
        for tA_ind=1:numTiles
            tileA = tile_array.imgs{tA_ind};
            [hA, wA] = size(tileA);
            tile_array.prevrows{size(tile_array.prevrows, 2)+1} = tile_array.rows;
            tile_array.prevcols{size(tile_array.prevcols, 2)+1} = tile_array.cols;
            tAr = tile_array.rows{tA_ind}; tAc = tile_array.cols{tA_ind};
            [tAr_min, tAr_max, tAc_min, tAc_max ] = ...
                getRowColLimits(tile_array, tA_ind, r_maxoffset, c_maxoffset, r_range, c_range);

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
            %imm2 = imstitchTile(tile_array);
            %figure(1)
            %imshow(imm2)
            %filename = ['E:/Gabriella/LivemRNA\Data\DynamicsResults\',...
             %           Prefix,'\TilingImages\StitchedMid.png'];
            %title(['Iteration: ', num2str((iter-1)*numTiles+tA_ind)])
            %saveas(gcf, filename)


        end
        numIters = length(tile_array.prevrows);
        oldrows = [tile_array.prevrows{numIters-(numTiles-1)}{:}];
        oldcols = [tile_array.prevcols{numIters-(numTiles-1)}{:}];
        newrows = [tile_array.rows{:}];
        newcols = [tile_array.cols{:}];
        if mod(iter, 10) == 0
            %fprintf('%d ', i); 
            fprintf('%d/%d', [iter, NIterations]);
            %fprintf([num2str(iter),'/',num2str(NIterations)]);
        else
            fprintf('.');
        end
        
        if (isequal(oldrows, newrows) && isequal(oldcols, newcols))
            imm = imstitchTile(tile_array);
            figure            
            imshow(imm)
            %mkdir(
            filename = ['E:/Gabriella/LivemRNA\Data\DynamicsResults\',...
                        Prefix,'\TilingImages\StitchedMid.png'];
            title(['Iteration: ', num2str(iter)])
            saveas(gcf, filename)
            fprintf('Done!\n');
            break
        end
        
    end
    if ~(isequal(oldrows, newrows) || ~isequal(oldcols, newcols))
       disp('Tile Positions failed to stabilize.');
    end
    %Save the information
    saveVars = {};
    saveVars = [saveVars, 'tile_array'];
    save([DropboxFolder,filesep,Prefix,filesep,ID,'TileStitch.mat'],saveVars{:});
    imm2 = imstitchTile(tile_array);
    mkdir([DropboxFolder,filesep,Prefix,filesep,'APDetection'])
    imwrite(imm2,[DropboxFolder,filesep,Prefix,filesep,'StitchedEmbryoImages',filesep,'FullEmbryo',ID,'.tif'],'compression','none');
    imwrite(imm2,[DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'FullEmbryo',ID,'.tif'],'compression','none');
end

%% 


function [r_min, r_max, c_min, c_max] = getRowColLimits(tile_array, tA_ind, r_maxoffset, c_maxoffset, r_range, c_range)
    r_min = tile_array.rows{tA_ind} - r_range; 
    r_max = tile_array.rows{tA_ind} + r_range; 
    c_min = tile_array.cols{tA_ind} - c_range; 
    c_max = tile_array.cols{tA_ind} + c_range; 
    numImages = length(tile_array.imgs);
    if isequal(tile_array.grid_positions{tA_ind}, [1,1])
        for n=1:numImages
            if isequal(tile_array.grid_positions{n}, [1,2])
                 [c_min, c_max] = ...
                     correctColLimitsRightTile(tile_array, tA_ind,...
                     n, c_min, c_max, c_maxoffset);
            elseif isequal(tile_array.grid_positions{n}, [2,1])
                [r_min, r_max] = ...
                    correctRowLimitsLowerTile(tile_array, tA_ind, n,...
                    r_min, r_max, r_maxoffset);
            end    
        end
    elseif isequal(tile_array.grid_positions{tA_ind}, [1,2])
        for n=1:numImages
            if isequal(tile_array.grid_positions{n}, [1,1])
                [c_min, c_max] = ...
                     correctColLimitsLeftTile(tile_array, tA_ind,...
                     n, c_min, c_max, c_maxoffset);
            elseif isequal(tile_array.grid_positions{n}, [1,3])
                [c_min, c_max] = ...
                     correctColLimitsRightTile(tile_array, tA_ind,...
                     n, c_min, c_max, c_maxoffset);
            elseif isequal(tile_array.grid_positions{n}, [2,2])
                [r_min, r_max] = ...
                    correctRowLimitsLowerTile(tile_array, tA_ind, n,...
                    r_min, r_max, r_maxoffset);
            end    
        end   
    elseif isequal(tile_array.grid_positions{tA_ind}, [1,3])
        for n=1:numImages
            if isequal(tile_array.grid_positions{n}, [1,2])
                [c_min, c_max] = ...
                     correctColLimitsLeftTile(tile_array, tA_ind,...
                     n, c_min, c_max, c_maxoffset);
            elseif isequal(tile_array.grid_positions{n}, [2,3])
                [r_min, r_max] = ...
                    correctRowLimitsLowerTile(tile_array, tA_ind, n,...
                    r_min, r_max, r_maxoffset);
            end    
        end
    elseif isequal(tile_array.grid_positions{tA_ind}, [2,1])
        for n=1:numImages
            if isequal(tile_array.grid_positions{n}, [2,2])
                [c_min, c_max] = ...
                     correctColLimitsRightTile(tile_array, tA_ind,...
                     n, c_min, c_max, c_maxoffset);
            elseif isequal(tile_array.grid_positions{n}, [3,1])
                [r_min, r_max] = ...
                    correctRowLimitsLowerTile(tile_array, tA_ind, n,...
                    r_min, r_max, r_maxoffset);
            elseif isequal(tile_array.grid_positions{n}, [1,1])
                [r_min, r_max] = ...
                    correctRowLimitsUpperTile(tile_array, tA_ind, n,...
                    r_min, r_max, r_maxoffset);
            end    
        end
    elseif isequal(tile_array.grid_positions{tA_ind}, [2,2])
        for n=1:numImages
            if isequal(tile_array.grid_positions{n}, [2,3])
                [c_min, c_max] = ...
                     correctColLimitsRightTile(tile_array, tA_ind,...
                     n, c_min, c_max, c_maxoffset);
            elseif isequal(tile_array.grid_positions{n}, [2,1])
                [c_min, c_max] = ...
                     correctColLimitsLeftTile(tile_array, tA_ind,...
                     n, c_min, c_max, c_maxoffset);
            elseif isequal(tile_array.grid_positions{n}, [3,2])
                [r_min, r_max] = ...
                    correctRowLimitsLowerTile(tile_array, tA_ind, n,...
                    r_min, r_max, r_maxoffset);
            elseif isequal(tile_array.grid_positions{n}, [1,2])
                [r_min, r_max] = ...
                    correctRowLimitsUpperTile(tile_array, tA_ind, n,...
                    r_min, r_max, r_maxoffset);
            end    
        end
    elseif isequal(tile_array.grid_positions{tA_ind}, [2,3])
        for n=1:numImages
            if isequal(tile_array.grid_positions{n}, [2,2])
                [c_min, c_max] = ...
                     correctColLimitsLeftTile(tile_array, tA_ind,...
                     n, c_min, c_max, c_maxoffset);
            elseif isequal(tile_array.grid_positions{n}, [3,3])
                [r_min, r_max] = ...
                    correctRowLimitsLowerTile(tile_array, tA_ind, n,...
                    r_min, r_max, r_maxoffset);
            elseif isequal(tile_array.grid_positions{n}, [1,3])
                [r_min, r_max] = ...
                    correctRowLimitsUpperTile(tile_array, tA_ind, n,...
                    r_min, r_max, r_maxoffset);
            end    
        end
    elseif isequal(tile_array.grid_positions{tA_ind}, [3,1])
        for n=1:numImages
            if isequal(tile_array.grid_positions{n}, [3,2])
                [c_min, c_max] = ...
                     correctColLimitsRightTile(tile_array, tA_ind,...
                     n, c_min, c_max, c_maxoffset);
            elseif isequal(tile_array.grid_positions{n}, [2,1])
                [r_min, r_max] = ...
                    correctRowLimitsUpperTile(tile_array, tA_ind, n,...
                    r_min, r_max, r_maxoffset);
            end    
        end
    elseif isequal(tile_array.grid_positions{tA_ind}, [3,2])
        for n=1:numImages
            if isequal(tile_array.grid_positions{n}, [3,3])
                [c_min, c_max] = ...
                     correctColLimitsRightTile(tile_array, tA_ind,...
                     n, c_min, c_max, c_maxoffset);
            elseif isequal(tile_array.grid_positions{n}, [3,1])
                [c_min, c_max] = ...
                     correctColLimitsLeftTile(tile_array, tA_ind,...
                     n, c_min, c_max, c_maxoffset);
            elseif isequal(tile_array.grid_positions{n}, [2,2])
                [r_min, r_max] = ...
                    correctRowLimitsUpperTile(tile_array, tA_ind, n,...
                    r_min, r_max, r_maxoffset);
            end    
        end
    elseif isequal(tile_array.grid_positions{tA_ind}, [3,3])
        for n=1:numImages
            if isequal(tile_array.grid_positions{n}, [3,2])
                [c_min, c_max] = ...
                     correctColLimitsLeftTile(tile_array, tA_ind,...
                     n, c_min, c_max, c_maxoffset);
            elseif isequal(tile_array.grid_positions{n}, [2,3])
                [r_min, r_max] = ...
                    correctRowLimitsUpperTile(tile_array, tA_ind, n,...
                    r_min, r_max, r_maxoffset);
            end    
        end        
    else
        error('need to write other cases')
    end

end

function [r_min, r_max] = correctRowLimitsLowerTile(tile_array, tA_ind, n, r_min, r_max, r_maxoffset)
    hA = size(tile_array.imgs{tA_ind}, 1);
    r_min = max([r_min, (tile_array.rows{n} - hA +1)]);
    r_max = min([r_max, (tile_array.rows{n} - hA+r_maxoffset)]);
end

function [r_min, r_max] = correctRowLimitsUpperTile(tile_array, tA_ind, n, r_min, r_max, r_maxoffset)
    h = size(tile_array.imgs{n}, 1);
    r_min = max([r_min, (tile_array.rows{n} + h - r_maxoffset)]);
    r_max = min([r_max, (tile_array.rows{n} + h-1)]);
end

function [c_min, c_max] = correctColLimitsRightTile(tile_array, tA_ind, n, c_min, c_max, c_maxoffset)
    wA = size(tile_array.imgs{tA_ind}, 2);
    c_min = max([c_min, (tile_array.cols{n} - wA +1)]);
    c_max = min([c_max, (tile_array.cols{n} - wA+c_maxoffset)]);
end

function [c_min, c_max] = correctColLimitsLeftTile(tile_array, tA_ind, n, c_min, c_max, c_maxoffset)
    w = size(tile_array.imgs{n}, 2);
    c_min = max([c_min, (tile_array.cols{n} + w - c_maxoffset)]);
    c_max = min([c_max, (tile_array.cols{n} + w-1)]);
end

% 
%% 







