function [tile_array]=ManualStitchingCorrection(Prefix, ID)
% author: Gabriella Martini
% date created: 12/29/19
% date last modified: 12/31/19

%m - Choose the tile to move with the mouse
%t - Choose the tile to move using row and column indices
%right arrow - Move the selected tile to the right by one column
%left arrow - Move the selected tile to the left by one column 
%up arrow - Move the selected tile up by one row
%down arrow - Move the selected tile down by one row 
%> - Move the selected tile to the right by 5 columns
%< - Move the selected tile to the left by 5 columns
%+ - Move the selected tile up by 5 rows
%- - Move the selected tile down by 5 rows
%s - save
%x - close
%q - Save and close

%Make an overlay of the zoomed in and zoomed out real
%images as well as of a quickly segmented nuclear mask

%Close existing images
close all


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
%% Load folder information and stored TileArray Data 

% relevant folder info 
[SourcePath,FISHPath,DropboxFolder,MS2CodePath]=...
        DetermineLocalFolders(Prefix);
stitchingDataFolder = [DropboxFolder,filesep,Prefix,filesep,'FullEmbryoStitching'];

% load TileArray with stitching information 
if exist([stitchingDataFolder, filesep, ID, 'TileArray.mat'], 'file')
     load([stitchingDataFolder, filesep, ID, 'TileArray.mat']);
else
     error('No TileArray data stored. Seed a new TileArray using "NewTileArrayFromMetadata".')  
end
%% Load stitching information into memory 
rows = [tile_array.rows{:}];
cols = [tile_array.cols{:}];
tiles = tile_array.imgs;
NTiles = length(tiles);
heights = [];
widths = [];
gridrows = [];
gridcolumns = [];
for t=1:NTiles
    heights(t) = size(tiles{t}, 1);
    widths(t) = size(tiles{t},2);
    gridrows(t) = tile_array.grid_positions{t}(1);
    gridcolumns(t) = tile_array.grid_positions{t}(2);
end
rmaxs = rows+heights-1;
cmaxs = cols+widths-1;
grid_positions = tile_array.grid_positions;

%% Copy TileArray info into temporary values used to store intermediate stitiching values

temp_array = tile_array;
temprows = rows;
tempcols = cols;
temprmaxs = rmaxs;
tempcmaxs = cmaxs;
rshift = zeros(NTiles);
cshift = zeros(NTiles);

%% Generate a figure for viewing current stitching information 

TempFigure = figure;
axesTemp = axes(TempFigure);




cc=1;
TileSelected = 0;
Changed = 0;
%Overlay the zoom in and zoom out images
while (cc~='x' && cc~='q')
    imm2 = imstitchTileColor(temp_array);

    %imshow(APImage,DisplayRange)
    %imshow(imm2, 'Parent', axesTemp)
    hImage = image( imm2, 'Parent', axesTemp );
    axis(axesTemp, 'image')
    set(gca,'xtick',[],'ytick',[],'xlabel',[],'ylabel',[]);
    %axis image
    %axis off
    shiftinfo = {''};
    lcounter = 1;
    for t=1:NTiles
        if t == 1
            shiftinfo{lcounter} = [shiftinfo{lcounter},'$$T_{', num2str(grid_positions{t}(1)),...
                num2str(grid_positions{t}(2)), '}$$ shift: (',...
                num2str(rshift(t)), ',',num2str(cshift(t)),')'];
        else
            shiftinfo{lcounter} = [shiftinfo{lcounter},...
                ' $$T_{', num2str(grid_positions{t}(1)),...
                num2str(grid_positions{t}(2)), '}$$ shift: (',...
                num2str(rshift(t)), ',',num2str(cshift(t)),')'];
            if mod(t,4) == 0
                lcounter = lcounter + 1;
                shiftinfo{lcounter} = '';
            end
        end
    end
    if TileSelected == 0
        plottitle = {'No tile selected'};
    else
        plottitle = {['Selected Tile: (',...
            num2str(temp_array.grid_positions{TileSelected}(1)),...
            ',', num2str(temp_array.grid_positions{TileSelected}(2)),...
            ')']};
    end
    for pl=1:length(shiftinfo)
        plottitle{length(plottitle)+1} = shiftinfo{pl};
    end
  
    title(axesTemp, plottitle, 'Interpreter', 'latex');
    
    
    ct=waitforbuttonpress;
    cc=get(TempFigure,'currentcharacter');
    cc_value = double(cc);
    %disp(num2str(cc_value));
    cm=get(gca,'CurrentPoint');
    
    
    if (ct~=0)&&(cc=='m')        
        [ci,ri] = ginput(1);
        for t=1:NTiles
            if ((ri >= temprows(t)) && (ri <= temprmaxs(t)) && (ci >= tempcols(t)) && (ci <= tempcmaxs(t)))
                TileSelected = t;
                break
            end
        end
    elseif (ct~=0)&&(cc=='t')        %Select tile to move using keyboard inputs
        disp(['Enter a row number between ',...
            num2str(min(gridrows)),...
            ' and ', num2str(max(gridrows)),...
            '(press c to cancel): '])
        ct2  = waitforbuttonpress;
        gri=get(TempFigure,'currentcharacter');
        gri_value = double(gri);
        if gri_value == 99 % 'c'
            disp('Cancelled tile selection.')
            continue
        end
        while (gri_value - 48 > max(gridrows)) ||...
              (gri_value - 48 < min(gridrows))
            fprintf('%s not a valid row selection.', gri);
            fprintf('Enter a row number  between %d and %d (press "c" to cancel): \n',...
            min(gridrows), max(gridrows));
            ct2  =waitforbuttonpress;
            gri=get(TempFigure,'currentcharacter');
            gri_value = double(gri);
            if gri_value == 99
                disp('Cancelled tile selection.')
                continue
            end
        end
        if gri_value == 99
            break
        end
        fprintf('Row number %d selected. Enter a column number between %d and %d (press "c" to cancel): \n',...
            (gri_value-48), min(gridcolumns), max(gridcolumns))
        ct3 = waitforbuttonpress;
        gci=get(TempFigure,'currentcharacter');
        gci_value = double(gci);
        if gci_value == 99 % 'c' has a double value of 99
            continue
        end
        while (gci_value - 48 < min(gridcolumns)) ||...
              (gci_value - 48 > max(gridcolumns))
              
            %(sum(gridcolumns == (gci_value-48)) == 0
            fprintf('%s not a valid column selection.', gci);
            fprintf('Enter a column number  between %d and %d (press "c" to cancel): \n',...
            min(gridcolumns), max(gridcolumns))
           
            ct3  =waitforbuttonpress;
            gci=get(TempFigure,'currentcharacter');
            gci_value = double(gci);
            if gci_value == 99
                disp('Cancelled tile selection.')
                break
            end
        end
        if gci_value == 99
            continue
        end
        disp(['Column number ', num2str(gci_value-48), ' selected.']);
        for t=1:NTiles
            if ((gri_value-48 == grid_positions{t}(1)) && (gci_value-48 == grid_positions{t}(2)))
                TileSelected = t;
                break
            end
        end    
    elseif (ct~=0) && (cc_value == 30)
        %disp('left arrow');
        if TileSelected == 0
            disp('No tile selected. Select tile using "m" or "t" options.');
            continue
        end
        temprows(TileSelected) = temprows(TileSelected)-1;
        temprmaxs(TileSelected) = temprmaxs(TileSelected)-1;
        rshift(TileSelected) = rshift(TileSelected)-1;
    elseif (ct~=0) && (cc_value == 31)
        %disp('right arrow')
        if TileSelected == 0
            disp('No tile selected. Select tile using "m" or "t" options.');
            continue
        end
        temprows(TileSelected) = temprows(TileSelected)+1;
        temprmaxs(TileSelected) = temprmaxs(TileSelected)+1;
        rshift(TileSelected) = rshift(TileSelected)+1;
    elseif (ct~=0) && (cc_value == 28)
        %disp('up arrow')
        if TileSelected == 0
            disp('No tile selected. Select tile using "m" or "t" options.');
            continue
        end
        tempcols(TileSelected) = tempcols(TileSelected)-1;
        tempcmaxs(TileSelected) = tempcmaxs(TileSelected)-1;
        cshift(TileSelected) = cshift(TileSelected)-1;
    elseif (ct~=0) && (cc_value == 29)
        %disp('down arrow')
        if TileSelected == 0
            disp('No tile selected. Select tile using "m" or "t" options.');
            continue
        end
        tempcols(TileSelected) = tempcols(TileSelected)+1;
        tempcmaxs(TileSelected) = tempcmaxs(TileSelected)+1;
        cshift(TileSelected) = cshift(TileSelected)+1;
    elseif (ct~=0) && (cc_value == 43)
        %disp('-')
        if TileSelected == 0
            disp('No tile selected. Select tile using "m" or "t" options.');
            continue
        end
        temprows(TileSelected) = temprows(TileSelected)-5;
        temprmaxs(TileSelected) = temprmaxs(TileSelected)-5;
        rshift(TileSelected) = rshift(TileSelected)-5;
    elseif (ct~=0) && (cc_value == 45)
        %disp('+')
        if TileSelected == 0
            disp('No tile selected. Select tile using "m" or "t" options.');
            continue
        end
        temprows(TileSelected) = temprows(TileSelected)+5;
        temprmaxs(TileSelected) = temprmaxs(TileSelected)+5;
        rshift(TileSelected) = rshift(TileSelected)+5;
    elseif (ct~=0) && (cc_value == 62)
        %disp('>')
        if TileSelected == 0
            disp('No tile selected. Select tile using "m" or "t" options.');
            continue
        end
        tempcols(TileSelected) = tempcols(TileSelected)+5;
        tempcmaxs(TileSelected) = tempcmaxs(TileSelected)+5;
        cshift(TileSelected) = cshift(TileSelected)+5;
    elseif (ct~=0) && (cc_value == 60)
        %disp('<')
        if TileSelected == 0
            disp('No tile selected. Select tile using "m" or "t" options.');
            continue
        end
        tempcols(TileSelected) = tempcols(TileSelected)-5;
        tempcmaxs(TileSelected) = tempcmaxs(TileSelected)-5;
        cshift(TileSelected) = cshift(TileSelected)-5;
    elseif (ct~=0) && (cc=='s')
        %Save the information
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
    rmin = min(temprows);
    cmin = min(temprows);
    temprows = temprows-rmin + 1;
    tempcols = tempcols-cmin+1;
    temprmaxs = temprmaxs-rmin+1;
    tempcmaxs = tempcmaxs-cmin+1;
    for t=1:NTiles
        temp_array.rows{t} = temprows(t);
        temp_array.cols{t} = tempcols(t);
        
    end

    if all(temprows == rows) && all(tempcols == cols)
        Changed = 0;
    else 
        Changed = 1;
    end


end

if Changed == 1
    tile_array.prevrows{size(tile_array.prevrows, 2)+1} = tile_array.rows;
    tile_array.prevcols{size(tile_array.prevcols, 2)+1} = tile_array.cols;
    for t=1:NTiles
        tile_array.rows{t} = temprows(t);
        tile_array.cols{t} = tempcols(t);
        
    end
end

if cc=='q'
    %Save the information
    %Create the output folder if it doesn't exist 
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

close(TempFigure)
end
