function [stitchOrder]=ManualStitchOrder(Prefix, ID)
% author: Gabriella Martini
% date created: 7/26/20
% date last modified: 7/26/20

%m - Choose the tile to move with the mouse
%t - Choose the tile to move using row and column indices
%right arrow - Move the selected tile to the right by one column
%left arrow - Move the selected tile to the left by one column 
%up arrow - Move the selected tile up by one row
%down arrow - Move the selected tile down by one row 
%> - Move the selected tile to the right by one column
%< - Move the selected tile to the left by one column
%enter - Move to next tile
%q - quit without saving

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

%% Generate a figure for viewing current stitching information 

TempFigure = figure;
axesTemp = axes(TempFigure);




cc=1;
TileSelected = 1;
gr = gridrows(TileSelected); gc = gridcolumns(TileSelected);
counter = 0;
stitchOrder = [];
plottitle= 'Use arrows to move selected tile. Press enter to choose selected tile.';
order_labels = [];

%Overlay the zoom in and zoom out images
while (cc~='x' && cc~='q')
    hold off
    rmin = rows(TileSelected);
    height = rmaxs(TileSelected)-rmin+1;
    cmin = cols(TileSelected);
    width = cmaxs(TileSelected)-cmin+1;
    displayNumbers = [];
    imm2 = imstitchTile(tile_array);

    %imshow(APImage,DisplayRange)
    %imshow(imm2, 'Parent', axesTemp)
    hImage = imagesc( imm2, 'Parent', axesTemp );
    axis(axesTemp, 'image')
    set(gca,'xtick',[],'ytick',[],'xlabel',[],'ylabel',[]);
    hold on
    curr_rect = rectangle('Position',[cmin,rmin,width,height],...
      'EdgeColor', 'r',...
      'LineWidth', 3,...
      'LineStyle','-');
    if length(order_labels) > 0
        for i=1:length(order_labels)
            label_info = order_labels{i};
            text(label_info{1}, label_info{2}, label_info{3}, 'Color', 'yellow', 'FontSize', 14)
        end
    end
    

    %axis image
    %axis off
    
    title(axesTemp, plottitle, 'Interpreter', 'latex');

    if counter >= 10
        break
    end
    ct=waitforbuttonpress;
    cc=get(TempFigure,'currentcharacter');
    cc_value = double(cc);
    %disp(num2str(cc_value));
    cm=get(gca,'CurrentPoint');
    
    
    if (ct~=0) && (cc_value == 29)
        %disp('right arrow')
        for tn=1:NTiles
            if isequal(tile_array.grid_positions{tn}, [gr, gc+1])
                TileSelected = tn;
                gc = gc+1;
                break
            end
        end  
    elseif (ct~=0) && (cc_value == 28)
        %disp('left arrow')
        for tn=1:NTiles
            if isequal(tile_array.grid_positions{tn}, [gr, gc-1])
                TileSelected = tn;
                gc = gc-1;
                break
            end
        end
    elseif (ct~=0) && (cc_value == 30)
        %disp('up arrow')
        for tn=1:NTiles
            if isequal(tile_array.grid_positions{tn}, [gr-1, gc])
                TileSelected = tn;
                gr = gr-1;
                break
            end
        end
    elseif (ct~=0) && (cc_value == 31)
        %disp('down arrow')
        for tn=1:NTiles
            if isequal(tile_array.grid_positions{tn}, [gr+1, gc])
                TileSelected = tn;
                gr = gr+1;
                break
            end
        end
    elseif (ct~=0) && (cc_value == 13)
        %disp('enter')
        if ~ismember(TileSelected, stitchOrder)
            stitchOrder(length(stitchOrder)+1) = TileSelected; 
            counter = counter + 1;
            order_labels{length(order_labels)+1} ={cmin+width*.8-1,rmin+height*.8-1,num2str(counter)};
            TileSelected = 1;
            gr = gridrows(TileSelected); gc = gridcolumns(TileSelected);
        end

    end

end



close(TempFigure)
end
