function StoreSpotInfo(Prefix)
liveExperiment = LiveExperiment(Prefix);
Spots = getSpots(liveExperiment);
snippet_size = 13;
xDim = liveExperiment.xDim;
yDim = liveExperiment.yDim;
NumFrames = length(Spots);
TifList = dir([liveExperiment.preFolder, filesep, '*ch01.tif']);
ProbList = dir([liveExperiment.procFolder,'dogs',  filesep, 'prob*']);
MovieCells = cell(1, NumFrames);
ProbCells = cell(1, NumFrames);
NumSlices = 0;
for i = 1:NumFrames
    MovieCells{i} = imreadStack([TifList(i).folder, filesep,  TifList(i).name]);
    ProbCells{i} = imreadStack([ProbList(i).folder, filesep,  ProbList(i).name]);
    NumSlices = NumSlices + size(MovieCells{i}, 3);
end
%%



%%
FrameIndex = 1;
SpotIndex = 1;
if isfile([liveExperiment.resultsFolder, filesep, 'StoreSpotInfo.mat'])
    load([liveExperiment.resultsFolder, filesep, 'StoreSpotInfo.mat']);
else
    GoodIndices = cell(1, length(Spots));
    BadIndices = cell(1, length(Spots));
    MultiIndices = cell(1, length(Spots));
end
%%
theta = 0 : (2 * pi / 10000) : (2 * pi);
pline_x = snippet_size/2 * cos(theta) + snippet_size+1;
pline_y = snippet_size/2  * sin(theta) + snippet_size+1;

close all
TempFigure = figure;
set(TempFigure, 'units', 'normalized', 'position', [0.2, 0.2, .2, .25]);
axesTemp = axes(TempFigure);
axis off
CurrentZ = Spots(FrameIndex).Fits(SpotIndex).z(1);
plot_title= {strrep(TifList(FrameIndex).name, '_', '\_'), ['z: ', num2str(CurrentZ)]};
SnippetImage2 = zeros(2*snippet_size+1,2*snippet_size+1, 'double');

xmin = Spots(FrameIndex).Fits(SpotIndex).xDoG-snippet_size;
xmax = Spots(FrameIndex).Fits(SpotIndex).xDoG+snippet_size;
ymin = Spots(FrameIndex).Fits(SpotIndex).yDoG-snippet_size;
ymax = Spots(FrameIndex).Fits(SpotIndex).yDoG+snippet_size;
if xmin> 0 & ymin > 0 & xmax <= xDim & ymax <= yDim
    SnippetImage2 = MovieCells{FrameIndex}(ymin:ymax,xmin:xmax, CurrentZ);
else
    if xmin <= 0
        xSnip = 1-(xmin);
    else
        xSnip =1;
    end
    if ymin <= 0
        ySnip = 1-(ymin);
    else
        ySnip=1;
    end
    xmaxSnip = xSnip + (xmax-xmin);
    ymaxSnip = ySnip + (ymax-ymin);
    SnippetImage2(ySnip:ymaxSnip, xSnip:xmaxSnip) = MovieCells{FrameIndex}(ymin:ymax,xmin:xmax, CurrentZ);
end
imagesc(axesTemp, SnippetImage2);
colormap('gray') 
hold on 
plot(pline_x, pline_y, '.r');
axis off
title(plot_title)

if isempty(BadIndices{FrameIndex})
    BadIndices{FrameIndex} = [];
end

if isempty(GoodIndices{FrameIndex})
    GoodIndices{FrameIndex} = [];
end

if isempty(MultiIndices{FrameIndex})
    MultiIndices{FrameIndex} = [];
end
if ismember(SpotIndex, BadIndices{FrameIndex})
     set(TempFigure,'Color','r')
elseif ismember(SpotIndex, GoodIndices{FrameIndex})
     set(TempFigure,'Color','g')
elseif ismember(SpotIndex, MultiIndices{FrameIndex})
     set(TempFigure,'Color','b')
else
    set(TempFigure,'Color','y')
end

hold off

%%
currentCharacter = 0;
set(0, 'CurrentFigure', TempFigure)
    
    
while (currentCharacter~='x')
    CurrentZ = Spots(FrameIndex).Fits(SpotIndex).z(1);
    plot_title= {strrep(TifList(FrameIndex).name, '_', '\_'), ['z: ', num2str(CurrentZ)]};
    %Load subsequent images
    
    if isempty(BadIndices{FrameIndex})
        BadIndices{FrameIndex} = [];
    end
    
    if isempty(GoodIndices{FrameIndex})
        GoodIndices{FrameIndex} = [];
    end
    
    if isempty(MultiIndices{FrameIndex})
        MultiIndices{FrameIndex} = [];
    end
   
    SnippetImage = zeros(2*snippet_size+1,2*snippet_size+1, 'double');
    xmin = Spots(FrameIndex).Fits(SpotIndex).xDoG-snippet_size;
    xmax = Spots(FrameIndex).Fits(SpotIndex).xDoG+snippet_size;
    ymin = Spots(FrameIndex).Fits(SpotIndex).yDoG-snippet_size;
    ymax = Spots(FrameIndex).Fits(SpotIndex).yDoG+snippet_size;
    if xmin> 0 & ymin > 0 & xmax <= xDim & ymax <= yDim
        SnippetImage = MovieCells{FrameIndex}(ymin:ymax,xmin:xmax,CurrentZ);%.*double(SpotMat(xmin:xmax,ymin:ymax,ReorderedIndex));
    else
        if xmin <= 0
            xSnip = 1-(xmin);
            xmin = 1;
        else
            xSnip =1;
        end
        if ymin <= 0
            ySnip = 1-(ymin);
            ymin = 1;
        else
            ySnip=1;
        end
        if xmax > xDim
            xmax = xDim;
        end
        if ymax > yDim
            ymax = yDim;
        end
        xmaxSnip = xSnip + (xmax-xmin);
        ymaxSnip = ySnip + (ymax-ymin);
        SnippetImage(ySnip:ymaxSnip,xSnip:xmaxSnip) = MovieCells{FrameIndex}(ymin:ymax,xmin:xmax,CurrentZ);%.*double(SpotMat(xmin:xmax,ymin:ymax,ReorderedIndex));
    end
    clim = [min(min(min(MovieCells{FrameIndex}(:,:,CurrentZ)))), 20];
    imagesc(axesTemp, SnippetImage, clim);
    colormap('gray') 
    hold on
    plot(pline_x, pline_y,'.r');
    row_pixels= Spots(FrameIndex).Fits(SpotIndex).bwPixelList(:,1).'-double(xmin)+1;
    col_pixels= Spots(FrameIndex).Fits(SpotIndex).bwPixelList(:,2).'-double(ymin)+1;
    scatter(row_pixels,col_pixels, 'r*', 'MarkerFaceAlpha', .5, 'MarkerEdgeAlpha', .5);
    [high_row, high_cols] = find(SnippetImage > 20);
    scatter(high_cols,high_row,  'y*',  'MarkerFaceAlpha', .5, 'MarkerEdgeAlpha', .5);
    axis off
    title(plot_title)
    
    if isempty(BadIndices{FrameIndex})
        BadIndices{FrameIndex} = [];
    end
    
    if isempty(GoodIndices{FrameIndex})
        GoodIndices{FrameIndex} = [];
    end
    
    if isempty(MultiIndices{FrameIndex})
        MultiIndices{FrameIndex} = [];
    end
    if ismember(SpotIndex, BadIndices{FrameIndex})
        set(TempFigure,'Color','r')
    elseif ismember(SpotIndex, GoodIndices{FrameIndex})
        set(TempFigure,'Color','g')
    elseif ismember(SpotIndex, MultiIndices{FrameIndex})
     set(TempFigure,'Color','b')
    else
        set(TempFigure,'Color','y')
    end
    
    hold off
    %%

    ct=waitforbuttonpress;
    currentCharacter=get(TempFigure,'currentcharacter');
    currentMouse=get(TempFigure,'CurrentPoint');

    if (ct~=0)&(currentCharacter=='.')&(SpotIndex<length(Spots(FrameIndex).Fits))
        SpotIndex=SpotIndex+1;
    elseif (ct~=0)&(currentCharacter==',')&(SpotIndex>1)
        SpotIndex=SpotIndex-1;
        
    elseif (ct~=0)&(currentCharacter=='>')&(FrameIndex<length(Spots))
        FrameIndex=FrameIndex+1;
        SpotIndex = 1;
  
    elseif (ct~=0)&(currentCharacter=='<')&(FrameIndex>1)
        FrameIndex=FrameIndex-1;
        SpotIndex = 1;

   
    elseif (ct~=0)&(currentCharacter=='y')
        GoodIndices{FrameIndex} =[GoodIndices{FrameIndex} SpotIndex];
        if ismember(SpotIndex, BadIndices{FrameIndex})
            BadIndices{FrameIndex} = BadIndices{FrameIndex}(BadIndices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, MultiIndices{FrameIndex})
            MultiIndices{FrameIndex} = MultiIndices{FrameIndex}(MultiIndices{FrameIndex} ~= SpotIndex);
        end
    elseif (ct~=0)&(currentCharacter=='n')
        BadIndices{FrameIndex} =[BadIndices{FrameIndex} SpotIndex];
        if ismember(SpotIndex, GoodIndices{FrameIndex})
            GoodIndices{FrameIndex} = GoodIndices{FrameIndex}(GoodIndices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, MultiIndices{FrameIndex})
            MultiIndices{FrameIndex} = MultiIndices{FrameIndex}(MultiIndices{FrameIndex} ~= SpotIndex);
        end
   elseif (ct~=0)&(currentCharacter=='j')
       try
           iJump = inputdlg('Number of Spots to jump forward:', ...
               'Move forward');
           iJump = str2double(iJump{1});
       catch
           iJump = 0;
       end
       if (SpotIndex + iJump <= length(Spots(FrameIndex).Fits))
           SpotIndex = SpotIndex+iJump;
       else
           SpotIndex = length(Spots(FrameIndex).Fits);
       end
    elseif (ct~=0)&(currentCharacter=='k')
       try
           iJump = inputdlg('Number of Spots to jump backward:', ...
               'Move backward');
           iJump = str2double(iJump{1});
       catch
           iJump = 0;
       end
       if (SpotIndex - iJump >= 1)
           SpotIndex = SpotIndex-iJump;
       else
           SpotIndex = 1;
       end     
   elseif (ct~=0)&(currentCharacter=='o')
        MultiIndices{FrameIndex} =[MultiIndices{FrameIndex} SpotIndex];
        if ismember(SpotIndex, GoodIndices{FrameIndex})
            GoodIndices{FrameIndex} = GoodIndices{FrameIndex}(GoodIndices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, BadIndices{FrameIndex})
            BadIndices{FrameIndex} = BadIndices{FrameIndex}(BadIndices{FrameIndex} ~= SpotIndex);
        end
    elseif (ct~=0)&(currentCharacter=='r')
        if ismember(SpotIndex, MultiIndices{FrameIndex})
            MultiIndices{FrameIndex} = MultiIndices{FrameIndex}(MultiIndices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, GoodIndices{FrameIndex})
            GoodIndices{FrameIndex} = GoodIndices{FrameIndex}(GoodIndices{FrameIndex} ~= SpotIndex);
        end
        if ismember(SpotIndex, BadIndices{FrameIndex})
            BadIndices{FrameIndex} = BadIndices{FrameIndex}(BadIndices{FrameIndex} ~= SpotIndex);
        end
    elseif (ct~=0)&(currentCharacter=='p')
        disp(['Frame Good Spots: ', num2str(length(GoodIndices{FrameIndex})),...
           ', Frame Bad Spots: ', num2str(length(BadIndices{FrameIndex})),...
           ', Frame Multi Spots: ', num2str(length(MultiIndices{FrameIndex}))]);
     elseif (ct~=0)&(currentCharacter=='s')
         GoodSpots = Spots;
         BadSpots= Spots;
         MultiSpots = Spots;
         for i = 1:length(Spots)
             GoodSpots(i).Fits = Spots(i).Fits(GoodIndices{i}) ;
             BadSpots(i).Fits = Spots(i).Fits(BadIndices{i}) ;
             MultiSpots(i).Fits = Spots(i).Fits(MultiIndices{i}) ;
         end


         outpath = [liveExperiment.resultsFolder, filesep, 'StoreSpotInfo.mat'];
         save(outpath, 'GoodSpots', 'BadSpots', 'MultiSpots', 'GoodIndices', 'BadIndices', 'MultiIndices');

         
    end
end
close all
GoodSpots = Spots;
BadSpots= Spots;
MultiSpots = Spots;
for i = 1:length(Spots)
    GoodSpots(i).Fits = Spots(i).Fits(GoodIndices{i}) ;
    BadSpots(i).Fits = Spots(i).Fits(BadIndices{i}) ;
    MultiSpots(i).Fits = Spots(i).Fits(MultiIndices{i}) ;
end


outpath = [liveExperiment.resultsFolder, filesep, 'StoreSpotInfo.mat'];
save(outpath, 'GoodSpots', 'BadSpots', 'MultiSpots', 'GoodIndices', 'BadIndices', 'MultiIndices');


%%