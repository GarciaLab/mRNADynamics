function UseSliceInfo = CheckStandardCandleSlices(Prefix)
liveExperiment = LiveExperiment(Prefix);
SeriesMatches = CheckSettings(Prefix);

if isfile([liveExperiment.resultsFolder, filesep, 'Spots.mat']);
    Spots = getSpots(liveExperiment);   
end

NumFrames = liveExperiment.nFrames;

TifList = dir([liveExperiment.preFolder, filesep, '*ch01.tif']);
ProbList = dir([liveExperiment.procFolder,'dogs',  filesep, 'prob*']);
MovieCells = cell(1, NumFrames);
SpotMat = cell(1, NumFrames);
NumSlices = 0;
for i = 1:NumFrames
    MovieCells{i} = imreadStack([TifList(i).folder, filesep,  TifList(i).name]);
    NumSlices = NumSlices + size(MovieCells{i}, 3);
end
%%
StackIDs = zeros(1, NumSlices, 'uint16');
zIDs = zeros(1, NumSlices, 'uint16');
MovieMat = zeros(liveExperiment.xDim,  liveExperiment.yDim, NumSlices);
SpotMat = zeros(liveExperiment.xDim,  liveExperiment.yDim, NumSlices, 'uint16');
z_counter = 0;
UseSeriesSliceInfo = ones(1, NumSlices, 'logical');
for i = 1:NumFrames
    MovieMat(:,:,(z_counter+1):(z_counter+size(MovieCells{i},3))) = ...
      MovieCells{i};
    SpotMat(:,:,(z_counter+1):(z_counter+size(MovieCells{i},3))) = ...
      imreadStack([ProbList(i).folder, filesep,  ProbList(i).name]);
    StackIDs((z_counter+1):(z_counter+size(MovieCells{i},3))) = i;
    zIDs((z_counter+1):(z_counter+size(MovieCells{i},3))) = 1:size(MovieCells{i},3);
    UseSeriesSliceInfo((z_counter+1):(z_counter+size(MovieCells{i},3))) = SeriesMatches(i);
    z_counter = z_counter + size(MovieCells{i},3);
  
end
    

%%
if isfile([liveExperiment.resultsFolder,filesep,'UseSliceInfo.mat'])
    load([liveExperiment.resultsFolder,filesep,'UseSliceInfo.mat'],'UseSliceInfo');
else
    UseSliceInfo = true(1, NumSlices);
end

UseSliceInfo = UseSliceInfo & UseSeriesSliceInfo;

SliceIndex = 1;

close all
TempFigure = figure;
set(TempFigure, 'units', 'normalized', 'position', [0.2, 0.2, .6, .5]);
axesTemp = axes(TempFigure);
axis off
plot_title= {strrep(TifList(StackIDs(SliceIndex)).name, '_', '\_'), ['z: ', num2str(zIDs(SliceIndex))]};
FigAx{1} = subplot(1,2,1);
imagesc(MovieMat(:,:,SliceIndex));
axis off
FigAx{2} = subplot(1,2,2);
imagesc(SpotMat(:,:,SliceIndex));
hold on 
if exist('Spots', 'var')
xVector = [Spots(StackIDs(SliceIndex)).Fits(:).xDoG];
yVector = [Spots(StackIDs(SliceIndex)).Fits(:).yDoG];
zSpots = [Spots(StackIDs(SliceIndex)).Fits(:).z] == zIDs(SliceIndex);
if sum(zSpots) == 0
  UseSliceInfo(SliceIndex) = false;
end
scatter(xVector(zSpots), yVector(zSpots), 'r.')
end
axis off 
sgtitle(plot_title)
if UseSliceInfo(SliceIndex)
     set(TempFigure,'Color','g')
else
     set(TempFigure,'Color','r')
end


currentCharacter = 0;
set(0, 'CurrentFigure', TempFigure)
while (currentCharacter~='x')
    
    %Load subsequent images
    if SliceIndex ~= NumSlices
        plot_title= {strrep(TifList(StackIDs(SliceIndex)).name, '_', '\_'), ['z: ', num2str(zIDs(SliceIndex))]};
    else
        plot_title= {strrep(TifList(StackIDs(SliceIndex)).name, '_', '\_'), ['z: ', num2str(zIDs(SliceIndex)), ' (Last Frame)']};
    end
    imagesc(FigAx{1} , MovieMat(:,:,SliceIndex));
    set(FigAx{1}, 'visible', 'off')
    imagesc(FigAx{2}, SpotMat(:,:,SliceIndex));
    axis off
    hold on
    if exist('Spots', 'var')
    xVector = [Spots(StackIDs(SliceIndex)).Fits(:).xDoG];
    yVector = [Spots(StackIDs(SliceIndex)).Fits(:).yDoG];
    zSpots = [Spots(StackIDs(SliceIndex)).Fits(:).z] == zIDs(SliceIndex);
    if sum(zSpots) == 0
        UseSliceInfo(SliceIndex) = false;
    end
    scatter(FigAx{2}, xVector(zSpots), yVector(zSpots), 'r.')
    end
    set(FigAx{2}, 'visible', 'off')
    sgtitle(plot_title)
    if UseSliceInfo(SliceIndex)
        set(TempFigure,'Color','g')
    else
        set(TempFigure,'Color','r')
    end
    %%

    ct=waitforbuttonpress;
    currentCharacter=get(TempFigure,'currentcharacter');
    currentMouse=get(TempFigure,'CurrentPoint');

    if (ct~=0)&(currentCharacter=='.')&(SliceIndex<NumSlices)
        SliceIndex=SliceIndex+1;
    elseif (ct~=0)&(currentCharacter==',')&(SliceIndex>1)
        SliceIndex=SliceIndex-1;
    elseif (ct~=0)&(currentCharacter=='s')
        save([liveExperiment.resultsFolder,filesep,'UseSliceInfo.mat'],'UseSliceInfo', '-v6');
        disp('UseSliceInfo saved.');
    elseif (ct~=0)&(currentCharacter=='y')
        UseSliceInfo(SliceIndex) = true;
    elseif (ct~=0)&(currentCharacter=='n')
        UseSliceInfo(SliceIndex) = false;
   
        
    end
end

close all
save([liveExperiment.resultsFolder,filesep,'UseSliceInfo.mat'],'UseSliceInfo', 'StackIDs','zIDs', '-v6');
disp('UseSliceInfo saved.');
%%