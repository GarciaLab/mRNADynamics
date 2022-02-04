function chooseCustomMaxProjections(Prefix)
liveExperiment = LiveExperiment(Prefix);

DropboxFolder = liveExperiment.userResultsFolder;
PreProcPath = liveExperiment.userPreFolder;

ProjectionMatPath = [liveExperiment.resultsFolder, filesep, 'ProjectionMat.Mat'];
movieMat = getFirstRepMat(liveExperiment);

xSize = size(movieMat,2);
ySize = size(movieMat,1);
zSize = size(movieMat,3);
NEmbryos = size(movieMat,4);

SubplotDims = numSubplots(zSize+1);


Channel1 = liveExperiment.Channel1;
Channel2 = liveExperiment.Channel2;
Channel3 = liveExperiment.Channel3;
Channels = {Channel1, Channel2, Channel3};
HisChannel= 1;
for i = 1:length(Channels)
    if contains(lower(Channels{i}), 'his')
        HisChannel = i;
    end
end



hisStacks = squeeze(movieMat(:,:,:,:,HisChannel));

if max(max(max(max(hisStacks)))) < 256
    hisPrecision = 'uint8';
else
    
    hisPrecision = 'uint16';
end

hisMovie = zeros(ySize, xSize, NEmbryos, 'double');
if isfile(ProjectionMatPath)
    load(ProjectionMatPath, 'ProjectionMat');
else
    ProjectionMat = ones( NEmbryos, zSize, 'logical');
end

for i = 1:NEmbryos
    EmbryoImageStack =hisStacks(:,:,:,i);
    hisMovie(:,:,i) = max(EmbryoImageStack(:,:,ProjectionMat(i,:)), [], 3);
    hisMovie(:, :, i)  = mat2gray(hisMovie(:, :, i));
    hisMovie(:, :, i) = hisMovie(:, :, i) *255;
end
MaxFrames = cell(1,NEmbryos);

CurrentEmbryo = 1;
EmbryoImageStack =hisStacks(:,:,:,CurrentEmbryo);
DisplayRange=[min(min(min(EmbryoImageStack))),max(max(max(EmbryoImageStack)))];
MaxImage = max(EmbryoImageStack(:,:,ProjectionMat(CurrentEmbryo,:)), [], 3);
MaxDisplayRange = [min(min(MaxImage)), max(max(MaxImage))];


%%
close all
EmbryoFigure=figure;
set(EmbryoFigure,'units', 'normalized', 'position',[0.01, .1, .7, .7]);

embryoAxes = axes(EmbryoFigure,'Units', 'normalized', 'Position', [0 0 1 1]);





% Show the first image
imEmbryo = imshow(MaxImage,MaxDisplayRange,'Parent',embryoAxes);




% axis image
axis off


IncludedZsString = 'z = ';
for i = 1:zSize
    if ProjectionMat(CurrentEmbryo,i)
        IncludedZsString= [IncludedZsString, num2str(i),', '];
    end
end

IncludedZsString = IncludedZsString(1:end-2);

FigureTitle=['Embryo: ',num2str(CurrentEmbryo),'/',num2str(NEmbryos),...
    ', ',IncludedZsString];



set(EmbryoFigure,'Name',FigureTitle)
zSlices = 1:zSize;
cc=1;
%%
while (cc~='x')
    
    
    imEmbryo.CData = MaxImage;
    try
        caxis(embryoAxes, MaxDisplayRange);
    end
    
    
    
    
    hold off
    
    
    IncludedZsString = 'z = ';
    for i = 1:zSize
        if ProjectionMat(CurrentEmbryo,i)
            IncludedZsString= [IncludedZsString, num2str(i),', '];
        end
    end
    
    IncludedZsString = IncludedZsString(1:end-2);
    
    FigureTitle=['Embryo: ',num2str(CurrentEmbryo),'/',num2str(NEmbryos),...
        ', ',IncludedZsString];
    
    set(EmbryoFigure,'Name',FigureTitle)
    
    figure(EmbryoFigure)
    ct=waitforbuttonpress;
    cc=get(EmbryoFigure,'currentcharacter');
    cm=get(gca,'CurrentPoint');
    
    
    if (ct~=0)&((cc=='1')|(cc=='2')|(cc=='3')|(cc=='4')|(cc=='5')|(cc=='6')|...
            (cc=='7')|(cc=='8')|(cc=='9'))
        cnum = str2num(cc);
        if cnum > zSize
            disp([cc, ' is larger than the number of slices (', num2str(zSize),')']);
        elseif all(~ProjectionMat(CurrentEmbryo,zSlices ~= cnum))
            disp(['Cannot unselect this slice because no other slices are currently selected.']);
        else
            ProjectionMat(CurrentEmbryo, cnum) = ~ProjectionMat(CurrentEmbryo, cnum);
            MaxImage = max(EmbryoImageStack(:,:,ProjectionMat(CurrentEmbryo,:)), [], 3);
            MaxDisplayRange = [min(min(MaxImage)), max(max(MaxImage))];
            hisMovie(:,:,CurrentEmbryo) = MaxImage;
            MaxFrames{CurrentEmbryo} = MaxImage;
        end
    elseif (ct~=0)&(cc=='c')
        ProjectionMat(CurrentEmbryo, :) = true;
        MaxImage = max(EmbryoImageStack(:,:,ProjectionMat(CurrentEmbryo,:)), [], 3);
        MaxDisplayRange = [min(min(MaxImage)), max(max(MaxImage))];
        hisMovie(:,:,CurrentEmbryo) = MaxImage;
        MaxFrames{CurrentEmbryo} = MaxImage;
    elseif (ct~=0)&(cc=='.')&(CurrentEmbryo < NEmbryos)    %Increase contrast
        CurrentEmbryo = CurrentEmbryo+1;
        EmbryoImageStack =hisStacks(:,:,:,CurrentEmbryo);
        MaxImage = max(EmbryoImageStack(:,:,ProjectionMat(CurrentEmbryo,:)), [], 3);
        MaxDisplayRange = [min(min(MaxImage)), max(max(MaxImage))];
        hisMovie(:,:,CurrentEmbryo) = MaxImage;
        MaxFrames{CurrentEmbryo} = MaxImage;
    elseif (ct~=0)&(cc==',')&(CurrentEmbryo > 1)    %Increase contrast
        CurrentEmbryo = CurrentEmbryo-1;
        EmbryoImageStack =hisStacks(:,:,:,CurrentEmbryo);
        MaxImage = max(EmbryoImageStack(:,:,ProjectionMat(CurrentEmbryo,:)), [], 3);
        MaxDisplayRange = [min(min(MaxImage)), max(max(MaxImage))];
        hisMovie(:,:,CurrentEmbryo) = MaxImage;
        MaxFrames{CurrentEmbryo} = MaxImage;
        
    elseif (ct~=0)&(cc=='m')    %Increase contrast
        MaxDisplayRange(2)=MaxDisplayRange(2)/2;
    elseif (ct~=0)&(cc=='n')    %Increase contrast
        MaxDisplayRange(2)=MaxDisplayRange(2)*2;
    elseif (ct~=0)&(cc=='j')
        try
            iJump = inputdlg('Embryo to jump to:', ...
                'Move to frame');
            iJump = str2double(iJump{1});
        catch
            iJump =CurrentEmbryo;
        end
        if (iJump >= 1) & (iJump <= NEmbryos)
            CurrentEmbryo= iJump;
        end
        EmbryoImageStack =hisStacks(:,:,:,CurrentEmbryo);
        MaxImage = max(EmbryoImageStack(:,:,ProjectionMat(CurrentEmbryo,:)), [], 3);
        MaxDisplayRange = [min(min(MaxImage)), max(max(MaxImage))];
        hisMovie(:,:,CurrentEmbryo) = MaxImage;
        MaxFrames{CurrentEmbryo} = MaxImage;
    elseif (ct~=0)&(cc=='r')    %Reset the contrast
        MaxDisplayRange = [min(min(MaxImage)), max(max(MaxImage))];
    elseif (ct~=0)&(cc=='s')
        save(ProjectionMatPath,'ProjectionMat');
    end
end
close all
for i = 1:NEmbryos
    EmbryoImageStack =hisStacks(:,:,:,i);
    hisMovie(:,:,i) = max(EmbryoImageStack(:,:,ProjectionMat(i,:)), [], 3);
    hisMovie(:, :, i)  = mat2gray(hisMovie(:, :, i));
    hisMovie(:, :, i) = hisMovie(:, :, i) *255;
end
save(ProjectionMatPath,'ProjectionMat');
saveNuclearProjection(hisMovie, [liveExperiment.preFolder, filesep, Prefix,'-CustomHis.tif']);
%%