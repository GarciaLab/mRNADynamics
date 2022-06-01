function [CompiledEmbryos] = DefineMembraneAxesAllEmbryos(Prefix)

%%
liveExperiment = LiveExperiment(Prefix);
FrameInfo = getFrameInfo(liveExperiment);
load([liveExperiment.resultsFolder, filesep, 'MarkAndFindInfo.Mat'], 'MarkAndFindInfo');
    PreProcFolder = liveExperiment.preFolder;
    
NEmbryos = MarkAndFindInfo.NSeries;
MembraneZoomResultsFolder = [liveExperiment.resultsFolder, filesep, 'ZoomMembraneInfo'];


MembraneZoomPixelPath = [MembraneZoomResultsFolder, filesep, 'MembraneZoomPixelSize.mat'];
load(MembraneZoomPixelPath,'PixelSize_um');


%%
MedNameSuffix = 'MedZoomMembrane';
    MedNewName = [Prefix,'-',MedNameSuffix, '.tif'];
    
    
    medMembraneFile = [PreProcFolder, filesep, MedNewName];

if exist(medMembraneFile, 'file')
    haveMemTifStack = true;
else
    haveMemTifStack = false;
end





if haveMemTifStack
    %load in sequential tif stacks
    membraneMat = imreadStack2(medMembraneFile);
    
    
else
    %load in individual tif slices
    error('Membrane Tif Stack does not exist. Run ExportDataForLivemRNA.')
    
end



%             end
membraneMat = double(membraneMat);
xDim = size(membraneMat, 2);
yDim = size(membraneMat, 1);
NEmbryos = size(membraneMat, 3);
%%
CompiledEmbryoPath = [liveExperiment.resultsFolder, filesep, 'CompiledEmbryos.Mat'];
if isfile(CompiledEmbryoPath)
    load(CompiledEmbryoPath,'CompiledEmbryos');
    
    if ~isfield(CompiledEmbryos, 'MemCoordAs')
        CompiledEmbryos.MemCoordAs = CompiledEmbryos.CoordAs/(PixelSize_um/liveExperiment.pixelSize_um);
        CompiledEmbryos.MemCoordPs =CompiledEmbryos.CoordPs/(PixelSize_um/liveExperiment.pixelSize_um);
        CompiledEmbryos.MemCoordDs = CompiledEmbryos.CoordDs/(PixelSize_um/liveExperiment.pixelSize_um);
        CompiledEmbryos.MemCoordVs = CompiledEmbryos.CoordVs/(PixelSize_um/liveExperiment.pixelSize_um);
    end
else
    error('Need to first run DefineAPAxesForAllEmbryos');
    
end






%%
close all
CurrentEmbryo = 1;
EmbryoImage = membraneMat(:,:,CurrentEmbryo);
DisplayRange=[min(min(EmbryoImage)),max(max(EmbryoImage))];


EmbryoFigure=figure;
set(EmbryoFigure,'units', 'normalized', 'position',[0.01, .1, .7, .7]);

embryoAxes = axes(EmbryoFigure,'Units', 'normalized', 'Position', [0 0 1 1]);
cc=1;

% Show the first image
imEmbryo = imshow(EmbryoImage,DisplayRange,'Parent',embryoAxes);
hold on

AnteriorPole = plot(CompiledEmbryos.MemCoordAs(CurrentEmbryo,1), CompiledEmbryos.MemCoordAs(CurrentEmbryo,2),...
    'g.','MarkerSize',20);
PosteriorPole = plot(CompiledEmbryos.MemCoordPs(CurrentEmbryo,1), CompiledEmbryos.MemCoordPs(CurrentEmbryo,2),...
    'r.','MarkerSize',20);

DorsalEdge = plot(CompiledEmbryos.MemCoordDs(CurrentEmbryo,1), CompiledEmbryos.MemCoordDs(CurrentEmbryo,2),...
    'b.','MarkerSize',20);
VentralEdge = plot(CompiledEmbryos.MemCoordVs(CurrentEmbryo,1), CompiledEmbryos.MemCoordVs(CurrentEmbryo,2),...
    'y.','MarkerSize',20);

if CompiledEmbryos.MemCoordPs(CurrentEmbryo,1) > 0 & CompiledEmbryos.MemCoordAs(CurrentEmbryo,1) > 0
    Slope = (CompiledEmbryos.MemCoordPs(CurrentEmbryo,2)-CompiledEmbryos.MemCoordAs(CurrentEmbryo,2))/(CompiledEmbryos.MemCoordPs(CurrentEmbryo,1)-CompiledEmbryos.MemCoordAs(CurrentEmbryo,1));
    Intercept = CompiledEmbryos.MemCoordAs(CurrentEmbryo,2)-Slope*CompiledEmbryos.MemCoordAs(CurrentEmbryo,1);
    MidPoint = [(CompiledEmbryos.MemCoordAs(CurrentEmbryo,1)+CompiledEmbryos.MemCoordPs(CurrentEmbryo,1))/2, (CompiledEmbryos.MemCoordAs(CurrentEmbryo,2)+CompiledEmbryos.MemCoordPs(CurrentEmbryo,2))/2];
    DVSlope = -1/Slope;
    DVIntercept = MidPoint(2)-DVSlope*MidPoint(1);
else
    Slope = 0;
    Intercept= 0;
    DVSlope=0;
    DVIntercept= 0;
end


APaxis = plot([0 xDim], Slope*[0 xDim]+Intercept,'c-');
DVaxis = plot([0 xDim], DVSlope*[0 xDim]+DVIntercept,'c-');
if (CompiledEmbryos.Approved(CurrentEmbryo))
    set(EmbryoFigure,'Color','g')
else
    set(EmbryoFigure,'Color','r')
end
% axis image
axis off
title('Anterior (green), posterior (red); original')
hold off

FigureTitle=['Embryo: ',num2str(CurrentEmbryo),'/',num2str(NEmbryos),...
    ', Flag: ', num2str(CompiledEmbryos.Flags(CurrentEmbryo)),', Anterior (green), Posterior (red)'];


set(EmbryoFigure,'Name',FigureTitle)
%%
while (cc~='x')
    EmbryoImage = membraneMat(:,:,CurrentEmbryo);
    
    imEmbryo.CData = EmbryoImage;
    try
        caxis(embryoAxes, DisplayRange);
    end
    AnteriorPole.XData = CompiledEmbryos.MemCoordAs(CurrentEmbryo,1);
    AnteriorPole.YData = CompiledEmbryos.MemCoordAs(CurrentEmbryo,2);
    PosteriorPole.XData = CompiledEmbryos.MemCoordPs(CurrentEmbryo,1);
    PosteriorPole.YData = CompiledEmbryos.MemCoordPs(CurrentEmbryo,2);
    DorsalEdge.XData = CompiledEmbryos.MemCoordDs(CurrentEmbryo,1);
    DorsalEdge.YData = CompiledEmbryos.MemCoordDs(CurrentEmbryo,2);
    VentralEdge.XData = CompiledEmbryos.MemCoordVs(CurrentEmbryo,1);
    VentralEdge.YData = CompiledEmbryos.MemCoordVs(CurrentEmbryo,2);
    
    %     imshow(imadjust(APImage), 'Parent', apAx)
    %imshow(APImage,DisplayRange)
    if CompiledEmbryos.MemCoordPs(CurrentEmbryo,1) > 0 & CompiledEmbryos.MemCoordAs(CurrentEmbryo,1) > 0
        Slope = (CompiledEmbryos.MemCoordPs(CurrentEmbryo,2)-CompiledEmbryos.MemCoordAs(CurrentEmbryo,2))/(CompiledEmbryos.MemCoordPs(CurrentEmbryo,1)-CompiledEmbryos.MemCoordAs(CurrentEmbryo,1));
        Intercept = CompiledEmbryos.MemCoordAs(CurrentEmbryo,2)-Slope*CompiledEmbryos.MemCoordAs(CurrentEmbryo,1);
        MidPoint = [(CompiledEmbryos.MemCoordAs(CurrentEmbryo,1)+CompiledEmbryos.MemCoordPs(CurrentEmbryo,1))/2, (CompiledEmbryos.MemCoordAs(CurrentEmbryo,2)+CompiledEmbryos.MemCoordPs(CurrentEmbryo,2))/2];
        DVSlope = -1/Slope;
        DVIntercept = MidPoint(2)-DVSlope*MidPoint(1);
    else
        Slope = 0;
        Intercept= 0;
        DVSlope=0;
        DVIntercept= 0;
    end
    APaxis.YData = Slope*[0 xDim]+Intercept;
    DVaxis.YData =DVSlope*[0 xDim]+DVIntercept;
    
    
    if (CompiledEmbryos.Approved(CurrentEmbryo))
        set(EmbryoFigure,'Color','g')
    else
        set(EmbryoFigure,'Color','r')
    end
    
    hold off
    
    
    FigureTitle=['Embryo: ',num2str(CurrentEmbryo),'/',num2str(NEmbryos),...
        ', Flag: ', num2str(CompiledEmbryos.Flags(CurrentEmbryo)),', Anterior (green), Posterior (red)'];
    
    
    set(EmbryoFigure,'Name',FigureTitle)
    
    figure(EmbryoFigure)
    ct=waitforbuttonpress;
    cc=get(EmbryoFigure,'currentcharacter');
    cm=get(gca,'CurrentPoint');
    
    
    if (ct~=0)&(cc=='c')        %Clear all AP information
        CompiledEmbryos.MemCoordAs(CurrentEmbryo,:)=0;
        CompiledEmbryos.MemCoordPs(CurrentEmbryo,:)=0;
        CompiledEmbryos.MemCoordDs(CurrentEmbryo,:)=0;
        CompiledEmbryos.MemCoordVs(CurrentEmbryo,:)=0;
    elseif (ct~=0)&(cc=='a')	%Select anterior end
        [CoordAx,CoordAy]=ginputc(1,'Color',[1,1,1]);
        CompiledEmbryos.MemCoordAs(CurrentEmbryo,:) = [CoordAx,CoordAy];
    elseif (ct~=0)&(cc=='p')    %Select posterior end
        [CoordPx,CoordPy]=ginputc(1,'Color',[1,1,1]);
        CompiledEmbryos.MemCoordPs(CurrentEmbryo,:) = [CoordPx,CoordPy];
    elseif (ct~=0)&(cc=='d')	%Select dorsal edge
        [CoordDx,CoordDy]=ginputc(1,'Color',[1,1,1]);
        CompiledEmbryos.MemCoordDs(CurrentEmbryo,:) = [CoordDx,CoordDy];
    elseif (ct~=0)&(cc=='v')    %Select ventral edge
        [CoordVx,CoordVy]=ginputc(1,'Color',[1,1,1]);
        CompiledEmbryos.MemCoordVs(CurrentEmbryo,:) = [CoordVx,CoordVy];
    elseif (ct~=0)&(cc=='.')&(CurrentEmbryo < NEmbryos)    %Increase contrast
        CurrentEmbryo = CurrentEmbryo+1;
         EmbryoImage = membraneMat(:,:,CurrentEmbryo);
        DisplayRange=[min(min(EmbryoImage)),max(max(EmbryoImage))];
    elseif (ct~=0)&(cc==',')&(CurrentEmbryo > 1)    %Increase contrast
        CurrentEmbryo = CurrentEmbryo-1;
         EmbryoImage = membraneMat(:,:,CurrentEmbryo);
        DisplayRange=[min(min(EmbryoImage)),max(max(EmbryoImage))];
    elseif (ct~=0)&(cc=='>')&(CurrentEmbryo < NEmbryos)  
        NewIndex = find(1:NEmbryos > CurrentEmbryo & ~CompiledEmbryos.Checked, 1);
        if isempty(NewIndex)
            CurrentEmbryo = NEmbryos;
        else
            CurrentEmbryo = NewIndex;
        end
        EmbryoImage = membraneMat(:,:,CurrentEmbryo);
        DisplayRange=[min(min(EmbryoImage)),max(max(EmbryoImage))];
    elseif (ct~=0)&(cc=='<')&(CurrentEmbryo > 1)  
        NewIndex = find(1:NEmbryos < CurrentEmbryo & ~CompiledEmbryos.Checked, 1,'last');
        if isempty(NewIndex)
            CurrentEmbryo = 1;
        else
            CurrentEmbryo = NewIndex;
        end
        EmbryoImage = membraneMat(:,:,CurrentEmbryo);
        DisplayRange=[min(min(EmbryoImage)),max(max(EmbryoImage))];
    elseif (ct~=0)&(cc=='m')    %Increase contrast
        DisplayRange(2)=DisplayRange(2)/2;
    elseif (ct~=0)&(cc=='n')    %Increase contrast
        DisplayRange(2)=DisplayRange(2)*2;
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
    elseif (ct~=0)&(cc=='r')    %Reset the contrast
        DisplayRange=[min(min(EmbryoImage)),max(max(EmbryoImage))];
    elseif (ct~=0)&(cc == 'f')
        try
            [flag, flag_string]  = chooseFixedEmbryoFlag;
            disp(flag_string)
            CompiledEmbryos.Flags(CurrentEmbryo) =flag;
        catch
            disp('No Flag Selected')
        end
    elseif (ct~=0)&(cc == 'q')
        CompiledEmbryos.Approved(CurrentEmbryo) = 1;
    elseif (ct~=0)&(cc == 'w')
        CompiledEmbryos.Approved(CurrentEmbryo) = 0;
        
    elseif (ct~=0)&(cc=='s')

        save(CompiledEmbryoPath,'CompiledEmbryos');
        disp('Axis information saved.')
        
        
    end
end
CompiledEmbryos = UpdateAPAxisInfo(Prefix, CompiledEmbryos);
save(CompiledEmbryoPath,'CompiledEmbryos');
disp('Axis information saved.')
close all
%%
RotateMemEmbryoImages(Prefix);