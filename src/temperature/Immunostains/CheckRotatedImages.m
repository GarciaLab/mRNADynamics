function CheckRotatedImages(Prefix)

liveExperiment = LiveExperiment(Prefix);
PixelSize_um = liveExperiment.pixelSize_um;
xDim = liveExperiment.xDim;
yDim = liveExperiment.yDim;
FrameInfo = getFrameInfo(liveExperiment);
load([liveExperiment.resultsFolder, filesep, 'MarkAndFindInfo.Mat'], 'MarkAndFindInfo');

NEmbryos = MarkAndFindInfo.NSeries;

%%
RotatedMembraneMat = imreadStack2([liveExperiment.preFolder, filesep,...
    liveExperiment.Prefix, '-Membrane_Rotated.tif'], liveExperiment.yDim,...
    liveExperiment.xDim, liveExperiment.nFrames);


CompiledEmbryoPath = [liveExperiment.resultsFolder, filesep, 'CompiledEmbryos.Mat'];
load(CompiledEmbryoPath);

%%

close all
CurrentEmbryo = 1;
EmbryoImage = RotatedMembraneMat(:,:,CurrentEmbryo);
DisplayRange=[min(min(EmbryoImage)),max(max(EmbryoImage))];


EmbryoFigure=figure;
set(EmbryoFigure,'units', 'normalized', 'position',[0.01, .1, .7, .7]);

embryoAxes = axes(EmbryoFigure,'Units', 'normalized', 'Position', [0 0 1 1]);
cc=1;

% Show the first image
imEmbryo = imshow(EmbryoImage,DisplayRange,'Parent',embryoAxes);
hold on

AnteriorPole = plot(CompiledEmbryos.RotatedCoordAs(CurrentEmbryo,1), CompiledEmbryos.RotatedCoordAs(CurrentEmbryo,2),...
    'g.','MarkerSize',20);
PosteriorPole = plot(CompiledEmbryos.RotatedCoordPs(CurrentEmbryo,1), CompiledEmbryos.RotatedCoordPs(CurrentEmbryo,2),...
    'r.','MarkerSize',20);

DorsalEdge = plot(CompiledEmbryos.RotatedCoordDs(CurrentEmbryo,1), CompiledEmbryos.RotatedCoordDs(CurrentEmbryo,2),...
    'b.','MarkerSize',20);
VentralEdge = plot(CompiledEmbryos.RotatedCoordVs(CurrentEmbryo,1), CompiledEmbryos.RotatedCoordVs(CurrentEmbryo,2),...
    'y.','MarkerSize',20);

if CompiledEmbryos.RotatedCoordPs(CurrentEmbryo,1) > 0 & CompiledEmbryos.RotatedCoordAs(CurrentEmbryo,1) > 0
    Slope = (CompiledEmbryos.RotatedCoordPs(CurrentEmbryo,2)-CompiledEmbryos.RotatedCoordAs(CurrentEmbryo,2))/(CompiledEmbryos.RotatedCoordPs(CurrentEmbryo,1)-CompiledEmbryos.RotatedCoordAs(CurrentEmbryo,1));
    Intercept = CompiledEmbryos.RotatedCoordAs(CurrentEmbryo,2)-Slope*CompiledEmbryos.RotatedCoordAs(CurrentEmbryo,1);
    MidPoint = [(CompiledEmbryos.RotatedCoordAs(CurrentEmbryo,1)+CompiledEmbryos.RotatedCoordPs(CurrentEmbryo,1))/2, (CompiledEmbryos.RotatedCoordAs(CurrentEmbryo,2)+CompiledEmbryos.RotatedCoordPs(CurrentEmbryo,2))/2];
    DVSlope = -1/Slope;
    DVIntercept = MidPoint(2)-DVSlope*MidPoint(1);
else
    Slope = 0;
    Intercept= 0;
    DVSlope=0;
    DVIntercept= 0;
    MidPoint = [xDim/2, xDim/2];
end


APaxis = plot([0 xDim], Slope*[0 xDim]+Intercept,'c-');
if abs(Slope)  > 10^(-6);
    DVaxis = plot([0 xDim], DVSlope*[0 xDim]+DVIntercept,'c-');
else
    DVaxis = plot([MidPoint(1) MidPoint(1)], [0 yDim],'c-');
end
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

while (cc~='x')
    EmbryoImage = RotatedMembraneMat(:,:,CurrentEmbryo);
    
    imEmbryo.CData = EmbryoImage;
    try
        caxis(embryoAxes, DisplayRange);
    end
    AnteriorPole.XData = CompiledEmbryos.RotatedCoordAs(CurrentEmbryo,1);
    AnteriorPole.YData = CompiledEmbryos.RotatedCoordAs(CurrentEmbryo,2);
    PosteriorPole.XData = CompiledEmbryos.RotatedCoordPs(CurrentEmbryo,1);
    PosteriorPole.YData = CompiledEmbryos.RotatedCoordPs(CurrentEmbryo,2);
    DorsalEdge.XData = CompiledEmbryos.RotatedCoordDs(CurrentEmbryo,1);
    DorsalEdge.YData = CompiledEmbryos.RotatedCoordDs(CurrentEmbryo,2);
    VentralEdge.XData = CompiledEmbryos.RotatedCoordVs(CurrentEmbryo,1);
    VentralEdge.YData = CompiledEmbryos.RotatedCoordVs(CurrentEmbryo,2);
    
    %     imshow(imadjust(APImage), 'Parent', apAx)
    %imshow(APImage,DisplayRange)
    if CompiledEmbryos.RotatedCoordPs(CurrentEmbryo,1) > 0 & CompiledEmbryos.RotatedCoordAs(CurrentEmbryo,1) > 0
        Slope = (CompiledEmbryos.RotatedCoordPs(CurrentEmbryo,2)-CompiledEmbryos.RotatedCoordAs(CurrentEmbryo,2))/(CompiledEmbryos.RotatedCoordPs(CurrentEmbryo,1)-CompiledEmbryos.RotatedCoordAs(CurrentEmbryo,1));
        Intercept = CompiledEmbryos.RotatedCoordAs(CurrentEmbryo,2)-Slope*CompiledEmbryos.RotatedCoordAs(CurrentEmbryo,1);
        MidPoint = [(CompiledEmbryos.RotatedCoordAs(CurrentEmbryo,1)+CompiledEmbryos.RotatedCoordPs(CurrentEmbryo,1))/2, (CompiledEmbryos.RotatedCoordAs(CurrentEmbryo,2)+CompiledEmbryos.RotatedCoordPs(CurrentEmbryo,2))/2];
        DVSlope = -1/Slope;
        DVIntercept = MidPoint(2)-DVSlope*MidPoint(1);
    else
        Slope = 0;
        Intercept= 0;
        DVSlope=0;
        DVIntercept= 0;
        MidPoint = [xDim/2, xDim/2];
    end
    APaxis.YData = Slope*[0 xDim]+Intercept;
    if abs(Slope)  > 10^(-6);
        DVaxis.YData =DVSlope*[0 xDim]+DVIntercept;
        DVaxis.XData =[0 xDim];
    else
        DVaxis.YData =[0 yDim]
        DVaxis.XData =[MidPoint(1) MidPoint(1)];
        
    end
    
    
    
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
        CompiledEmbryos.RotatedCoordAs(CurrentEmbryo,:)=0;
        CompiledEmbryos.RotatedCoordPs(CurrentEmbryo,:)=0;
        CompiledEmbryos.RotatedCoordDs(CurrentEmbryo,:)=0;
        CompiledEmbryos.RotatedCoordVs(CurrentEmbryo,:)=0;
        CompiledEmbryos.CoordAs(CurrentEmbryo,:)=0;
        CompiledEmbryos.CoordPs(CurrentEmbryo,:)=0;
        CompiledEmbryos.CoordDs(CurrentEmbryo,:)=0;
        CompiledEmbryos.CoordVs(CurrentEmbryo,:)=0;
    elseif (ct~=0)&(cc=='a')	%Select anterior end
        [CoordAx,CoordAy]=ginputc(1,'Color',[1,1,1]);
        CompiledEmbryos.RotatedCoordAs(CurrentEmbryo,:) = [CoordAx,CoordAy];
    elseif (ct~=0)&(cc=='p')    %Select posterior end
        [CoordPx,CoordPy]=ginputc(1,'Color',[1,1,1]);
        CompiledEmbryos.RotatedCoordPs(CurrentEmbryo,:) = [CoordPx,CoordPy];
    elseif (ct~=0)&(cc=='d')	%Select dorsal edge
        [CoordDx,CoordDy]=ginputc(1,'Color',[1,1,1]);
        CompiledEmbryos.RotatedCoordDs(CurrentEmbryo,:) = [CoordDx,CoordDy];
    elseif (ct~=0)&(cc=='v')    %Select ventral edge
        [CoordVx,CoordVy]=ginputc(1,'Color',[1,1,1]);
        CompiledEmbryos.RotatedCoordVs(CurrentEmbryo,:) = [CoordVx,CoordVy];
    elseif (ct~=0)&(cc=='.')&(CurrentEmbryo < NEmbryos)    %Increase contrast
        CurrentEmbryo = CurrentEmbryo+1;
        EmbryoImage = RotatedMembraneMat(:,:,CurrentEmbryo);
        DisplayRange=[min(min(EmbryoImage)),max(max(EmbryoImage))];
    elseif (ct~=0)&(cc==',')&(CurrentEmbryo > 1)    %Increase contrast
        CurrentEmbryo = CurrentEmbryo-1;
        EmbryoImage = RotatedMembraneMat(:,:,CurrentEmbryo);
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
        CompiledEmbryos.CoordAs = TransformToOriginalImageCoordinates(CompiledEmbryos.RotatedCoordAs, Prefix,CompiledEmbryos);
        CompiledEmbryos.CoordPs = TransformToOriginalImageCoordinates(CompiledEmbryos.RotatedCoordPs, Prefix,CompiledEmbryos);
        CompiledEmbryos.CoordDs = TransformToOriginalImageCoordinates(CompiledEmbryos.RotatedCoordDs, Prefix,CompiledEmbryos);
        CompiledEmbryos.CoordVs = TransformToOriginalImageCoordinates(CompiledEmbryos.RotatedCoordVs, Prefix,CompiledEmbryos);
        CompiledEmbryos = UpdateAPAxisInfo(Prefix, CompiledEmbryos);
        
        
        save(CompiledEmbryoPath,'CompiledEmbryos');
        disp('Axis information saved.')
    end
end

CompiledEmbryos.CoordAs = TransformToOriginalImageCoordinates(CompiledEmbryos.RotatedCoordAs, Prefix,CompiledEmbryos);
CompiledEmbryos.CoordPs = TransformToOriginalImageCoordinates(CompiledEmbryos.RotatedCoordPs, Prefix,CompiledEmbryos);
CompiledEmbryos.CoordDs = TransformToOriginalImageCoordinates(CompiledEmbryos.RotatedCoordDs, Prefix,CompiledEmbryos);
CompiledEmbryos.CoordVs = TransformToOriginalImageCoordinates(CompiledEmbryos.RotatedCoordVs, Prefix,CompiledEmbryos);
CompiledEmbryos = UpdateAPAxisInfo(Prefix, CompiledEmbryos);


save(CompiledEmbryoPath,'CompiledEmbryos');
disp('Axis information saved.')
close all
RotateEmbryoImages(Prefix);



