function [CompiledEmbryos] = DefineAPAxesForAllEmbryos(Prefix, MarkPV)
%%
if ~exist('MarkPV', 'var')
    MarkPV = false;
end

%%
liveExperiment = LiveExperiment(Prefix);
PixelSize_um = liveExperiment.pixelSize_um;
xDim = liveExperiment.xDim;
yDim = liveExperiment.yDim;
try
    FrameInfo = getFrameInfo(liveExperiment);
    load([liveExperiment.resultsFolder, filesep, 'MarkAndFindInfo.Mat'], 'MarkAndFindInfo');
catch 
    FrameInfo = {};
    MarkAndFindInfo = {};
end




%%
MembraneMat = getMembraneMat(liveExperiment);
if ~isempty(MarkAndFindInfo)
    NEmbryos = MarkAndFindInfo.NSeries;
else
    NEmbryos = size(MembraneMat, 3);
end

xDim = size(MembraneMat, 2);
yDim = size(MembraneMat, 1);

CompiledEmbryoPath = [liveExperiment.resultsFolder, filesep, 'CompiledEmbryos.Mat'];
if isfile(CompiledEmbryoPath)
    load(CompiledEmbryoPath,'CompiledEmbryos');

    
else
    CompiledEmbryos = {};
    CompiledEmbryos.CoordAs = zeros(NEmbryos,2);
    CompiledEmbryos.CoordPs =zeros(NEmbryos,2);
    CompiledEmbryos.CoordDs = zeros(NEmbryos,2);
    CompiledEmbryos.CoordVs = zeros(NEmbryos,2);
    CompiledEmbryos.Flags = zeros(1,NEmbryos);
    CompiledEmbryos.Approved = ones(1, NEmbryos,'logical');
end

if MarkPV & ~isfield(CompiledEmbryos, 'PVCoordAs')
    CompiledEmbryos.PVCoordAs = zeros(NEmbryos,2);
    CompiledEmbryos.PVCoordPs = zeros(NEmbryos,2);
    CompiledEmbryos.PVCoordDs = zeros(NEmbryos,2);
    CompiledEmbryos.PVCoordVs = zeros(NEmbryos,2);
end

CompiledEmbryos = UpdateAPAxisInfo(Prefix, CompiledEmbryos);




%%
close all
CurrentEmbryo = 1;
EmbryoImage = MembraneMat(:,:,CurrentEmbryo);
DisplayRange=[min(min(EmbryoImage)),max(max(EmbryoImage))];


EmbryoFigure=figure;
set(EmbryoFigure,'units', 'normalized', 'position',[0.01, .1, .7, .7]);

embryoAxes = axes(EmbryoFigure,'Units', 'normalized', 'Position', [0 0 1 1]);
cc=1;

% Show the first image
imEmbryo = imshow(EmbryoImage,DisplayRange,'Parent',embryoAxes);
hold on

AnteriorPole = plot(CompiledEmbryos.CoordAs(CurrentEmbryo,1), CompiledEmbryos.CoordAs(CurrentEmbryo,2),...
    'g.','MarkerSize',20);
PosteriorPole = plot(CompiledEmbryos.CoordPs(CurrentEmbryo,1), CompiledEmbryos.CoordPs(CurrentEmbryo,2),...
    'r.','MarkerSize',20);

DorsalEdge = plot(CompiledEmbryos.CoordDs(CurrentEmbryo,1), CompiledEmbryos.CoordDs(CurrentEmbryo,2),...
    'b.','MarkerSize',20);
VentralEdge = plot(CompiledEmbryos.CoordVs(CurrentEmbryo,1), CompiledEmbryos.CoordVs(CurrentEmbryo,2),...
    'y.','MarkerSize',20);
if MarkPV
PVAnteriorPole = plot(CompiledEmbryos.PVCoordAs(CurrentEmbryo,1), CompiledEmbryos.PVCoordAs(CurrentEmbryo,2),...
    'gp','MarkerFaceColor','g', 'MarkerSize',5);
PVPosteriorPole = plot(CompiledEmbryos.PVCoordPs(CurrentEmbryo,1), CompiledEmbryos.PVCoordPs(CurrentEmbryo,2),...
    'rp','MarkerFaceColor','r','MarkerSize',5);

PVDorsalEdge = plot(CompiledEmbryos.PVCoordDs(CurrentEmbryo,1), CompiledEmbryos.PVCoordDs(CurrentEmbryo,2),...
    'bp','MarkerSize',5,'MarkerFaceColor','b');
PVVentralEdge = plot(CompiledEmbryos.PVCoordVs(CurrentEmbryo,1), CompiledEmbryos.PVCoordVs(CurrentEmbryo,2),...
    'yp','MarkerSize',5,'MarkerFaceColor','y');
end
if CompiledEmbryos.CoordPs(CurrentEmbryo,1) > 0 & CompiledEmbryos.CoordAs(CurrentEmbryo,1) > 0
    Slope = (CompiledEmbryos.CoordPs(CurrentEmbryo,2)-CompiledEmbryos.CoordAs(CurrentEmbryo,2))/(CompiledEmbryos.CoordPs(CurrentEmbryo,1)-CompiledEmbryos.CoordAs(CurrentEmbryo,1));
    Intercept = CompiledEmbryos.CoordAs(CurrentEmbryo,2)-Slope*CompiledEmbryos.CoordAs(CurrentEmbryo,1);
    MidPoint = [(CompiledEmbryos.CoordAs(CurrentEmbryo,1)+CompiledEmbryos.CoordPs(CurrentEmbryo,1))/2, (CompiledEmbryos.CoordAs(CurrentEmbryo,2)+CompiledEmbryos.CoordPs(CurrentEmbryo,2))/2];
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
    EmbryoImage = MembraneMat(:,:,CurrentEmbryo);
    
    imEmbryo.CData = EmbryoImage;
    try
        caxis(embryoAxes, DisplayRange);
    end
    AnteriorPole.XData = CompiledEmbryos.CoordAs(CurrentEmbryo,1);
    AnteriorPole.YData = CompiledEmbryos.CoordAs(CurrentEmbryo,2);
    PosteriorPole.XData = CompiledEmbryos.CoordPs(CurrentEmbryo,1);
    PosteriorPole.YData = CompiledEmbryos.CoordPs(CurrentEmbryo,2);
    DorsalEdge.XData = CompiledEmbryos.CoordDs(CurrentEmbryo,1);
    DorsalEdge.YData = CompiledEmbryos.CoordDs(CurrentEmbryo,2);
    VentralEdge.XData = CompiledEmbryos.CoordVs(CurrentEmbryo,1);
    VentralEdge.YData = CompiledEmbryos.CoordVs(CurrentEmbryo,2);
    
    if MarkPV
    PVAnteriorPole.XData = CompiledEmbryos.PVCoordAs(CurrentEmbryo,1);
    PVAnteriorPole.YData = CompiledEmbryos.PVCoordAs(CurrentEmbryo,2);
    PVPosteriorPole.XData = CompiledEmbryos.PVCoordPs(CurrentEmbryo,1);
    PVPosteriorPole.YData = CompiledEmbryos.PVCoordPs(CurrentEmbryo,2);
    PVDorsalEdge.XData = CompiledEmbryos.PVCoordDs(CurrentEmbryo,1);
    PVDorsalEdge.YData = CompiledEmbryos.PVCoordDs(CurrentEmbryo,2);
    PVVentralEdge.XData = CompiledEmbryos.PVCoordVs(CurrentEmbryo,1);
    PVVentralEdge.YData = CompiledEmbryos.PVCoordVs(CurrentEmbryo,2);
    end
    %     imshow(imadjust(APImage), 'Parent', apAx)
    %imshow(APImage,DisplayRange)
    if CompiledEmbryos.CoordPs(CurrentEmbryo,1) > 0 & CompiledEmbryos.CoordAs(CurrentEmbryo,1) > 0
        Slope = (CompiledEmbryos.CoordPs(CurrentEmbryo,2)-CompiledEmbryos.CoordAs(CurrentEmbryo,2))/(CompiledEmbryos.CoordPs(CurrentEmbryo,1)-CompiledEmbryos.CoordAs(CurrentEmbryo,1));
        Intercept = CompiledEmbryos.CoordAs(CurrentEmbryo,2)-Slope*CompiledEmbryos.CoordAs(CurrentEmbryo,1);
        MidPoint = [(CompiledEmbryos.CoordAs(CurrentEmbryo,1)+CompiledEmbryos.CoordPs(CurrentEmbryo,1))/2, (CompiledEmbryos.CoordAs(CurrentEmbryo,2)+CompiledEmbryos.CoordPs(CurrentEmbryo,2))/2];
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
        CompiledEmbryos.CoordAs(CurrentEmbryo,:)=0;
        CompiledEmbryos.CoordPs(CurrentEmbryo,:)=0;
        CompiledEmbryos.CoordDs(CurrentEmbryo,:)=0;
        CompiledEmbryos.CoordVs(CurrentEmbryo,:)=0;
    elseif (ct~=0)&(cc=='a')	%Select anterior end
        [CoordAx,CoordAy]=ginputc(1,'Color',[1,1,1]);
        CompiledEmbryos.CoordAs(CurrentEmbryo,:) = [CoordAx,CoordAy];
    elseif (ct~=0)&(cc=='A') & MarkPV	%Select anterior end
        [CoordAx,CoordAy]=ginputc(1,'Color',[1,1,1]);
        CompiledEmbryos.PVCoordAs(CurrentEmbryo,:) = [CoordAx,CoordAy];
    elseif (ct~=0)&(cc=='p')    %Select posterior end
        [CoordPx,CoordPy]=ginputc(1,'Color',[1,1,1]);
        CompiledEmbryos.CoordPs(CurrentEmbryo,:) = [CoordPx,CoordPy];
    elseif (ct~=0)&(cc=='P')& MarkPV		%Select anterior end
        [CoordPx,CoordPy]=ginputc(1,'Color',[1,1,1]);
        CompiledEmbryos.PVCoordPs(CurrentEmbryo,:) = [CoordPx,CoordPy];
    elseif (ct~=0)&(cc=='d')	%Select dorsal edge
        [CoordDx,CoordDy]=ginputc(1,'Color',[1,1,1]);
        CompiledEmbryos.CoordDs(CurrentEmbryo,:) = [CoordDx,CoordDy];
    elseif (ct~=0)&(cc=='D') & MarkPV	%Select dorsal edge
        [CoordDx,CoordDy]=ginputc(1,'Color',[1,1,1]);
        CompiledEmbryos.PVCoordDs(CurrentEmbryo,:) = [CoordDx,CoordDy];
    elseif (ct~=0)&(cc=='v')    %Select ventral edge
        [CoordVx,CoordVy]=ginputc(1,'Color',[1,1,1]);
        CompiledEmbryos.CoordVs(CurrentEmbryo,:) = [CoordVx,CoordVy];
    elseif (ct~=0)&(cc=='V')& MarkPV	    %Select ventral edge
        [CoordVx,CoordVy]=ginputc(1,'Color',[1,1,1]);
        CompiledEmbryos.PVCoordVs(CurrentEmbryo,:) = [CoordVx,CoordVy];
    elseif (ct~=0)&(cc=='.')&(CurrentEmbryo < NEmbryos)    %Increase contrast
        CurrentEmbryo = CurrentEmbryo+1;
         EmbryoImage = MembraneMat(:,:,CurrentEmbryo);
        DisplayRange=[min(min(EmbryoImage)),max(max(EmbryoImage))];
    elseif (ct~=0)&(cc==',')&(CurrentEmbryo > 1)    %Increase contrast
        CurrentEmbryo = CurrentEmbryo-1;
         EmbryoImage = MembraneMat(:,:,CurrentEmbryo);
        DisplayRange=[min(min(EmbryoImage)),max(max(EmbryoImage))];
    elseif (ct~=0)&(cc=='>')&(CurrentEmbryo < NEmbryos)  
        NewIndex = find(1:NEmbryos > CurrentEmbryo & ~CompiledEmbryos.Checked, 1);
        if isempty(NewIndex)
            CurrentEmbryo = NEmbryos;
        else
            CurrentEmbryo = NewIndex;
        end
        EmbryoImage = MembraneMat(:,:,CurrentEmbryo);
        DisplayRange=[min(min(EmbryoImage)),max(max(EmbryoImage))];
    elseif (ct~=0)&(cc=='<')&(CurrentEmbryo > 1)  
        NewIndex = find(1:NEmbryos < CurrentEmbryo & ~CompiledEmbryos.Checked, 1,'last');
        if isempty(NewIndex)
            CurrentEmbryo = 1;
        else
            CurrentEmbryo = NewIndex;
        end
        EmbryoImage = MembraneMat(:,:,CurrentEmbryo);
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
        CompiledEmbryos = UpdateAPAxisInfo(Prefix, CompiledEmbryos);
        
        
        save(CompiledEmbryoPath,'CompiledEmbryos');
        disp('Axis information saved.')
        
        
    end
end
CompiledEmbryos = UpdateAPAxisInfo(Prefix, CompiledEmbryos);
save(CompiledEmbryoPath,'CompiledEmbryos');
disp('Axis information saved.')
close all
%%
if ~MarkPV
RotateEmbryoImages(Prefix);
end