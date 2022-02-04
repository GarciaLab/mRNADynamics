function StraightenEmbryoImages(Prefix)

close all
cleanupObj = onCleanup(@myCleanupFun);

liveExperiment = LiveExperiment(Prefix);

FrameInfo = getFrameInfo(liveExperiment);
FrameTimes_seconds = [FrameInfo(:).Time];
FrameTimes_minutes = FrameTimes_seconds/60;

ProcPath = liveExperiment.userProcFolder;
DropboxFolder = liveExperiment.userResultsFolder;

Channel1 = liveExperiment.Channel1;
anaphaseFrames = liveExperiment.anaphaseFrames;
nc9 = anaphaseFrames(1);
nc10 = anaphaseFrames(2);
nc11 = anaphaseFrames(3);
nc12 = anaphaseFrames(4);
nc13 = anaphaseFrames(5);
nc14 = anaphaseFrames(6);

xSize = liveExperiment.xDim;
ySize = liveExperiment.yDim;
PixelSize_um = liveExperiment.pixelSize_um;

load([DropboxFolder,filesep,Prefix,filesep,'EmbryoOrientationInfo.mat'])


membraneFile = [liveExperiment.preFolder, filesep,liveExperiment.Prefix, '-Membrane.tif'];

membraneMat = getMembraneMat(liveExperiment);
%%
rotatedMembraneMat = imrotate(membraneMat, -BoundaryAngleDegrees, 'bilinear', 'crop');

nFrames = size(membraneMat, 3);
%Get information about the image size
% HisImage=imread([PreProcPath,filesep,Prefix,filesep,D(1).name]);


%%
counter= 1;
CurrentFrame=nFrames;
MembraneImage = rotatedMembraneMat(:,:,CurrentFrame);
DisplayRange=[min(min(MembraneImage)),max(max(MembraneImage))];

Overlay=figure;
set(Overlay,'units', 'normalized', 'position',[0.01, .1, .8, .8]);

overlayAxes = axes(Overlay,'Units', 'normalized', 'Position', [0 0 1 1]);

tb = axtoolbar(overlayAxes);
tb.Visible = 'off';


currentCharacter=1;



% Show the first image
imOverlay = imshow(MembraneImage,DisplayRange,'Parent',overlayAxes);



hold(overlayAxes, 'on')


plot(overlayAxes, 0:xSize, BoundaryPosition*ones(1,xSize+1), 'b-');
rotatedMembraneFile = [liveExperiment.preFolder, filesep,liveExperiment.Prefix, '-Membrane_Rotated.tif'];

saveNuclearProjection(rotatedMembraneMat, rotatedMembraneFile)


