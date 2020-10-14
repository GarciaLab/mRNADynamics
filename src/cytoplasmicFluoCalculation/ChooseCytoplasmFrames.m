function ChooseCytoplasmFrames(Prefix, varargin)
%%
% DESCRIPTION
%
% PARAMETERS
% Prefix: Prefix of the dataset being analyzed
%
% OPTIONS
%
% GUI COMMANDS
% .  - Move a frame forward
% ,  - Move a frame backwards
% >  - Move 5 frames forward
% <  - Move 5 frames backwards
% j  - Jump to a frame
% q  - Move a cycle forward
% w  - Move a cycle backwards
% a - Set current frame to first included in cytoplasmFrames from this NuclearCycle 
% b - Set current frame to last included in cytoplasmFrames from this
% nuclear cycle 
% s  - Save current analysis
% ~  - Create a different projection for the nuclear image
% m  - Increase contrast
% n  - Decrease contrast
% r  - Reset contrast setting
% x  - Exit and save
% 9  - Debug mode
%
%
% OUTPUT
% Ellipses.mat: saved to the folder 'Dropbox\Prefix\'
%
%
% Author (contact): uknown (hggarcia@berkeley.edu)
% Created: XXXX-XX-XX
% Last Updated: XXXX-XX-XX
% Documented by: Meghan Turner (meghan_turner@berkeley.edu)
%

cleanupObj = onCleanup(@myCleanupFun);




% for k = 1:length(varargin)
% end

liveExperiment = LiveExperiment(Prefix);

ProcPath = liveExperiment.userProcFolder;
DropboxFolder = liveExperiment.userResultsFolder;

Channel1 = liveExperiment.Channel1;
Channel2 = liveExperiment.Channel2;
Channel3 = liveExperiment.Channel3;
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

%Get the nuclei segmentation data
%Load the reference histogram for the fake histone channel

load('ReferenceHist.mat', 'ReferenceHist')
%% 
if exist([DropboxFolder, filesep, Prefix, filesep, 'cytoplasmFrames.mat'], 'file')
    load([DropboxFolder, filesep, Prefix, filesep, 'cytoplasmFrames.mat'])
else
    cytoplasmFrames = [];
end


Channels = {Channel1, Channel2, Channel3};


hisMat = getHisMat(liveExperiment);


nFrames = size(hisMat, 3);
%Get information about the image size
% HisImage=imread([PreProcPath,filesep,Prefix,filesep,D(1).name]);
HisImage = hisMat(:,:,1);
DisplayRange=[min(min(HisImage)),max(max(HisImage))];

nc = [];

%Make a vector containing the nc corresponding to each frame
for k=1:nFrames
    if k<nc9
        nc(k)=8;
    elseif (k>=nc9)&(k<nc10 || isnan(nc10))
        nc(k)=9;
    elseif (k>=nc10)&(k<nc11 || isnan(nc11))
        nc(k)=10;
    elseif (k>=nc11)&(k<nc12 || isnan(nc12))
        nc(k)=11;
    elseif (k>=nc12)& (k<nc13 || isnan(nc13))
        nc(k)=12;
    elseif (k>=nc13)&( k<nc14 || isnan(nc14) ) %#ok<*AND2>
        nc(k)=13;
    elseif k>=nc14
        nc(k)=14;
    end
end


ncCytoStartFrames = zeros(1, 14, 'uint16');
ncCytoEndFrames = zeros(1, 14, 'uint16');


%%
Overlay=figure;
set(Overlay,'units', 'normalized', 'position',[0.01, .5, .4, .4]);

overlayAxes = axes(Overlay,'Units', 'normalized', 'Position', [0 0 1 1]);

%%

tb = axtoolbar(overlayAxes);
tb.Visible = 'off';


CurrentFrame=1;
currentCharacter=1;

% Show the first image
imOverlay = imshow(HisImage,DisplayRange,'Border','Tight','Parent',overlayAxes);
projFlag = false;
set(0, 'CurrentFigure', Overlay)
CurrentNC = min(nc);
while (currentCharacter~='x')
    
    %Load subsequent images
    if ~projFlag
        HisImage = hisMat(:, :, CurrentFrame);
    else
        HisImage = Projection(:, :,CurrentFrame);
    end
    if ncCytoStartFrames(CurrentNC) == CurrentFrame
        HisOverlay=cat(3,mat2gray(HisImage),...
            mat2gray(HisImage)+2,...
            mat2gray(HisImage));
    elseif ncCytoEndFrames(CurrentNC) == CurrentFrame
        HisOverlay=cat(3,mat2gray(HisImage),...
            mat2gray(HisImage),...
            mat2gray(HisImage)+2);
    elseif (ncCytoStartFrames(CurrentNC) < CurrentFrame) & (ncCytoEndFrames(CurrentNC) > CurrentFrame)
        HisOverlay=cat(3,mat2gray(HisImage),...
            mat2gray(HisImage),...
            mat2gray(HisImage)+2);
    else 
        HisOverlay=cat(3,mat2gray(HisImage),...
            mat2gray(HisImage),...
            mat2gray(HisImage));
    end
    
    
    imOverlay.CData = HisOverlay;
    DisplayRange(2) = DisplayRange(2)+2;
    try
        caxis(overlayAxes, DisplayRange);
    end
    

    
    try
        FigureTitle=['Frame: ',num2str(CurrentFrame),'/',num2str(nFrames),...
            ', nc: ',num2str(nc(CurrentFrame))];
    catch
        FigureTitle=['Frame: ',num2str(CurrentFrame),'/',num2str(nFrames)];
    end
    
    set(Overlay,'Name',FigureTitle)
    

     
    
    
    
    %%
    
    tb = axtoolbar(overlayAxes);
    tb.Visible = 'off';

    
    ct=waitforbuttonpress;
    currentCharacter=get(Overlay,'currentcharacter');
    currentMouse=get(overlayAxes,'CurrentPoint');
    
    
    
    
    if (ct~=0)&(currentCharacter=='.')&(CurrentFrame<nFrames)
        CurrentFrame=CurrentFrame+1;
    elseif (ct~=0)&(currentCharacter==',')&(CurrentFrame>1)
        CurrentFrame=CurrentFrame-1;
    elseif (ct~=0)&(currentCharacter=='>')&(CurrentFrame+5<nFrames)
        CurrentFrame=CurrentFrame+5;
    elseif (ct~=0)&(currentCharacter=='<')&(CurrentFrame-4>1)
        CurrentFrame=CurrentFrame-5;
    elseif (ct~=0)&(currentCharacter=='s')
        save([DropboxFolder,filesep,Prefix,filesep,'cytoplasmFrames.mat'],'cytoplasmFrames', '-v6')
        disp('cytoplasmFrames saved.')
    elseif (ct~=0)&(currentCharacter=='a')
        ncCytoStartFrames(CurrentNC) = CurrentFrame;
    elseif (ct~=0)&(currentCharacter=='b')
        ncCytoEndFrames(CurrentNC) = CurrentFrame;
   
        
    elseif (ct~=0)&(currentCharacter=='j')
        iJump=input('Frame to jump to: ');
        if (floor(iJump)>0)&(iJump<=nFrames)
            CurrentFrame=iJump;
        else
            disp('Frame out of range.');
        end
        
    elseif (ct~=0)&(currentCharacter=='m')&(CurrentNC < max(nc))   %Change nuclear Cycle 
        CurrentNC = CurrentNC + 1;
        CurrentFrame = find(nc == CurrentNC, 1);
        
    elseif (ct~=0)&(currentCharacter=='n')&(CurrentNC-1 >= min(nc))     %Decrease contrast
        CurrentNC = CurrentNC - 1;
        CurrentFrame = find(nc == CurrentNC, 1);
        
    elseif (ct~=0)&(currentCharacter=='r')    %Reset the contrast
        DisplayRange=[min(min(HisImage)),max(max(HisImage))];
        
  
        
    elseif (ct~=0)&(currentCharacter=='~')
        
        ProjectionType = 'midsumprojection';
        
        movieMat = getMovieMat(liveExperiment); 
        [~, ~, Projection] = chooseNuclearChannels2(...
            movieMat, 'ProjectionType', ProjectionType,'Channels',...
            Channels,'ReferenceHist', ReferenceHist);
        
%         DisplayRange = [mean(mean(Projection(:, :, CurrentFrame))),...
%             max(max(Projection(:, :, CurrentFrame))) ];
            DisplayRange = [0,...
            max(max(Projection(:, :, CurrentFrame))) ];
        
        projFlag = true;
        
        disp('changed projection');
        
   
    elseif (ct~=0)&(currentCharacter=='q') %go to next nc
        nextncframes = find(nc == (nc(CurrentFrame)+1));
        if ~isempty(nextncframes)
            CurrentFrame = nextncframes(1);
        end
    elseif (ct~=0)&(currentCharacter=='w') %go to previous nc
        previousncframes = find(nc == (nc(CurrentFrame)-1));
        if ~isempty(previousncframes)
            CurrentFrame = previousncframes(1);
        end
  

    elseif (ct~=0)&(currentCharacter=='0')    %Debug mode
        keyboard
        
    end
end
close(Overlay)
%% 

cytoplasmFrames = [];
for nc=9:14
    if (nc==14) & (ncCytoEndFrames(nc) == 0)
        ncCytoEndFrames(nc) = nFrames;
    end
    if (ncCytoEndFrames(nc) > 0) & (ncCytoStartFrames(nc) > 0)
        cytoplasmFrames = [cytoplasmFrames ncCytoStartFrames(nc):ncCytoEndFrames(nc)];
    end
end


save([DropboxFolder,filesep,Prefix,filesep,'cytoplasmFrames.mat'],'cytoplasmFrames', '-v6')
fid = fopen([DropboxFolder,filesep,Prefix,filesep,'cytoplasmFrames.csv'],'w');
for i = 1:length(cytoplasmFrames)
    frame_string = num2str(cytoplasmFrames(i));
    if length(frame_string) == 1
        frame_string = ['00', frame_string];
    elseif length(frame_string) == 2
        frame_string = ['0', frame_string];
    end
     fprintf(fid,'%s,',frame_string);
end

fclose(fid)

writematrix([DropboxFolder,filesep,Prefix,filesep,'cytoplasmFrames.csv'], 'cytoplasmFrames');

end

