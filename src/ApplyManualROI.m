function ManualROI(Prefix, varargin)
%%  ManualROI(Prefix, varargin)
%
% DESCRIPTION
% This function allows for the manual specification of a region of
% interest to be included in CompiledParticles or CompiledNuclei. 
%
% ARGUMENTS
% varargin: A cell in which the first element is the prefix string of the data set
%           to analyze. Subsequent elements can be the options below.
%
% MANUAL: 
%
%c: clear all AP information
%s: select corners of the bounding rectangle for the ROI
%.: Increase contrast
%,: Decrease contrast
%r: Reset the contrast
%right click: Remove the previous selection 
%m: Manual stitching mode
%x: Exit and save

% Author (contact): Gabriella Martini (martini@berkeley.edu)
% Created: 11/29/19
% Last Updated: 11/29/19
%
% Documented by: Gabriella Martini (martini@berkeley.edu)
%%  Extract all Data
close all

%Load the folder information
[SourcePath, FISHPath, DefaultDropboxFolder, DropboxFolder, MS2CodePath, PreProcPath,...
configValues, movieDatabasePath] = DetermineAllLocalFolders(Prefix);

[Date, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution,...
Channel1, Channel2, Objective, Power, DataFolder, DropboxFolderName, Comments,...
nc9, nc10, nc11, nc12, nc13, nc14, CF] = getExperimentDataFromMovieDatabase(Prefix, DefaultDropboxFolder)


%This checks whether all ncs have been defined
ncCheck=[nc9,nc10,nc11,nc12,nc13,nc14];
if length(ncCheck)~=6
    error('Check the nc frames in the MovieDatabase entry. Some might be missing')
end

%Do we need to convert any NaN chars into doubles?
if strcmp(lower(nc14),'nan')
    nc14=nan;
end
if strcmp(lower(nc13),'nan')
    nc13=nan;
end
if strcmp(lower(nc12),'nan')
    nc12=nan;
end
if strcmp(lower(nc11),'nan')
    nc11=nan;
end
if strcmp(lower(nc10),'nan')
    nc10=nan;
end
if strcmp(lower(nc9),'nan')
    nc9=nan;
end


ncs=[nc9,nc10,nc11,nc12,nc13,nc14];

if (length(find(isnan(ncs)))==length(ncs))||(length(ncs)<6)
    error('Have the ncs been defined in MovieDatabase?')
end

% Load AP markers - This must be run after AP position is assigned. Add
% error message. 
% if exist([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'], 'file')
%     load([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'], 'coordA', 'coordP');
% end

%% Extract images 
DataFolder=[DropboxFolder,filesep,Prefix];
FilePrefix=[DataFolder(length(DropboxFolder)+2:end),'_'];

%% GUI

keeptracking=true;
ncexist=find(ncs~=0 & ~isnan(ncs)); % Logical array to tell if the Nuclear cycle exists
Current_Nuclear_Cycle=1;

% The user can choose to cycle between Nuclear Cycles or Frames, so these
% variables keep track of the Current Nuclear Cycle and the Current Frame
CurrentFrame=ncs(ncexist(Current_Nuclear_Cycle));%-1; # GM: Is this just done incorrectly? Why would it try to pull up frame 000 when there is no 000 frame 


while keeptracking
    displayFrame(CurrentFrame,[]);
    % Current_NC is only for display. It will not be used computationally
    [~,Current_Nuclear_Cycle]=min(abs(CurrentFrame-ncs(ncexist)));
    Current_NC=Current_Nuclear_Cycle+7+min(ncexist);
    display(Current_NC);
    axis image
    axis off
    title('Anterior (green), posterior (red); original')
    hold on
    if exist('coordA', 'var')
        plot(coordA(1),coordA(2),'g.','MarkerSize',20);
    end
    if exist('coordP', 'var')
        plot(coordP(1),coordP(2),'r.','MarkerSize',20);
    end
    if exist('coord1', 'var')
        plot(coord1(1),coord1(2),'.','MarkerSize',20);
    end
    if exist('coord2', 'var')
        plot(coord2(1),coord2(2),'c.','MarkerSize',20);
    end
    
    % This stores the variable to see if the current frame is the default
    % frame - it is used to toggle whether the user is moving between
    % frames
    frameReset=true;
    ct=waitforbuttonpress;
    cc=get(gcf,'currentcharacter');
    cm=get(gcf,'CurrentPoint');
    
    
    if (ct~=0)&(cc=='1')	%Select 1st corner of bounding rectangle
        [coord1x,coord1y]=ginputc(1,'Color',[1,1,1]);
        coord1 = [coord1x,coord1y];
    elseif (ct~=0)&(cc=='2')	%Select 2nd opposing corner of bounding rectangle
        [coord2x,coord2y]=ginputc(1,'Color',[0,1,1]);
        coord2 = [coord2x,coord2y];
    elseif (ct~=0)&(cc=='a')	%Select 2nd opposing corner of bounding rectangle
        [coordAx,coordAy]=ginputc(1,'Color',[0,1,1]);
        coordA = [coordAx,coordAy];
    elseif (ct~=0)&(cc=='p')	%Select 2nd opposing corner of bounding rectangle
        [coordPx,coordPy]=ginputc(1,'Color',[0,1,1]);
        coordP = [coordPx,coordPy];
    elseif cc=='x'
        keeptracking=false;
        % Moving between frames
    elseif cc=='.'
        CurrentFrame=CurrentFrame+1;
        frameReset=false;
    elseif cc=='j'
        jump=input('Which Frame do you want to jump to?');
        CurrentFrame=jump;
        frameReset=false;
    elseif cc=='<'
        CurrentFrame=CurrentFrame-5;
        frameReset=false;
    elseif cc=='>'
        CurrentFrame=CurrentFrame+5;
        frameReset=false;
    elseif cc==','
        CurrentFrame=CurrentFrame-1;
        frameReset=false;
    elseif cc=='q'
        CurrentFrame=CurrentFrame-1;
        frameReset=false;
    elseif cc=='t'
        Current_Nuclear_Cycle=Current_Nuclear_Cycle+1;
        if Current_Nuclear_Cycle>length(ncexist)
            display('Already in Last NC');
            Current_Nuclear_Cycle=length(ncexist);
        end
    elseif cc=='r'
        Current_Nuclear_Cycle=Current_Nuclear_Cycle-1;
        if Current_Nuclear_Cycle<1
            display('Already in First NC');
            Current_Nuclear_Cycle=1;
        end
    end
    
    if keeptracking && frameReset
        CurrentFrame=ncs(ncexist(Current_Nuclear_Cycle))-1;
    end
    CurrentFrame=max(1,CurrentFrame); 

end


%% Display
    function displayFrame(CurrentFrame,varargin)
        ImageHis=imread([PreProcPath,filesep,FilePrefix(1:end-1),filesep,...
            FilePrefix(1:end-1),'-His_',iIndex(CurrentFrame,3),'.tif']);
        %ImageHis=insertShape(ImageHis, 'circle', circles,'Color',colour)
        set(gcf,'Name',['Current Frame: ',num2str(CurrentFrame)]);
        imshow(ImageHis,'DisplayRange',[],'Border','Tight')
        %    set(gcf,'units', 'normalized', 'position',[.1   .55   .4   .35])
        set(gcf,'MenuBar','none','ToolBar','none') 
    end
close all;
end