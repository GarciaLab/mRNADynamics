function CheckParticleTracking(varargin)

%The point of this function is to check the tracking of particles. The
%logic of this function should be similar to schnitzcells: We want to be
%able to correct both the segmentation and tracking.

%V5: Modified to handle Laurent's schnitzcells
%V4: Modified the code to use and correct the schnitzcells tracking


%Usage:

%Frame specific:
%. Move a frame forward
%, Move a frame backwards
%> Move five frames forward
%, Move five frames backwards
%a Move up in Z
%z Move down in z
%j Jump to a specified frame
%g,b Increase/decrease histone channel contrast



%Particle specific:
%m Move to the next particle
%n Move to the previous particle,
%k Jump to a specified particle
%c Connect two particle traces
%d Separate traces. If this is done on a particle with only one frame then
%  it disconnects it from its nucleus.
%q Cycle between approved status: green - approved; yellow - approved but
%  with conditions (drift of nucleus, for example)
%w Disapproove a trace
%p Identify a particle. It will also tell you the particle associated with
%  the clicked nucleus.
%e Approve/Disapproove  a frame within a trace
%u Move a particle detected with Threshold2 into the our structure.
%i Move a particle detected with Threshold2 into the our structure and
%  connect it to the current particle. This is a combination of "u" and
%  "c".


%Nuclear tracking specific:
%l Split a nucleus and select one or two daughter nuclei or stop the
%  lineage. Usage:
%       Click on one new nucleus + ENTER: Continue the schnitz with that nucleus.
%       Click on the current nucleus + ENTER: Split the schnitz. This time
%           point will be the first frame of the new schnitz.
%       Click on two nuclei: Split the current nucleus into two daughter
%       nuclei.
%       Click on the same nucleus twice: Split the current nucleus, but
%       with only one daughter nucleus.
%2 set parent of current nucleus
%p Find the particle associated with the clicked nucleus. It will also tell
%  you the closest particle associated you clicked on.
%9 check for nuclear tracking consistencies. This is useful while we're
%  getting the code to work well.
%1 give the nucleus number in the schnitzcell segmentation structure. This
%  only  works for troubleshooting and you need to be online and on the
%  Princeton network/VPN for now.


%General:
%t Show/hide particles from the second threshold
%s Save the current Particles structure
%h Show non-approved particles yellow or dissapproved particles
%y Input the frame/nc information again. This only works in the absence of
%  the histone channel
%r Reorder the particles according to initial frame
%f Redo tracking. It only gets done on the non-approved particles.
%o Zoom in/out around the particle's first frame.
%-/= Change the zoom factor when in zoom mode.
%0 Enter debug mode to fix things manually

close all


%% Information about about folders

[SourcePath,FISHPath,DefaultDropboxFolder,MS2CodePath,SchnitzcellsFolder]=...
    DetermineLocalFolders;

%Also get the computer name. We'll use this later

%Find out which computer this is. That will determine the folder structure.
[ret, name] = system('hostname');  
if ret ~= 0,  
   if ispc  
      name = getenv('COMPUTERNAME');  
   else  
      name = getenv('HOSTNAME');  
   end  
end  
name = lower(name); 


if isempty(varargin)
    DataFolder=uigetdir(DefaultDropboxFolder,'Select data set to analyze');
else
    [SourcePath,FISHPath,DropboxFolder,MS2CodePath,SchnitzcellsFolder]=...
        DetermineLocalFolders(varargin{1});
    DataFolder=[DropboxFolder,filesep,varargin{1}];
end
    
%Flag to sort or not particles according to their starting frame
NoSort=0;
%Flag to just save the data. This is good for CompileAll
ForCompileAll=0;
if length(varargin)>1
    for i=2:length(varargin)
        if strcmp(varargin{i},'NoSort')
            NoSort=1;
        elseif strcmp(varargin{i},'ForCompileAll')
            ForCompileAll=1;
        end
    end
end

% ES 2014-01-07: Allows for the code to automatically accept defaults
% without requiring user input
if sum(ismember(varargin, 'AutoMode'))
    AutoFlag = 1;
else
    AutoFlag = 0;
end

    
%%


FilePrefix=[DataFolder(length(DropboxFolder)+2:end),'_'];

%Now get the actual folders
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,SchnitzcellsFolder]=...
    DetermineLocalFolders(FilePrefix(1:end-1));

load([DataFolder,filesep,'Particles.mat'])

%Check that FrameInfo exists
if exist([DataFolder,filesep,'FrameInfo.mat'])
    load([DataFolder,filesep,'FrameInfo.mat'])
else
    warning('No FrameInfo.mat found. Trying to continue')
    %Adding frame information
    DHis=dir([FISHPath,filesep,'Data',filesep,FilePrefix(1:end-1),filesep,'*His*.tif']);
    FrameInfo(length(DHis)).nc=[];
    %Adding information

    Dz=dir([FISHPath,filesep,'Data',filesep,FilePrefix(1:end-1),filesep,FilePrefix(1:end-1),'*001*.tif']);
    NumberSlices=length(Dz)-1;
    
    for i=1:length(FrameInfo)
        FrameInfo(i).NumberSlices=NumberSlices;
    end
    
end


%Some parameters:
SnippetSize=13;     %Size of the snippets generated by Michael's code


%Add FramesApproved where necessary
if ~isfield(Particles,'FrameApproved')
    for i=1:length(Particles)
        Particles(i).FrameApproved=logical(ones(size(Particles(i).Frame)));
    end
else
    for i=1:length(Particles)
        if isempty(Particles(i).FrameApproved)
            Particles(i).FrameApproved=logical(ones(size(Particles(i).Frame)));
        end
    end
end



%Check if we have the histone channel and we have done the nuclear
%segmentation.
if exist([FISHPath,filesep,'Data',filesep,FilePrefix(1:end-1),filesep,...
        FilePrefix(1:end-1),'-His_',iIndex(1,3),'.tif'])|...
        exist([FISHPath,filesep,'Data',filesep,FilePrefix(1:end-1),filesep,...
        FilePrefix(1:end-1),'_His_',iIndex(1,3),'.tif'])
    load([DropboxFolder,filesep,FilePrefix(1:end-1),filesep,'Ellipses.mat'])
    UseHistoneOverlay=1;
else
    UseHistoneOverlay=0;
end


%Check that we have the nuclear tracking done using schnitzcells
if exist([DropboxFolder,filesep,FilePrefix(1:end-1),filesep,FilePrefix(1:end-1),'_lin.mat'])
    UseSchnitz=1;
else
    UseSchnitz=0;
end
    





%Determine division times
%Load the information about the nc from the XLS file
[Num,Txt,XLSRaw]=xlsread([DefaultDropboxFolder,filesep,'MovieDatabase.xlsx']);
XLSHeaders=Txt(1,:);
Txt=Txt(2:end,:);

ExperimentTypeColumn=find(strcmp(XLSRaw(1,:),'ExperimentType'));
ExperimentAxisColumn=find(strcmp(XLSRaw(1,:),'ExperimentAxis'));

DataFolderColumn=find(strcmp(XLSRaw(1,:),'DataFolder'));
Dashes=findstr(FilePrefix,'-');
PrefixRow=find(strcmp(XLSRaw(:,DataFolderColumn),[FilePrefix(1:Dashes(3)-1),'\',FilePrefix(Dashes(3)+1:end-1)]));

ExperimentType=XLSRaw{PrefixRow,ExperimentTypeColumn};
ExperimentAxis=XLSRaw{PrefixRow,ExperimentAxisColumn};

%Find the different columns.
DataFolderColumn=find(strcmp(XLSRaw(1,:),'DataFolder'));
nc9Column=find(strcmp(XLSRaw(1,:),'nc9'));
nc10Column=find(strcmp(XLSRaw(1,:),'nc10'));
nc11Column=find(strcmp(XLSRaw(1,:),'nc11'));
nc12Column=find(strcmp(XLSRaw(1,:),'nc12'));
nc13Column=find(strcmp(XLSRaw(1,:),'nc13'));
nc14Column=find(strcmp(XLSRaw(1,:),'nc14'));
CFColumn=find(strcmp(XLSRaw(1,:),'CF'));
Channel2Column=find(strcmp(XLSRaw(1,:),'Channel2'));


%Find the corresponding entry in the XLS file
if (~isempty(findstr(FilePrefix,'Bcd')))&(isempty(findstr(FilePrefix,'BcdE1')))&...
        (isempty(findstr(FilePrefix(1:end-1),'NoBcd')))
    warning('This step in CheckParticleTracking will most likely have to be modified to work')
    XLSEntry=find(strcmp(XLSRaw(:,DataFolderColumn),...
        [Date,'\BcdGFP-HisRFP']));
else
    XLSEntry=find(strcmp(XLSRaw(:,DataFolderColumn),...
        [FilePrefix(1:Dashes(3)-1),'\',FilePrefix(Dashes(3)+1:end-1)]));
end


if strcmp(XLSRaw(XLSEntry,Channel2Column),'His-RFP')
    nc9=XLSRaw{XLSEntry,nc9Column};
    nc10=XLSRaw{XLSEntry,nc10Column};
    nc11=XLSRaw{XLSEntry,nc11Column};
    nc12=XLSRaw{XLSEntry,nc12Column};
    nc13=XLSRaw{XLSEntry,nc13Column};
    nc14=XLSRaw{XLSEntry,nc14Column};
    CF=XLSRaw{XLSEntry,CFColumn};
    
    
    for i=1:length(FrameInfo)
        if i<nc9
            FrameInfo(i).nc=8;
        elseif (i>=nc9)&(i<nc10)
            FrameInfo(i).nc=9;
        elseif (i>=nc10)&(i<nc11)
            FrameInfo(i).nc=10;
        elseif (i>=nc11)&(i<=nc12)
            FrameInfo(i).nc=11;
        elseif (i>=nc12)&(i<=nc13)
            FrameInfo(i).nc=12;
        elseif (i>=nc13)&(i<=nc14)
            FrameInfo(i).nc=13;
        elseif i>=nc14
            FrameInfo(i).nc=14;
        end
    end
else
    error('nc information not define in MovieDatabase.xlsx')
end




% 
% %Find the different columns.
% DataFolderColumn=7;
% 
% %Convert the prefix into the string used in the XLS file
% Dashes=findstr(FilePrefix,'-');
% 
% %Find the corresponding entry in the XLS file
% XLSEntry=find(strcmp(Txt(:,DataFolderColumn),...
%     ['D:\Hernan\LivemRNA\Analysis\MS2\MCPNoNLS+MS2\',...
%     FilePrefix(1:Dashes(3)-1),filesep,FilePrefix(Dashes(3)+1:end-1)]));
% 
% if isempty(XLSEntry)
%     warning('The name format is a mess. I had to do this for KITP')
%     XLSEntry=find(strcmp(Txt(:,DataFolderColumn),...
%         ['D:\Hernan\LivemRNA\Analysis\MS2\MCPNoNLS+MS2\2013-',FilePrefix(1:Dashes(2)-1),...
%         filesep,FilePrefix(Dashes(2)+1:end-1)]));
% end
% 
% 
% if strcmp(Txt(XLSEntry,4),'Yes')
% 
%     nc9=Num(XLSEntry,5);
%     nc10=Num(XLSEntry,6);
%     nc11=Num(XLSEntry,7);
%     nc12=Num(XLSEntry,8);
%     nc13=Num(XLSEntry,9);
%     nc14=Num(XLSEntry,10);
%     CF=Num(XLSEntry,11);
% 
% 
%     for i=1:length(FrameInfo)
%         if i<nc9
%             FrameInfo(i).nc=8;
%         elseif (i>=nc9)&(i<nc10)
%             FrameInfo(i).nc=9;
%         elseif (i>=nc10)&(i<nc11)
%             FrameInfo(i).nc=10;
%         elseif (i>=nc11)&(i<=nc12)
%             FrameInfo(i).nc=11;
%         elseif (i>=nc12)&(i<=nc13)
%             FrameInfo(i).nc=12;
%         elseif (i>=nc13)&(i<=nc14)
%             FrameInfo(i).nc=13;
%         elseif i>=nc14
%             FrameInfo(i).nc=14;
%         end
%     end
% else
%     display('No nc information provided in the XLS file.')
% end



%Check if we have already determined nc
if (~isfield(FrameInfo,'nc'))&(~UseHistoneOverlay)
    FrameInfo=DetermineNC(fad,Particles,FrameInfo);
elseif UseSchnitz
    load([DropboxFolder,filesep,FilePrefix(1:end-1),filesep,FilePrefix(1:end-1),'_lin.mat'])
    
    %Remove the schnitz fields that can give us problems potentially if
    %present. I don't know how this came to be, but it's for fields that
    %are not all that relevant. The fields are: approved, ang
    if isfield(schnitzcells,'approved')
        schnitzcells=rmfield(schnitzcells,'approved');
    end
    if isfield(schnitzcells,'ang')
        schnitzcells=rmfield(schnitzcells,'ang');
    end
end


%Order particles by the earliest frame they appear at. This makes the
%tracking a lot easier!
if ~NoSort
    for i=1:length(Particles)
        FirstFrame(i)=Particles(i).Frame(1);
    end
    [Dummy,Permutations]=sort(FirstFrame);
    Particles=Particles(Permutations);
end

% HG: Do this later like compileschnitz
% 
% %Now, get the particle positions (if they're not there already). Notice
% %that the code pulls out the position information from fad. This is because
% %of historical reasons mostly.
% if ~isfield(Particles,'xPos')
%     for i=1:length(Particles)
%         for j=1:length(Particles(i).Frame)
%             [x,y]=fad2xyzFit(Particles(i).Frame(j),fad, 'addMargin'); 
%             Particles(i).xPos(j)=x(Particles(i).Index(j));
%             Particles(i).yPos(j)=y(Particles(i).Index(j));
%         end
%     end
% end




TotalFrames=length(fad.channels);


%Some flags and initial parameters
ShowThreshold2=1;
HideApprovedFlag=0;
ParticleToFollow=[];
ZSlices=FrameInfo(1).NumberSlices+2;   %Note that the blank slices are included
CurrentZ=round(ZSlices/2);          
ManualZFlag=0;
CurrentParticle=1;
PreviousParticle=1;
CurrentFrameWithinParticle=1;


CurrentFrame=Particles(1).Frame(1);


DisplayRange=[];

ZoomMode=0;
ZoomRange=50;


Overlay=figure;
SnippetFig=figure;
ZProfileFig=figure;
TraceFig=figure;

if UseHistoneOverlay
    HisOverlayFig=figure;
end


cc=1;

%See if we just want to save the data
if ForCompileAll
    cc=13;
end

while (cc~=13)
    EllipseHandle=[];
    EllipseHandleYellow=[];
    EllipseHandleBlue=[];
    EllipseHandleWhite=[];
    EllipseHandleGreen=[];

    
    %Get the coordinates taking the margins into account
    [x,y]=fad2xyzFit(CurrentFrame,fad, 'addMargin'); 
    
    %If the approved field does not exist create it
    if ~isfield(Particles,'Approved')
        for i=1:length(Particles)
            Particles(i).Approved=0;
        end
    end
    ApprovedParticles=[Particles.Approved];
    
    %Pull out the right particle if it exists in this frame
    CurrentParticleIndex=Particles(CurrentParticle).Index(Particles(CurrentParticle).Frame==CurrentFrame);
    
    %This is the position of the current particle
    xTrace=x(CurrentParticleIndex);
    yTrace=y(CurrentParticleIndex);
    
    %These are the positions of all the approved and disapprooved particles
    %Find the particles in this frame
    IndexApprovedParticles=[];
    for i=1:length(Particles)
        if sum(Particles(i).Frame==CurrentFrame)&sum(Particles(i).Approved==1)
            IndexApprovedParticles=[IndexApprovedParticles,...
                Particles(i).Index(Particles(i).Frame==CurrentFrame)];
        end
    end
    xApproved=x(IndexApprovedParticles);
    yApproved=y(IndexApprovedParticles);

    IndexDisapprovedParticles=[];
    for i=1:length(Particles)
        if sum(Particles(i).Frame==CurrentFrame)&sum(Particles(i).Approved==-1)
            IndexDisapprovedParticles=[IndexDisapprovedParticles,...
                Particles(i).Index(Particles(i).Frame==CurrentFrame)];
        end
    end
    xDisapproved=x(IndexDisapprovedParticles);
    yDisapproved=y(IndexDisapprovedParticles);
    
    
    if (~isempty(xTrace))&(~ManualZFlag)
        CurrentZ=...
            fad.channels(CurrentFrame).fits.z(CurrentParticleIndex);
        ManualZFlag=0;
    end
    
    
    try
        Image=imread([FISHPath,filesep,'Data',filesep,FilePrefix(1:end-1),filesep,...
            FilePrefix,iIndex(CurrentFrame,3),'_z',iIndex(CurrentZ,2),'.tif']);
    catch
        display(['Warning: Could not load file: ',...
            FilePrefix,iIndex(CurrentFrame,3),'_z',iIndex(CurrentZ,2),'.tif'])
    end
    

    figure(Overlay)
    imshow(Image,[],'Border','Tight')
    hold on
    plot(x,y,'or')
    plot(xApproved,yApproved,'ob')
    plot(xTrace,yTrace,'og')
    plot(xDisapproved,yDisapproved,'^r')
    hold off
    if ~UseHistoneOverlay
        set(gcf,'Position',[10    89   677   599]);
        set(gcf,'MenuBar','none','ToolBar','none')
    else
        if strcmp('albert-pc',name(1:end-1))
            set(gcf,'Position',[7   513   835   532]);
            % halfscreen [7   513   835   532]
            % fullscreen [9         594        1106         451]
        else
            set(gcf,'Position',[10   414   677   351]);
        end
        set(gcf,'MenuBar','none','ToolBar','none')
    end
    
    set(gcf,'Name',['Particle: ',num2str(CurrentParticle),'/',num2str(length(Particles)),...
        ', Frame: ',num2str(CurrentFrame),'/',num2str(TotalFrames),...
        ', Z: ',num2str(CurrentZ),'/',num2str(ZSlices),' nc: ', num2str(FrameInfo(CurrentFrame).nc)])
    
    if UseSchnitz
        %Show all the nuclei
        hold on       
        EllipseHandle=notEllipse(Ellipses{CurrentFrame}(:,3),...
            Ellipses{CurrentFrame}(:,4),...
            Ellipses{CurrentFrame}(:,5),...
            Ellipses{CurrentFrame}(:,1)+1,...
            Ellipses{CurrentFrame}(:,2)+1,'r',50);
        hold off
            
%         hold on
%         [NEllipses,Dummy]=size(Ellipses{CurrentFrame});
%         for i=1:NEllipses
%             EllipseHandle=[EllipseHandle,ellipse(Ellipses{CurrentFrame}(i,3),...
%                 Ellipses{CurrentFrame}(i,4),...
%                 Ellipses{CurrentFrame}(i,5),...
%                 Ellipses{CurrentFrame}(i,1)+1,...
%                 Ellipses{CurrentFrame}(i,2)+1)];
%         end
%         set(EllipseHandle,'Color','r')
%         hold off
        
        
        %Show the ones that have been approved
        
        hold on
        schnitzCellNo=[];
        for i=1:length(Particles)
            if Particles(i).Approved==1
                schnitzIndex=find((schnitzcells(Particles(i).Nucleus).frames)==CurrentFrame);
                schnitzCellNo=[schnitzCellNo,schnitzcells(Particles(i).Nucleus).cellno(schnitzIndex)];
            end
        end
        
        EllipseHandleBlue=notEllipse(Ellipses{CurrentFrame}(schnitzCellNo,3),...
            Ellipses{CurrentFrame}(schnitzCellNo,4),...
            Ellipses{CurrentFrame}(schnitzCellNo,5),...
            Ellipses{CurrentFrame}(schnitzCellNo,1)+1,...
            Ellipses{CurrentFrame}(schnitzCellNo,2)+1,'b',50);
        
        
        hold off
        
        
        %Show the corresponding nucleus
        if ~isempty(Particles(CurrentParticle).Nucleus)
            SchnitzIndex=find(schnitzcells(Particles(CurrentParticle).Nucleus).frames==CurrentFrame);
            NucleusIndex=schnitzcells(Particles(CurrentParticle).Nucleus).cellno(SchnitzIndex);

            if ~isempty(NucleusIndex)
                hold on
                EllipseHandleGreen=ellipse(Ellipses{CurrentFrame}(NucleusIndex,3),...
                    Ellipses{CurrentFrame}(NucleusIndex,4),...
                    Ellipses{CurrentFrame}(NucleusIndex,5),...
                    Ellipses{CurrentFrame}(NucleusIndex,1)+1,...
                    Ellipses{CurrentFrame}(NucleusIndex,2)+1);
                set(EllipseHandleGreen,'Color','g')
                hold off
            else
                %('Error: Particle without an associated nucleus?')
            end
        
        
        
            %Show the daughter nuclei if applicable
            DaughterE=schnitzcells(Particles(CurrentParticle).Nucleus).E;
            DaughterD=schnitzcells(Particles(CurrentParticle).Nucleus).D;

            
            if DaughterE~=0
                SchnitzIndex=find(schnitzcells(DaughterE).frames==CurrentFrame);
                NucleusIndex=schnitzcells(DaughterE).cellno(SchnitzIndex);

                if ~isempty(NucleusIndex)
                    hold on
                    EllipseHandleWhite=[EllipseHandleWhite,ellipse(Ellipses{CurrentFrame}(NucleusIndex,3),...
                        Ellipses{CurrentFrame}(NucleusIndex,4),...
                        Ellipses{CurrentFrame}(NucleusIndex,5),...
                        Ellipses{CurrentFrame}(NucleusIndex,1)+1,...
                        Ellipses{CurrentFrame}(NucleusIndex,2)+1)];
                    
                    hold off
                else
                    %('Error: Particle without an associated nucleus?')
                end
            end

            if DaughterD~=0
                SchnitzIndex=find(schnitzcells(DaughterD).frames==CurrentFrame);
                NucleusIndex=schnitzcells(DaughterD).cellno(SchnitzIndex);

                if ~isempty(NucleusIndex)
                    hold on
                    EllipseHandleWhite=[EllipseHandleWhite,ellipse(Ellipses{CurrentFrame}(NucleusIndex,3),...
                        Ellipses{CurrentFrame}(NucleusIndex,4),...
                        Ellipses{CurrentFrame}(NucleusIndex,5),...
                        Ellipses{CurrentFrame}(NucleusIndex,1)+1,...
                        Ellipses{CurrentFrame}(NucleusIndex,2)+1)];
                    hold off
                else
                    %('Error: Particle without an associated nucleus?')
                end
            end
            
            if ~isempty(EllipseHandleWhite)
                set(EllipseHandleWhite,'Color','w')
            end
            
            %Show the mother nucleus if applicable
            Mother=schnitzcells(Particles(CurrentParticle).Nucleus).P;

            if Mother~=0
                SchnitzIndex=find(schnitzcells(Mother).frames==CurrentFrame);
                NucleusIndex=schnitzcells(Mother).cellno(SchnitzIndex);

                if ~isempty(NucleusIndex)
                    hold on
                    EllipseHandleYellow=ellipse(Ellipses{CurrentFrame}(NucleusIndex,3),...
                        Ellipses{CurrentFrame}(NucleusIndex,4),...
                        Ellipses{CurrentFrame}(NucleusIndex,5),...
                        Ellipses{CurrentFrame}(NucleusIndex,1)+1,...
                        Ellipses{CurrentFrame}(NucleusIndex,2)+1);
                    set(EllipseHandleYellow,'Color','y')
                    hold off
                else
                    %('Error: Particle without an associated nucleus?')
                end
            end
            
        else
            warning('Warning: This particle does not have an associated nucleus')
        end
        
            
    end
    

    
    if ApprovedParticles(CurrentParticle)==1
        set(gcf,'Color','g')
    elseif ApprovedParticles(CurrentParticle)==-1
        set(gcf,'Color','r')
    elseif ApprovedParticles(CurrentParticle)==2
        set(gcf,'Color','y')
    else
        set(gcf,'Color','default')
    end
    
    if ShowThreshold2
        [x2,y2]=fad2xyzFit(CurrentFrame,fad2, 'addMargin'); 
        hold on
        plot(x2,y2,'sr')
        hold off
    end
    
    if ZoomMode
        %Find the closest frame
        [Dummy,MinIndex]=min((Particles(CurrentParticle).Frame-CurrentFrame).^2);
        if length(MinIndex)>1
            MinIndex=MinIndex(1);
        end
        [xForZoom,yForZoom]=fad2xyzFit(Particles(CurrentParticle).Frame(MinIndex),fad, 'addMargin'); 
        xForZoom=xForZoom(Particles(CurrentParticle).Index(MinIndex));
        yForZoom=yForZoom(Particles(CurrentParticle).Index(MinIndex));
       
 
        xlim([xForZoom-ZoomRange,xForZoom+ZoomRange])
        ylim([yForZoom-ZoomRange/2,yForZoom+ZoomRange/2])
    end
    

    if UseHistoneOverlay
        figure(HisOverlayFig)
        
        try
            ImageHis=imread([FISHPath,filesep,'Data',filesep,FilePrefix(1:end-1),filesep,...
                FilePrefix(1:end-1),'-His_',iIndex(CurrentFrame,3),'.tif']);
        catch %Had to do this for KITP
            ImageHis=imread([FISHPath,filesep,'Data',filesep,FilePrefix(1:end-1),filesep,...
                FilePrefix(1:end-1),'_His_',iIndex(CurrentFrame,3),'.tif']);
        end

        if isempty(DisplayRange)
            HisOverlayImage=cat(3,mat2gray(ImageHis)*2,mat2gray(Image),zeros(size(Image)));
        else
            HisOverlayImage=cat(3,mat2gray(ImageHis,double(DisplayRange))*2,mat2gray(Image),zeros(size(Image)));
        end
        imshow(HisOverlayImage,'Border','Tight')
        if strcmp('albert-pc',name(1:end-1))
            set(gcf,'Position',[6    88   691   379])
            %Fullscreen [8         110        1107         440]
        else
            set(gcf,'Position',[-7    60   691   311])
        end
        set(gcf,'MenuBar','none','ToolBar','none')
        
       
        hold on
        plot(x,y,'ow')
        plot(xApproved,yApproved,'ob')
        plot(xTrace,yTrace,'og')
        hold off
        
        if ShowThreshold2
            [x2,y2]=fad2xyzFit(CurrentFrame,fad2, 'addMargin'); 
            hold on
            plot(x2,y2,'sw')
            hold off
        end
  
        
        if UseSchnitz
            
            copyobj(EllipseHandle,gca)
            copyobj(EllipseHandleBlue,gca)
            copyobj(EllipseHandleGreen,gca)
            copyobj(EllipseHandleWhite,gca)
            copyobj(EllipseHandleYellow,gca)

            
            
            
%             %Show all the nuclei
%             hold on
%             [NEllipses,Dummy]=size(Ellipses{CurrentFrame});
%             for i=1:NEllipses
%                 EllipseHandle=ellipse(Ellipses{CurrentFrame}(i,3),...
%                     Ellipses{CurrentFrame}(i,4),...
%                     Ellipses{CurrentFrame}(i,5),...
%                     Ellipses{CurrentFrame}(i,1)+1,...
%                     Ellipses{CurrentFrame}(i,2)+1);
%                 set(EllipseHandle,'Color','r')
% 
%             end
%             hold off
%             
%             
%             %Show the ones that have been approved
%             hold on
%             for i=1:length(Particles)
%                 if Particles(i).Approved==1
%                     schnitzIndex=find((schnitzcells(Particles(i).Nucleus).frames-1)==CurrentFrame);
%                     schnitzCellNo=schnitzcells(Particles(i).Nucleus).cellno(schnitzIndex);
%                     EllipseHandle=ellipse(Ellipses{CurrentFrame}(schnitzCellNo,3),...
%                         Ellipses{CurrentFrame}(schnitzCellNo,4),...
%                         Ellipses{CurrentFrame}(schnitzCellNo,5),...
%                         Ellipses{CurrentFrame}(schnitzCellNo,1)+1,...
%                         Ellipses{CurrentFrame}(schnitzCellNo,2)+1);
%                     set(EllipseHandle,'Color','b')
%                 end
%             end
%             hold off
%                 
%             
%             %Show the corresponding nucleus
%             if ~isempty(Particles(CurrentParticle).Nucleus)
%                 SchnitzIndex=find(schnitzcells(Particles(CurrentParticle).Nucleus).frames==CurrentFrame+1);
%                 NucleusIndex=schnitzcells(Particles(CurrentParticle).Nucleus).cellno(SchnitzIndex);
% 
%                 if ~isempty(NucleusIndex)
%                     hold on
%                     EllipseHandle=ellipse(Ellipses{CurrentFrame}(NucleusIndex,3),...
%                         Ellipses{CurrentFrame}(NucleusIndex,4),...
%                         Ellipses{CurrentFrame}(NucleusIndex,5),...
%                         Ellipses{CurrentFrame}(NucleusIndex,1)+1,...
%                         Ellipses{CurrentFrame}(NucleusIndex,2)+1);
%                     set(EllipseHandle,'Color','g')
%                     hold off
%                 else
%                     %('Error: Particle without an associated nucleus?')
%                 end
% 
%                 %Show the daughter nuclei if applicable
%                 DaughterE=schnitzcells(Particles(CurrentParticle).Nucleus).E;
%                 DaughterD=schnitzcells(Particles(CurrentParticle).Nucleus).D;
% 
%                 if DaughterE~=0
%                     SchnitzIndex=find(schnitzcells(DaughterE).frames==CurrentFrame+1);
%                     NucleusIndex=schnitzcells(DaughterE).cellno(SchnitzIndex);
% 
%                     if ~isempty(NucleusIndex)
%                         hold on
%                         EllipseHandle=ellipse(Ellipses{CurrentFrame}(NucleusIndex,3),...
%                             Ellipses{CurrentFrame}(NucleusIndex,4),...
%                             Ellipses{CurrentFrame}(NucleusIndex,5),...
%                             Ellipses{CurrentFrame}(NucleusIndex,1)+1,...
%                             Ellipses{CurrentFrame}(NucleusIndex,2)+1);
%                         set(EllipseHandle,'Color','w')
%                         hold off
%                     else
%                         %('Error: Particle without an associated nucleus?')
%                     end
%                 end
% 
%                 if DaughterD~=0
%                     SchnitzIndex=find(schnitzcells(DaughterD).frames==CurrentFrame+1);
%                     NucleusIndex=schnitzcells(DaughterD).cellno(SchnitzIndex);
% 
%                     if ~isempty(NucleusIndex)
%                         hold on
%                         EllipseHandle=ellipse(Ellipses{CurrentFrame}(NucleusIndex,3),...
%                             Ellipses{CurrentFrame}(NucleusIndex,4),...
%                             Ellipses{CurrentFrame}(NucleusIndex,5),...
%                             Ellipses{CurrentFrame}(NucleusIndex,1)+1,...
%                             Ellipses{CurrentFrame}(NucleusIndex,2)+1);
%                         set(EllipseHandle,'Color','w')
%                         hold off
%                     else
%                         %('Error: Particle without an associated nucleus?')
%                     end
%                 end
%                 
%                 %Show the mother nucleus if applicable
%                 Mother=schnitzcells(Particles(CurrentParticle).Nucleus).P;
% 
%                 if Mother~=0
%                     SchnitzIndex=find(schnitzcells(Mother).frames==CurrentFrame+1);
%                     NucleusIndex=schnitzcells(Mother).cellno(SchnitzIndex);
% 
%                     if ~isempty(NucleusIndex)
%                         hold on
%                         EllipseHandle=ellipse(Ellipses{CurrentFrame}(NucleusIndex,3),...
%                             Ellipses{CurrentFrame}(NucleusIndex,4),...
%                             Ellipses{CurrentFrame}(NucleusIndex,5),...
%                             Ellipses{CurrentFrame}(NucleusIndex,1)+1,...
%                             Ellipses{CurrentFrame}(NucleusIndex,2)+1);
%                         set(EllipseHandle,'Color','y')
%                         hold off
%                     else
%                         %('Error: Particle without an associated nucleus?')
%                     end
%                 end
%                 
%             else
%                 warning('Warning: This particle does not have an associated nucleus')
%             end

        end
    
        set(gcf,'Name',['Particle: ',num2str(CurrentParticle),'/',num2str(length(Particles)),...
            ', Frame: ',num2str(CurrentFrame),'/',num2str(TotalFrames),...
            ', Z: ',num2str(CurrentZ),'/',num2str(ZSlices),' nc: ', num2str(FrameInfo(CurrentFrame).nc)])
        
        if ZoomMode
            %Find the closest frame
            [Dummy,MinIndex]=min((Particles(CurrentParticle).Frame-CurrentFrame).^2);
            if length(MinIndex)>1
                MinIndex=MinIndex(1);
            end
            [xForZoom,yForZoom]=fad2xyzFit(Particles(CurrentParticle).Frame(MinIndex),fad, 'addMargin'); 
            xForZoom=xForZoom(Particles(CurrentParticle).Index(MinIndex));
            yForZoom=yForZoom(Particles(CurrentParticle).Index(MinIndex));


            xlim([xForZoom-ZoomRange,xForZoom+ZoomRange])
            ylim([yForZoom-ZoomRange/2,yForZoom+ZoomRange/2])
        end
    end
    
    

    
    
    
    figure(SnippetFig)
    if (~isempty(xTrace))
        
        CurrentSnippet=mat2gray(fad.channels(CurrentFrame).fits.snippets(:,:,CurrentParticleIndex));
        IntegrationArea=bwperim(fad.channels(CurrentFrame).fits.maskUsedForTotalInt);
        SnippetOverlay=cat(3,IntegrationArea/2 + ...
            +CurrentSnippet,CurrentSnippet,CurrentSnippet);
        
    
        imshow(SnippetOverlay,...
            [],'Border','Tight','InitialMagnification',1000)
        if strcmp('albert-pc',name(1:end-1))
            set(gcf,'Position',[861   834   303   152])
            % fullscreen [1133          46         395         319]
        else
            set(gcf,'Position',[738    78   223   206])
        end
        hold on
        SnippetX=(SnippetSize-1)/2+1-...
            (single(fad.channels(CurrentFrame).fits.x(CurrentParticleIndex))-...
            fad.channels(CurrentFrame).fits.x_fit(CurrentParticleIndex));
        SnippetY=(SnippetSize-1)/2+1-...
            (single(fad.channels(CurrentFrame).fits.y(CurrentParticleIndex))-...
            fad.channels(CurrentFrame).fits.y_fit(CurrentParticleIndex));
        PlotHandle=ellipse(fad.channels(CurrentFrame).fits.r_max(CurrentParticleIndex),...
            fad.channels(CurrentFrame).fits.r_min(CurrentParticleIndex),...
            fad.channels(CurrentFrame).fits.theta(CurrentParticleIndex),SnippetX,SnippetY,'g');
        set(PlotHandle,'LineWidth',2.5)
        hold off
    else
        imshow(zeros(SnippetSize))
    end
    
    
    
    
    figure(ZProfileFig)
    if (~isempty(xTrace))
        ZProfile=fad.channels(CurrentFrame).fits.shadowsDog{CurrentParticleIndex};
        ZProfileRaw=fad.channels(CurrentFrame).fits.shadowsRaw{CurrentParticleIndex};

        [Dummy,MaxZ]=max(ZProfile);
        [Dummy,MaxZRaw]=max(ZProfileRaw);
        
        plot([1:length(ZProfile)]-MaxZ+...
            double(fad.channels(CurrentFrame).fits.z(CurrentParticleIndex)),...
            ZProfile,'.-k');
        
        
%         [Axis,Plot1,Plot2]=plotyy([1:length(ZProfile)]-MaxZ+...
%             double(fad.channels(CurrentFrame).fits.z(CurrentParticleIndex)),...
%             ZProfile,...
%             [1:length(ZProfile)]-MaxZ+...
%             double(fad.channels(CurrentFrame).fits.z(CurrentParticleIndex)),...
%             ZProfileRaw);
        hold on
        plot(CurrentZ,ZProfile(MaxZ),'ob')
        hold off
        
        
        set(gca,'XTick',[1:length(ZProfile)]-MaxZ+...
            double(fad.channels(CurrentFrame).fits.z(CurrentParticleIndex)))
        title('Z profile')
        if strcmp('albert-pc',name(1:end-1))
            set(gcf,'Position',[857   518   308   215])
            % fullscreen [1547          44         363         322]
        else
            set(gcf,'Position',[1029          56         216         234])
        end
    end
       
       
    
    figure(TraceFig)
    
    %Only update the trace information if we have switched particles
    if (CurrentParticle~=PreviousParticle)|~exist('Amp')
        PreviousParticle=CurrentParticle;
        [Frames,Amp]=PlotParticleTrace(CurrentParticle,Particles,fad,FISHPath,FilePrefix);
    end
    plot(Frames(Particles(CurrentParticle).FrameApproved),Amp(Particles(CurrentParticle).FrameApproved),'.-k')
    hold on
    plot(Frames(~Particles(CurrentParticle).FrameApproved),Amp(~Particles(CurrentParticle).FrameApproved),'.r')
    plot(Frames(Frames==CurrentFrame),Amp(Frames==CurrentFrame),'ob')
    hold off
    try
        xlim([min(Frames)-1,max(Frames)+1]);
    catch
    end
    xlabel('Frame')
    ylabel('Intensity')
    if strcmp('albert-pc',name(1:end-1))
        set(gcf,'Position',[716    88   439   321])
        % fullscreen [1135         471         776         515]
    else
        set(gcf,'Position',[703   386   564   327])
    end
    
    FigureTitle=['Particle: ',num2str(CurrentParticle),'/',num2str(length(Particles)),...
        ', Frame: ',num2str(CurrentFrame),'/',num2str(TotalFrames),...
        ', Z: ',num2str(CurrentZ),'/',num2str(ZSlices)];
    
    if HideApprovedFlag==1
        FigureTitle=[FigureTitle,', Showing non-flagged particles'];
    elseif HideApprovedFlag==2
        FigureTitle=[FigureTitle,', Showing disapproved particles'];
    end
    title(FigureTitle)
    
    

    
    
    
    figure(Overlay)
    ct=waitforbuttonpress;
    if ~AutoFlag
        cc=get(Overlay,'currentcharacter');
    else
        cc = 13;
    end
    % ES 2014-01-07
    cm=get(gca,'CurrentPoint');
    
    if (cc=='.')&(CurrentFrame<length(fad.channels))
        CurrentFrame=CurrentFrame+1;
        ManualZFlag=0;
        %DisplayRange=[];
    elseif (cc==',')&(CurrentFrame>1)
        CurrentFrame=CurrentFrame-1;
        ManualZFlag=0;
        %DisplayRange=[];
    elseif (cc=='>')&(CurrentFrame+5<length(fad.channels))
        CurrentFrame=CurrentFrame+5;
        ManualZFlag=0;
        %DisplayRange=[];
    elseif (cc=='<')&(CurrentFrame-5>1)
        CurrentFrame=CurrentFrame-5;
        ManualZFlag=0;
        %DisplayRange=[];
    elseif (cc=='a')&(CurrentZ<ZSlices)
        CurrentZ=CurrentZ+1;
        ManualZFlag=1;
    elseif (cc=='z')&(CurrentZ>1)
        CurrentZ=CurrentZ-1;
        ManualZFlag=1;
    elseif cc=='j'
        try
            iJump=input('Frame to jump to: ');
        catch
            iJump=CurrentFrame;
        end
        if (floor(iJump)>0)&(iJump<length(fad.channels))
            CurrentFrame=iJump;
            ManualZFlag=0;
        end
        DisplayRange=[];
    elseif cc=='k'
        try
            ParticleJump=input('Particle to jump to: ');
        catch
            ParticleJump=CurrentParticle;
        end
        if (floor(ParticleJump)>0)&(ParticleJump<=length(Particles))
            CurrentParticle=ParticleJump;
            CurrentFrame=Particles(CurrentParticle).Frame(1);
            ManualZFlag=0;
        end
        
        DisplayRange=[];
    elseif cc=='g'      %Increase histone channel contrast
        if isempty(DisplayRange)
            DisplayRange=[min(min(ImageHis)),max(max(ImageHis))/1.5];
        else
            DisplayRange=[DisplayRange(1),DisplayRange(2)/1.5];
        end
        
    elseif cc=='b'      %Decrease histone channel contrast
        DisplayRange=[min(min(ImageHis)),max(max(ImageHis))*1.5];

        
    elseif cc=='r'
        %Order particles by the earliest frame they appear at. This makes the
        %tracking a lot easier!
        clear FirstFrame
        for i=1:length(Particles)
            FirstFrame(i)=Particles(i).Frame(1);
        end
        [Dummy,Permutations]=sort(FirstFrame);
        Particles=Particles(Permutations);
        
    elseif cc=='f'
        Answer=input('Are you sure you want to redo the tracking?  (y/n) ','s');
        Answer=lower(Answer);
        if Answer=='y'
            %We need to save the data
            save([DataFolder,filesep,'FrameInfo.mat'],'FrameInfo')
            if UseHistoneOverlay
                save([DataFolder,filesep,'Particles.mat'],'Particles','fad','fad2','Threshold1','Threshold2')
                save([DropboxFolder,filesep,FilePrefix(1:end-1),filesep,FilePrefix(1:end-1),'_lin.mat'],'schnitzcells')
            else
                save([DataFolder,filesep,'Particles.mat'],'Particles','fad','fad2','Threshold1','Threshold2')            
            end
            
            
           [Particles,schnitzcells,fad,fad2]=TrackmRNADynamics(FilePrefix(1:end-1),...
               Threshold1,Threshold2); 
           %Check the FrameApproved field
            for i=1:length(Particles)
                if isempty(Particles(i).FrameApproved)
                    Particles(i).FrameApproved=logical(ones(size(Particles(i).Frame)));
                end
            end
        end
    elseif cc=='c'
        PreviousParticle=0;
        if ~sum(Particles(CurrentParticle).Frame==CurrentFrame)
            ConnectPosition=ginput(1);
            
            
            if ~isempty(ConnectPosition)
                [ParticleOutput,IndexOutput]=FindClickedParticle(ConnectPosition,CurrentFrame,fad,Particles);
               %Check that the clicked particle doesn't exist in a previous
                %frame, that there is no overlap of frames. If it does
                %exist in a previous frame we will have to disconnect it.
                if sum(Particles(ParticleOutput).Frame<CurrentFrame)
                    %Disconnect the clicked particle
                    Particles=SeparateParticleTraces(ParticleOutput,CurrentFrame,Particles);
                    ParticleOutput=ParticleOutput+1;
                end
                    
                Particles=JoinParticleTraces(CurrentParticle,ParticleOutput,Particles);
                %Do this in case the clicked particle comes before the current
                %particle in the structure
                if ParticleOutput<CurrentParticle
                    CurrentParticle=ParticleOutput;
                end
                %Sort the frames within the particle. This is useful if we
                %connected to a particle that came before.
                [SortedFrame,Permutations]=sort(Particles(CurrentParticle).Frame);
                Particles(CurrentParticle).Frame=Particles(CurrentParticle).Frame(Permutations);
                Particles(CurrentParticle).Index=Particles(CurrentParticle).Index(Permutations);
                Particles(CurrentParticle).FrameApproved=Particles(CurrentParticle).FrameApproved(Permutations);                

            end
            
        else
            ConnectPosition=ginput(1);
            [ParticleOutput,IndexOutput]=FindClickedParticle(ConnectPosition,CurrentFrame,fad,Particles);
            
            %If it's an independent particle swap it with the frame in the
            %current particle
            if (length(Particles(ParticleOutput).Frame)==1)&...
                    (sum(Particles(ParticleOutput).Frame==CurrentFrame)==1)
                
                ParticleTemp=Particles(ParticleOutput);
                
                %Copy the particle out
                Particles(ParticleOutput).Index=...
                    Particles(CurrentParticle).Index(Particles(CurrentParticle).Frame==CurrentFrame);

                %Copy the new particle in
                Particles(CurrentParticle).Index(Particles(CurrentParticle).Frame==CurrentFrame)=...
                    ParticleTemp.Index;
                Particles(CurrentParticle).FrameApproved(Particles(CurrentParticle).Frame==CurrentFrame)=1;
            else
                display('Cannnot connect to two particles!')
            end
            
            
        end
    elseif cc=='p' %Identify a particle. It will also tell you the particle associated with
                   %  the clicked nucleus.
        ConnectPosition=ginput(1);
        if ~isempty(ConnectPosition)
            %Find the closest particle
            [ParticleOutput,IndexOutput]=FindClickedParticle(ConnectPosition,CurrentFrame,fad,Particles);
            display(['Clicked particle: ',num2str(ParticleOutput)]);
            
            if UseHistoneOverlay
                %Find the closest nucleus
                NewNuclei=ConnectPosition;

                %Find which schnitz this corresponds to
                SchnitzSuspect=[];
                xPosSuspect=[];
                yPosSuspect=[];
                for j=1:length(schnitzcells)
                    if sum(schnitzcells(j).frames-1==CurrentFrame)
                        SchnitzSuspect=[SchnitzSuspect,j];
                        xPosSuspect=[xPosSuspect,...
                            schnitzcells(j).cenx(find((schnitzcells(j).frames)==CurrentFrame))];
                        yPosSuspect=[yPosSuspect,...
                            schnitzcells(j).ceny(find((schnitzcells(j).frames)==CurrentFrame))];
                    end
                end

                %Find the closest one to the point where we clicked
                Distance=sqrt((NewNuclei(1)-xPosSuspect).^2+(NewNuclei(2)-yPosSuspect).^2);
                [MinValue,ClosestNucleusIndex]=min(Distance);

                ClickedSchnitz=SchnitzSuspect(ClosestNucleusIndex);

                %Now, find its associated particle
                for i=1:length(Particles)
                    if ~isempty(Particles(i).Nucleus)
                        AssignedNuclei(i)=Particles(i).Nucleus;
                    else
                        AssignedNuclei(i)=nan;
                    end
                end
                AssociatedParticle=find(AssignedNuclei==ClickedSchnitz);

                if isempty(AssociatedParticle)
                    display(['Nucleus ',num2str(ClickedSchnitz),' does not have an associated particle'])
                else
                    display(['Particle ',num2str(AssociatedParticle),' is associate with nucleus ',...
                        num2str(ClickedSchnitz)])
                end
            end
        end
    elseif cc=='1'  %Give the nucleus number in the schnitzcell segmentation structure. This
                    %only  works for troubleshooting and you need to be online and on the
                    %Princeton network/VPN for now.
        if exist('Z:\FISHDrosophila\Analysis\schnitzcells\Analysis')
            ClickedPosition=round(ginput(1));
            if ~isempty(ClickedPosition)
                Dashes=strfind(FilePrefix,'-');
                FileDate=FilePrefix(1:Dashes(3)-1);
                load(['Z:\FISHDrosophila\Analysis\schnitzcells\Analysis\',...
                    FileDate,filesep,FilePrefix(1:end-1),'\segmentation\',...
                    FilePrefix(1:end-1),'seg',iIndex(CurrentFrame,3),'.mat'],...
                    '-mat','LN')
                %See which nucleus it corresponds to
                ClickedNucleus=LN(ClickedPosition(2),ClickedPosition(1));
                display(['Clicked nucleus: ',num2str(ClickedNucleus)])
                
                %Now check if there's any schnitzcell associated with this
                SchnitzFound=0;
                for i=1:length(schnitzcells)
                    if sum((schnitzcells(i).frames)==CurrentFrame)
                        IndexToSearch=find((schnitzcells(i).frames)==CurrentFrame);
                        try 
                            if schnitzcells(i).cellno(IndexToSearch)==ClickedNucleus
                                SchnitzFound=i;
                            end
                        catch
                            display(['WARNING: Consistency error in schnitz ',num2str(i)])
                        end
                    end
                end
                
                if SchnitzFound==0
                    display('No schnitz associated with it! Creating it...')
                    
                    Nschnitz=length(schnitzcells);

                    schnitzcells(Nschnitz+1).P=0;
                    schnitzcells(Nschnitz+1).D=0;
                    schnitzcells(Nschnitz+1).E=0;
                    schnitzcells(Nschnitz+1).frames=CurrentFrame;
                    schnitzcells(Nschnitz+1).cenx=Ellipses{CurrentFrame}(ClickedNucleus,1)+1;
                    schnitzcells(Nschnitz+1).ceny=Ellipses{CurrentFrame}(ClickedNucleus,2)+1;
                    %schnitzcells(Nschnitz+1).ang=Ellipses{CurrentFrame}(ClickedNucleus,5);
                    schnitzcells(Nschnitz+1).len=max([Ellipses{CurrentFrame}(ClickedNucleus,3),...
                        Ellipses{CurrentFrame}(ClickedNucleus,4)]);
                    schnitzcells(Nschnitz+1).cellno=ClickedNucleus;
                else
                    display(['It is associated with schnitz ',num2str(SchnitzFound)])
                end
            end
           
        else
            display('This feature only works on Hernan''s and if you''re connected to the VPN')
        end
                    
    elseif cc=='u'
         [x2,y2]=fad2xyzFit(CurrentFrame,fad2, 'addMargin'); 
         if ~isempty(x2)
            fad2Position=ginput(1);
            [fad,fad2,Particles]=...
                Integratefad2Particle(fad,fad2,fad2Position,Particles,CurrentFrame);
         end
    elseif cc=='i'
        PreviousParticle=0;
         [x2,y2]=fad2xyzFit(CurrentFrame,fad2, 'addMargin'); 
         if ~isempty(x2)
            fad2Position=ginput(1);
            if (~isempty(fad2Position))
                [fad,fad2,Particles]=...
                    Integratefad2Particle(fad,fad2,fad2Position,Particles,CurrentFrame);
            end
         end

         
         if (~sum(Particles(CurrentParticle).Frame==CurrentFrame))&(~isempty(fad2Position))
            ConnectPosition=fad2Position;
            [ParticleOutput,IndexOutput]=FindClickedParticle(ConnectPosition,CurrentFrame,fad,Particles);
            
            %Check that the clicked particle doesn't exist in a previous
            %frame, that there is no overlap of frams.  Maybe I can have
            %those in a different color.

           
            if sum(Particles(ParticleOutput).Frame<CurrentFrame)
                display(['Target particle (',num2str(ParticleOutput),') is already in a previous frame!']);
            else
                Particles=JoinParticleTraces(CurrentParticle,ParticleOutput,Particles);
                %Do this in case the clicked particle comes before the current
                %particle in the structure
                if ParticleOutput<CurrentParticle
                    CurrentParticle=ParticleOutput;
                end
                %Sort the frames within the particle. This is useful if we
                %connected to a particle that came before.
                [SortedFrame,Permutations]=sort(Particles(CurrentParticle).Frame);
                Particles(CurrentParticle).Frame=Particles(CurrentParticle).Frame(Permutations);
                Particles(CurrentParticle).Index=Particles(CurrentParticle).Index(Permutations);
                Particles(CurrentParticle).FrameApproved=Particles(CurrentParticle).FrameApproved(Permutations);
                
                if UseHistoneOverlay
                    %Check for consistency within schnitzcell 
                    [Particles,schnitzcells]=CheckSchnitzLineage(Particles,CurrentParticle,schnitzcells,CurrentFrame,...
                        Overlay);
                end
                
                
                
                
            end
            
        else
            display('Cannnot connect to two particles!')
        end
         
         
         
         
     elseif cc=='d'  %d Separate traces. The separated particle won't have a nucleus assigned!
         PreviousParticle=0;
        %Check that the particle does actually exist in this frame
        if ~(Particles(CurrentParticle).Frame(1)==CurrentFrame)
            if sum(Particles(CurrentParticle).Frame==CurrentFrame)
                Particles=SeparateParticleTraces(CurrentParticle,CurrentFrame,Particles);
            end
        elseif length(Particles(CurrentParticle).Frame)==1
            Particles(CurrentParticle).Nucleus=[];
        else
            display('Cannot divide a trace at the first time point')
        end
    elseif cc=='q'      %Approve a trace
        if Particles(CurrentParticle).Approved==1
            Particles(CurrentParticle).Approved=2;
        elseif Particles(CurrentParticle).Approved==0
            Particles(CurrentParticle).Approved=1;
        elseif Particles(CurrentParticle).Approved==2
            Particles(CurrentParticle).Approved=0;
        end
    elseif cc=='w'      %Disapproove a trace
        if Particles(CurrentParticle).Approved==-1
            Particles(CurrentParticle).Approved=0;
        else
            Particles(CurrentParticle).Approved=-1;
        end    
        
    elseif cc=='s'
        save([DataFolder,filesep,'FrameInfo.mat'],'FrameInfo')
        if UseHistoneOverlay
            save([DataFolder,filesep,'Particles.mat'],'Particles','fad','fad2','Threshold1','Threshold2')
            save([DropboxFolder,filesep,FilePrefix(1:end-1),filesep,FilePrefix(1:end-1),'_lin.mat'],'schnitzcells')
        else
            save([DataFolder,filesep,'Particles.mat'],'Particles','fad','fad2','Threshold1','Threshold2')            
        end
        display('Particles saved.')
    elseif cc=='t'
        ShowThreshold2=~ShowThreshold2;
    elseif (cc=='y')&(~UseHistoneOverlay)
            FrameInfo=DetermineNC(fad,Particles,FrameInfo);
    elseif cc=='h'
        if HideApprovedFlag==0
            HideApprovedFlag=1;         %Show only non-approved traces
        elseif HideApprovedFlag==1
            HideApprovedFlag=2;         %Show only yellow and red traces
        elseif HideApprovedFlag==2
            HideApprovedFlag=0;
        end
            
        %HideApprovedFlag=~HideApprovedFlag;
    elseif cc=='o'
        ZoomMode=~ZoomMode;        
    elseif (cc=='m')&(CurrentParticle<length(Particles))
        
        NextParticle=CurrentParticle+1;
        
        if NextParticle>length(Particles)
            NextParticle=length(Particles);
        end
        
        
        %Mode 1 - skip approved or flagged traces
        while (HideApprovedFlag)==1&(NextParticle<length(Particles))&...
                ((Particles(NextParticle).Approved==1)|(Particles(NextParticle).Approved==-1)|...
                (Particles(NextParticle).Approved==2))
                NextParticle=NextParticle+1;
        end
        
        %Mode 2 - skip approved traces
        while ((HideApprovedFlag)==2)&(NextParticle<length(Particles))&...
                ((Particles(NextParticle).Approved==1)|(Particles(NextParticle).Approved==2))
                NextParticle=NextParticle+1;
        end
        
        
        CurrentParticle=NextParticle;
%         if ~isempty(Particles(CurrentParticle).Nucleus)
%             CurrentFrame=schnitzcells(Particles(CurrentParticle).Nucleus).frames(1)-1;
%         else
%             CurrentFrame=Particles(CurrentParticle).Frame(1);
%         end
        CurrentFrame=Particles(CurrentParticle).Frame(1);

        ParticleToFollow=[];
        DisplayRange=[];
        
        
        display('Missing frames:')
        Particles(CurrentParticle).Frame(find(diff(Particles(CurrentParticle).Frame)>1))
        
    elseif (cc=='n')&(CurrentParticle>1)
        Approved=(find([Particles.Approved]));
        %NotApproved=(find(~[Particles.Approved]));
        
        NextParticle=CurrentParticle-1;
        
        
        
        %Mode 1 - show non-flagged traces
        while (HideApprovedFlag)==1&(NextParticle>1)&...
                ((Particles(NextParticle).Approved==1)|(Particles(NextParticle).Approved==-1)|...
                (Particles(NextParticle).Approved==2))
            NextParticle=NextParticle-1;
            if NextParticle<1
                NextParticle=1;
            end
        end
        
        
        %Mode 2 - show dissapproved traces
        while ((HideApprovedFlag)==2)&(NextParticle>1)&...
                ((Particles(NextParticle).Approved==1)|(Particles(NextParticle).Approved==2))
            NextParticle=NextParticle-1;
            if NextParticle<1
                NextParticle=1;
            end
        end
        
        
        if NextParticle<1
            NextParticle=CurrentParticle;
        end
        
        CurrentParticle=NextParticle;
        
%         if ~isempty(Particles(CurrentParticle).Nucleus)
%             CurrentFrame=schnitzcells(Particles(CurrentParticle).Nucleus).frames(1)-1;
%         else
%             CurrentFrame=Particles(CurrentParticle).Frame(1);
%         end
        CurrentFrame=Particles(CurrentParticle).Frame(1);

        
        ParticleToFollow=[];
        
        DisplayRange=[];
        
    elseif cc=='e'
        Particles(CurrentParticle).FrameApproved(Particles(CurrentParticle).Frame==CurrentFrame)=...
            ~Particles(CurrentParticle).FrameApproved(Particles(CurrentParticle).Frame==CurrentFrame);
    
    
    
    %Schnitzcells specific
    
    elseif cc=='l' %Split a nucleus and select one or two daughter cells or stop the lineage
        PreviousParticle=0;
        
        display('Select one/two daughter nuclei and press ENTER or just press ENTER to terminate lineage')
        NewNuclei=ginput(2);
        
        if isempty(NewNuclei)
            warning('Write the disconnect schnitz part')
        else
            [NClicks,Dummy]=size(NewNuclei);
        
            ClickedSchnitz=[];
            %Find the nuclei/schnitz that we clicked on
            for i=1:NClicks
                %Find which schnitz this corresponds to
                SchnitzSuspect=[];
                xPosSuspect=[];
                yPosSuspect=[];
                for j=1:length(schnitzcells)
                    if sum(schnitzcells(j).frames-1==CurrentFrame)
                        SchnitzSuspect=[SchnitzSuspect,j];
                        if (~isempty(schnitzcells(j).cenx(find((schnitzcells(j).frames)==CurrentFrame))))&...
                                (~isempty(schnitzcells(j).ceny(find((schnitzcells(j).frames)==CurrentFrame))))
                            xPosSuspect=[xPosSuspect,...
                                schnitzcells(j).cenx(find((schnitzcells(j).frames)==CurrentFrame))];
                            yPosSuspect=[yPosSuspect,...
                                schnitzcells(j).ceny(find((schnitzcells(j).frames)==CurrentFrame))];
                        else
                            xPosSuspect=[xPosSuspect,inf];
                            yPosSuspect=[yPosSuspect,inf];
                        end
                    end
                end

                %Find the closest one to the point where we clicked
                Distance=sqrt((NewNuclei(i,1)-xPosSuspect).^2+(NewNuclei(i,2)-yPosSuspect).^2);
                [MinValue,ClosestNucleusIndex]=min(Distance);

                ClickedSchnitz(i)=SchnitzSuspect(ClosestNucleusIndex);
            end
            
            %Now look at the different cases
            
            
%       Click on one nucleus + ENTER: Continue the schnitz with that nucleus.
%       Click on two nuclei: Split the current nucleus into two daughter
%       nuclei.
%       Click on the same nucleus twice: Split the current nucleus, but
%       with only one daughter nucleus.
            
            
            if length(ClickedSchnitz)==1
                if Particles(CurrentParticle).Nucleus==ClickedSchnitz %Split the lineage
                    [Particles,schnitzcells]=SplitSchnitz(Particles,schnitzcells,...
                            CurrentFrame,...
                            CurrentParticle);
                else
                    try
                        [Particles,schnitzcells]=JoinSchnitz(Particles,schnitzcells,Particles(CurrentParticle).Nucleus,...  
                            ClickedSchnitz,CurrentFrame);
                    catch
                        display('Error in JoinSchnitz')
                    end
                end
            elseif length(ClickedSchnitz)==2
                if ClickedSchnitz(1)~=ClickedSchnitz(2)
                    [Particles,schnitzcells]=SplitSchnitzDaughters(Particles,schnitzcells,...
                            CurrentFrame,...
                            Particles(CurrentParticle).Nucleus,ClickedSchnitz(1),ClickedSchnitz(2));
                else
                    [Particles,schnitzcells]=SplitSchnitzDaughters(Particles,schnitzcells,...
                        CurrentFrame,...
                        Particles(CurrentParticle).Nucleus,ClickedSchnitz(1),0);
                end
            else
                error('Too many cells selected')
            end
        end
        
        
    elseif cc=='2' %2 set parent of current nucleus
        display('Select the mother nucleus or press enter to delete mother information')
        NewNuclei=ginput(1);
        
        if isempty(NewNuclei)
            error('Write the disconnect schnitz part')
        else
            [NClicks,Dummy]=size(NewNuclei);
        
            ClickedSchnitz=[];
            %Find the nuclei/schnitz that we clicked on
            for i=1:NClicks
                %Find which schnitz this corresponds to
                SchnitzSuspect=[];
                xPosSuspect=[];
                yPosSuspect=[];
                for j=1:length(schnitzcells)
                    if sum(schnitzcells(j).frames==CurrentFrame)
                        SchnitzSuspect=[SchnitzSuspect,j];
                        if (~isempty(schnitzcells(j).cenx(find((schnitzcells(j).frames)==CurrentFrame))))&...
                                (~isempty(schnitzcells(j).ceny(find((schnitzcells(j).frames)==CurrentFrame))))
                            xPosSuspect=[xPosSuspect,...
                                schnitzcells(j).cenx(find((schnitzcells(j).frames)==CurrentFrame))];
                            yPosSuspect=[yPosSuspect,...
                                schnitzcells(j).ceny(find((schnitzcells(j).frames)==CurrentFrame))];
                        else
                            xPosSuspect=[xPosSuspect,inf];
                            yPosSuspect=[yPosSuspect,inf];
                        end
                    end
                end

                %Find the closest one to the point where we clicked
                Distance=sqrt((NewNuclei(i,1)-xPosSuspect).^2+(NewNuclei(i,2)-yPosSuspect).^2);
                [MinValue,ClosestNucleusIndex]=min(Distance);

                ClickedSchnitz(i)=SchnitzSuspect(ClosestNucleusIndex);
            end
            
            %Now look at the different cases. Note that I don't have a good
            %way to fix the parent nucleus itself. This might be a bad idea
            %after all
            
       
            
            if length(ClickedSchnitz)==1
                schnitzcells(Particles(CurrentParticle).Nucleus).P=ClickedSchnitz;

            elseif length(ClickedSchnitz)==0
                schnitzcells(Particles(CurrentParticle).Nucleus).P=0;
            end                
        end    
        
    elseif cc=='0'      %Debugging mode
        keyboard;
    elseif cc=='9'  %check for nuclear tracking consistencies. This is useful while we're
                    %  getting the code to work well.
        warning('This feature has been discontinued for now. Talk to HG.')
        [schnitzcells,Particles]=CheckSchnitzConsistency(schnitzcells,Particles,Ellipses);
        PreviousParticle=0;

    end
        
end


save([DataFolder,filesep,'FrameInfo.mat'],'FrameInfo')

if UseHistoneOverlay
    save([DataFolder,filesep,'Particles.mat'],'Particles','fad','fad2','Threshold1','Threshold2')
    save([DropboxFolder,filesep,FilePrefix(1:end-1),filesep,FilePrefix(1:end-1),'_lin.mat'],'schnitzcells')
else
    save([DataFolder,filesep,'Particles.mat'],'Particles','fad','fad2','Threshold1','Threshold2')            
end
display('Particles saved.')






%% Extra stuff that is useful in debug mode

%Reset approve status of all approved particles in a certain nc
% nc=13;
% 
% %Determine the start and end frame of the nc
% if nc==14
%     display('Do this')
% else
%     eval(['ncStart=nc',num2str(nc),';']);
%     eval(['ncEnd=nc',num2str(nc+1),';']);
% end
% 
% for i=1:length(Particles)
%     if Particles(i).Approved==1
%     
%         if (min(Particles(i).Frame(Particles(i).FrameApproved))>=ncStart)&...
%                 (min(Particles(i).Frame(Particles(i).FrameApproved))<ncEnd)
%             Particles(i).Approved=0;
%         end
%     end
% end


%This is if the reference from schnitz to nucleus is screwed up. It
%manifests itself in that cenx and ceny of the schnitz will not coincide
%with the actual ellipse. In this case we delete the schnitz (emptying it)
%and the nuclear reference and use '1' to recreate the schnitz.

% ParticlesToDecouple=[314:316,318,320:325,328,329,331:333,335,337,339,341,342,344:346]
% 
% for ParticleToDecouple=ParticlesToDecouple
% 
%     SchnitzToDecouple=Particles(ParticleToDecouple).Nucleus;
% 
%     %Delete the schnitz by emptying it
%     schnitzcells(SchnitzToDecouple).P=[];   
%     schnitzcells(SchnitzToDecouple).E=[];
%     schnitzcells(SchnitzToDecouple).D=[];
%     schnitzcells(SchnitzToDecouple).frames=[];
%     schnitzcells(SchnitzToDecouple).cenx=[];
%     schnitzcells(SchnitzToDecouple).ceny=[];
%     schnitzcells(SchnitzToDecouple).len=[];
%     schnitzcells(SchnitzToDecouple).cellno=[];
% 
%     %Decouple particle and schnitz
%     Particles(ParticleToDecouple).Nucleus=[];
% end
