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
%8 Change channels
%t Show/hide particles from the second threshold
%s Save the current Particles structure
%x Save and exit
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

[SourcePath,FISHPath,DefaultDropboxFolder,MS2CodePath,PreProcPath]=...
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
    [SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
        DetermineLocalFolders(varargin{1});
    DataFolder=[DropboxFolder,filesep,varargin{1}];
end
    
%Flag to sort or not particles according to their starting frame
NoSort=0;
%Flag to just save the data. This is good for CompileAll
ForCompileAll=0;
%Flag to plot only ellipses for current particle & save time
SpeedMode = 0;

if length(varargin)>1
    for i=2:length(varargin)
        if strcmp(varargin{i},'NoSort')
            NoSort=1;
        elseif strcmp(varargin{i},'ForCompileAll')
            ForCompileAll=1;
        elseif strcmpi(varargin{i}, 'speedmode')
            SpeedMode = 1;
        end
    end
end

%%


FilePrefix=[DataFolder(length(DropboxFolder)+2:end),'_'];

%Now get the actual folders
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders(FilePrefix(1:end-1));

load([DataFolder,filesep,'Particles.mat'])

%Check that FrameInfo exists
if exist([DataFolder,filesep,'FrameInfo.mat'])
    load([DataFolder,filesep,'FrameInfo.mat'])
else
    warning('No FrameInfo.mat found. Trying to continue')
    %Adding frame information
    DHis=dir([PreProcPath,filesep,FilePrefix(1:end-1),filesep,'*His*.tif']);
    FrameInfo(length(DHis)).nc=[];
    %Adding information

    Dz=dir([PreProcPath,filesep,FilePrefix(1:end-1),filesep,FilePrefix(1:end-1),'*001*.tif']);
    NumberSlices=length(Dz)-1;
    
    for i=1:length(FrameInfo)
        FrameInfo(i).NumberSlices=NumberSlices;
    end
    
end


%Some parameters:
SnippetSize=13;     %Size of the snippets generated by Michael's code


%Create the particle array. This is done so that we can support multiple
%channels. Also figure out the number of channels
if iscell(Particles)
    NChannels=length(Particles);
else
    Particles={Particles};
    NChannels=1;
end


%Add FramesApproved where necessary
if ~isfield(Particles{1},'FrameApproved')
    for NCh=1:NChannels
        for i=1:length(Particles{NCh})
            Particles{NCh}(i).FrameApproved=logical(ones(size(Particles{NCh}(i).Frame)));
        end
    end
else
    for NCh=1:NChannels
        for i=1:length(Particles{NCh})
            if isempty(Particles{NCh}(i).FrameApproved)
                Particles{NCh}(i).FrameApproved=logical(ones(size(Particles{NCh}(i).Frame)));
            end
        end
    end
end



%Check if we have the histone channel and we have done the nuclear
%segmentation.
if exist([PreProcPath,filesep,FilePrefix(1:end-1),filesep,...
        FilePrefix(1:end-1),'-His_',iIndex(1,3),'.tif'])|...
        exist([PreProcPath,filesep,FilePrefix(1:end-1),filesep,...
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
Channel1Column=find(strcmp(XLSRaw(1,:),'Channel1'));
Channel2Column=find(strcmp(XLSRaw(1,:),'Channel2'));
nc9Column=find(strcmp(XLSRaw(1,:),'nc9'));
nc10Column=find(strcmp(XLSRaw(1,:),'nc10'));
nc11Column=find(strcmp(XLSRaw(1,:),'nc11'));
nc12Column=find(strcmp(XLSRaw(1,:),'nc12'));
nc13Column=find(strcmp(XLSRaw(1,:),'nc13'));
nc14Column=find(strcmp(XLSRaw(1,:),'nc14'));
CFColumn=find(strcmp(XLSRaw(1,:),'CF'));

DataFolderColumn=find(strcmp(XLSRaw(1,:),'DataFolder'));
Dashes=findstr(FilePrefix,'-');
PrefixRow=find(strcmp(XLSRaw(:,DataFolderColumn),[FilePrefix(1:Dashes(3)-1),'\',FilePrefix(Dashes(3)+1:end-1)]));
if isempty(PrefixRow)
    PrefixRow=find(strcmp(XLSRaw(:,DataFolderColumn),[FilePrefix(1:Dashes(3)-1),'/',FilePrefix(Dashes(3)+1:end-1)]));
    if isempty(PrefixRow)
        error('Could not find data set in MovieDatabase.XLSX. Check if it is defined there.')
    end
end


ExperimentType=XLSRaw{PrefixRow,ExperimentTypeColumn};
ExperimentAxis=XLSRaw{PrefixRow,ExperimentAxisColumn};
Channel1=XLSRaw(PrefixRow,Channel1Column);
Channel2=XLSRaw(PrefixRow,Channel2Column);


%Find the corresponding entry in the XLS file
if (~isempty(findstr(FilePrefix,'Bcd')))&(isempty(findstr(FilePrefix,'BcdE1')))&...
        (isempty(findstr(FilePrefix(1:end-1),'NoBcd')))&...
        (isempty(findstr(FilePrefix(1:end-1),'Bcd1x')))
    warning('This step in CheckParticleTracking will most likely have to be modified to work')
    XLSEntry=find(strcmp(XLSRaw(:,DataFolderColumn),...
        [Date,'\BcdGFP-HisRFP']));
else
    XLSEntry=find(strcmp(XLSRaw(:,DataFolderColumn),...
        [FilePrefix(1:Dashes(3)-1),'\',FilePrefix(Dashes(3)+1:end-1)]));
    if isempty(XLSEntry)
        XLSEntry=find(strcmp(XLSRaw(:,DataFolderColumn),...
            [FilePrefix(1:Dashes(3)-1),'/',FilePrefix(Dashes(3)+1:end-1)]));
        if isempty(XLSEntry)
            error('Could not find data set in MovieDatabase.XLSX. Check if it is defined there.')
        end
    end
end


if strcmp(XLSRaw(XLSEntry,Channel2Column),'His-RFP')|strcmp(XLSRaw(XLSEntry,Channel1Column),'His-RFP')
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
    warning('Warning: no histone channel may result in strange behavior.');
    
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
end





%Check if we have already determined nc
if (~isfield(FrameInfo,'nc'))&&(~UseHistoneOverlay)
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
    for ChN=1:NChannels
        for i=1:length(Particles{ChN})
            FirstFrame(i)=Particles{ChN}(i).Frame(1);
        end
        [Dummy,Permutations]=sort(FirstFrame);
        Particles{ChN}=Particles{ChN}(Permutations);
        clear FirstFrame
    end
end

TotalFrames=length(fad(NChannels).channels);


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
CurrentChannel=1;
PreviousChannel=CurrentChannel;


CurrentFrame=Particles{1}(1).Frame(1);


DisplayRange=[];

ZoomMode=0;
ZoomRange=50;










%Determine the positions and size of the figures
ScreenSize=get( 0, 'ScreenSize' );
ScreenSize=ScreenSize(3:end);
ScreenRows=ScreenSize(2);
ScreenColumns=ScreenSize(1);
%Define the windows
Overlay=figure;
if UseHistoneOverlay 
    HisOverlayFig=figure;   
end
TraceFig=figure;
SnippetFig=figure;
ZProfileFig=figure;


cc=1;

%See if we just want to save the data
if ForCompileAll
    % Create the approved field if it does not exist
    if ~isfield(Particles{1},'Approved')
        for NCh=1:NChannels
            for i=1:length(Particles{NCh})
                Particles{NCh}(i).Approved=0;
            end
        end
    end    
    
    cc='x';
end


while (cc~='x') 
    EllipseHandle=[];
    EllipseHandleYellow=[];
    EllipseHandleBlue=[];
    EllipseHandleWhite=[];
    EllipseHandleGreen=[];

    
    %Get the coordinates taking the margins into account
    [x,y]=fad2xyzFit(CurrentFrame,fad(CurrentChannel), 'addMargin'); 
    
    %If the approved field does not exist create it
    if ~isfield(Particles{CurrentChannel},'Approved')
        for i=1:length(Particles{CurrentChannel})
            Particles{CurrentChannel}(i).Approved=0;
        end
    end
    ApprovedParticles=[Particles{CurrentChannel}.Approved];
    
    %Pull out the right particle if it exists in this frame
    CurrentParticleIndex=Particles{CurrentChannel}(CurrentParticle).Index(Particles{CurrentChannel}(CurrentParticle).Frame==CurrentFrame);
    
    %This is the position of the current particle
    xTrace=x(CurrentParticleIndex);
    yTrace=y(CurrentParticleIndex);
    
    %These are the positions of all the approved and disapproved particles
    %Find the particles in this frame
    IndexApprovedParticles=[];
    for i=1:length(Particles{CurrentChannel})
        if sum(Particles{CurrentChannel}(i).Frame==CurrentFrame)&sum(Particles{CurrentChannel}(i).Approved==1)
            IndexApprovedParticles=[IndexApprovedParticles,...
                Particles{CurrentChannel}(i).Index(Particles{CurrentChannel}(i).Frame==CurrentFrame)];
        end
    end
    xApproved=x(IndexApprovedParticles);
    yApproved=y(IndexApprovedParticles);

    IndexDisapprovedParticles=[];
    for i=1:length(Particles{CurrentChannel})
        if sum(Particles{CurrentChannel}(i).Frame==CurrentFrame)&sum(Particles{CurrentChannel}(i).Approved==-1)
            IndexDisapprovedParticles=[IndexDisapprovedParticles,...
                Particles{CurrentChannel}(i).Index(Particles{CurrentChannel}(i).Frame==CurrentFrame)];
        end
    end
    xDisapproved=x(IndexDisapprovedParticles);
    yDisapproved=y(IndexDisapprovedParticles);
    
    
    if (~isempty(xTrace))&(~ManualZFlag)
        CurrentZ=...
            fad(CurrentChannel).channels(CurrentFrame).fits.z(CurrentParticleIndex);
        ManualZFlag=0;
    end
    
    
    if (NChannels==1)&(~strcmp(lower(ExperimentType),'inputoutput'))
        try
            Image=imread([PreProcPath,filesep,FilePrefix(1:end-1),filesep,...
                FilePrefix,iIndex(CurrentFrame,3),'_z',iIndex(CurrentZ,2),'.tif']);
        catch
            display(['Warning: Could not load file: ',...
                FilePrefix,iIndex(CurrentFrame,3),'_z',iIndex(CurrentZ,2),'.tif'])
        end
    elseif (NChannels==1)&(strcmp(lower(ExperimentType),'inputoutput'))
        OutputChannelTemp1=strfind({lower(Channel1{1}),lower(Channel2{1})},'mcp');
        OutputChannelTemp2=strfind({lower(Channel1{1}),lower(Channel2{1})},'pcp');
        OutputChannelTemp1=~cellfun(@isempty,OutputChannelTemp1);
        OutputChannelTemp2=~cellfun(@isempty,OutputChannelTemp2);
        OutputChannel=find(OutputChannelTemp1|OutputChannelTemp2);
        
        Image=imread([PreProcPath,filesep,FilePrefix(1:end-1),filesep,...
                FilePrefix,iIndex(CurrentFrame,3),'_z',iIndex(CurrentZ,2),'_ch',iIndex(OutputChannel,2),'.tif']);
    else
        try
            Image=imread([PreProcPath,filesep,FilePrefix(1:end-1),filesep,...
                FilePrefix,iIndex(CurrentFrame,3),'_z',iIndex(CurrentZ,2),'_ch',iIndex(CurrentChannel,2),'.tif']);
        catch
            display(['Warning: Could not load file: ',...
                FilePrefix,iIndex(CurrentFrame,3),'_z',iIndex(CurrentZ,2),'_ch',iIndex(CurrentChannel,2),'.tif']);
        end
    end

    figure(Overlay)
    imshow(Image,[],'Border','Tight')
    hold on
    %Show all particles in regular mode
    if ~SpeedMode
        plot(x,y,'or')
        plot(xApproved,yApproved,'ob')
        plot(xDisapproved,yDisapproved,'^r')
    end
    %Always show current particle
    plot(xTrace,yTrace,'og')
    hold off

    
    set(gcf,'Name',['Particle: ',num2str(CurrentParticle),'/',num2str(length(Particles{CurrentChannel})),...
        ', Frame: ',num2str(CurrentFrame),'/',num2str(TotalFrames),...
        ', Z: ',num2str(CurrentZ),'/',num2str(ZSlices),' nc: ', num2str(FrameInfo(CurrentFrame).nc),...
        ', Ch: ',num2str(CurrentChannel)])
    
    if UseSchnitz
        %Show all the nuclei in regular mode
        if ~SpeedMode
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
            for i=1:length(Particles{CurrentChannel})
                if Particles{CurrentChannel}(i).Approved==1
                    schnitzIndex=find((schnitzcells(Particles{CurrentChannel}(i).Nucleus).frames)==CurrentFrame);
                    schnitzCellNo=[schnitzCellNo,schnitzcells(Particles{CurrentChannel}(i).Nucleus).cellno(schnitzIndex)];
                end
            end
            
            EllipseHandleBlue=notEllipse(Ellipses{CurrentFrame}(schnitzCellNo,3),...
                Ellipses{CurrentFrame}(schnitzCellNo,4),...
                Ellipses{CurrentFrame}(schnitzCellNo,5),...
                Ellipses{CurrentFrame}(schnitzCellNo,1)+1,...
                Ellipses{CurrentFrame}(schnitzCellNo,2)+1,'b',50);
            
            
            hold off
        end
        
        %Show the corresponding nucleus
        if ~isempty(Particles{CurrentChannel}(CurrentParticle).Nucleus)
            SchnitzIndex=find(schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).frames==CurrentFrame);
            NucleusIndex=schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).cellno(SchnitzIndex);

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
            DaughterE=schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).E;
            DaughterD=schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).D;

            
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
            Mother=schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).P;

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
        [x2,y2]=fad2xyzFit(CurrentFrame,fad2(CurrentChannel), 'addMargin'); 
        hold on
        plot(x2,y2,'sr')
        hold off
    end
    
    if ZoomMode
        %Find the closest frame
        [Dummy,MinIndex]=min((Particles{CurrentChannel}(CurrentParticle).Frame-CurrentFrame).^2);
        if length(MinIndex)>1
            MinIndex=MinIndex(1);
        end
        [xForZoom,yForZoom]=fad2xyzFit(Particles{CurrentChannel}(CurrentParticle).Frame(MinIndex),fad(CurrentChannel), 'addMargin'); 
        xForZoom=xForZoom(Particles{CurrentChannel}(CurrentParticle).Index(MinIndex));
        yForZoom=yForZoom(Particles{CurrentChannel}(CurrentParticle).Index(MinIndex));
       
 
        xlim([xForZoom-ZoomRange,xForZoom+ZoomRange])
        ylim([yForZoom-ZoomRange/2,yForZoom+ZoomRange/2])
    end
    

    if UseHistoneOverlay
        figure(HisOverlayFig)
        
        try
            ImageHis=imread([PreProcPath,filesep,FilePrefix(1:end-1),filesep,...
                FilePrefix(1:end-1),'-His_',iIndex(CurrentFrame,3),'.tif']);
        catch %Had to do this for KITP
            ImageHis=imread([PreProcPath,filesep,FilePrefix(1:end-1),filesep,...
                FilePrefix(1:end-1),'_His_',iIndex(CurrentFrame,3),'.tif']);
        end

        if isempty(DisplayRange)
            HisOverlayImage=cat(3,mat2gray(ImageHis)*2,mat2gray(Image),zeros(size(Image)));
        else
            HisOverlayImage=cat(3,mat2gray(ImageHis,double(DisplayRange))*2,mat2gray(Image),zeros(size(Image)));
        end
        imshow(HisOverlayImage,'Border','Tight')

       
        hold on
        if ~SpeedMode
            plot(x,y,'ow')
            plot(xApproved,yApproved,'ob')
        end
        plot(xTrace,yTrace,'og')
        hold off
        
        if ShowThreshold2
            [x2,y2]=fad2xyzFit(CurrentFrame,fad2(CurrentChannel), 'addMargin'); 
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
            
        end
    
        set(gcf,'Name',['Particle: ',num2str(CurrentParticle),'/',num2str(length(Particles{CurrentChannel})),...
            ', Frame: ',num2str(CurrentFrame),'/',num2str(TotalFrames),...
            ', Z: ',num2str(CurrentZ),'/',num2str(ZSlices),' nc: ', num2str(FrameInfo(CurrentFrame).nc),...
            ' Ch: ',num2str(CurrentChannel)])
        
        if ZoomMode
            %Find the closest frame
            [Dummy,MinIndex]=min((Particles{CurrentChannel}(CurrentParticle).Frame-CurrentFrame).^2);
            if length(MinIndex)>1
                MinIndex=MinIndex(1);
            end
            [xForZoom,yForZoom]=fad2xyzFit(Particles{CurrentChannel}(CurrentParticle).Frame(MinIndex),fad(CurrentChannel), 'addMargin'); 
            xForZoom=xForZoom(Particles{CurrentChannel}(CurrentParticle).Index(MinIndex));
            yForZoom=yForZoom(Particles{CurrentChannel}(CurrentParticle).Index(MinIndex));


            xlim([xForZoom-ZoomRange,xForZoom+ZoomRange])
            ylim([yForZoom-ZoomRange/2,yForZoom+ZoomRange/2])
        end
    end
    
    

    
    
    
    figure(SnippetFig)
    if (~isempty(xTrace))
        
        CurrentSnippet=mat2gray(fad(CurrentChannel).channels(CurrentFrame).fits.snippets(:,:,CurrentParticleIndex));
        IntegrationArea=bwperim(fad(CurrentChannel).channels(CurrentFrame).fits.maskUsedForTotalInt);
        SnippetOverlay=cat(3,IntegrationArea/2 + ...
            +CurrentSnippet,CurrentSnippet,CurrentSnippet);
        
    
        imshow(SnippetOverlay,...
            [],'Border','Tight','InitialMagnification',1000)

        hold on
        SnippetX=(SnippetSize-1)/2+1-...
            (single(fad(CurrentChannel).channels(CurrentFrame).fits.x(CurrentParticleIndex))-...
            fad(CurrentChannel).channels(CurrentFrame).fits.x_fit(CurrentParticleIndex));
        SnippetY=(SnippetSize-1)/2+1-...
            (single(fad(CurrentChannel).channels(CurrentFrame).fits.y(CurrentParticleIndex))-...
            fad(CurrentChannel).channels(CurrentFrame).fits.y_fit(CurrentParticleIndex));
        PlotHandle=ellipse(fad(CurrentChannel).channels(CurrentFrame).fits.r_max(CurrentParticleIndex),...
            fad(CurrentChannel).channels(CurrentFrame).fits.r_min(CurrentParticleIndex),...
            fad(CurrentChannel).channels(CurrentFrame).fits.theta(CurrentParticleIndex),SnippetX,SnippetY,'g');
        set(PlotHandle,'LineWidth',2.5)
        hold off
    else
        imshow(zeros(SnippetSize))
    end
    
    
    
    
    figure(ZProfileFig)
    if (~isempty(xTrace))
        ZProfile=fad(CurrentChannel).channels(CurrentFrame).fits.shadowsDog{CurrentParticleIndex};
        ZProfileRaw=fad(CurrentChannel).channels(CurrentFrame).fits.shadowsRaw{CurrentParticleIndex};

        [Dummy,MaxZ]=max(ZProfile);
        [Dummy,MaxZRaw]=max(ZProfileRaw);
        
        plot([1:length(ZProfile)]-MaxZ+...
            double(fad(CurrentChannel).channels(CurrentFrame).fits.z(CurrentParticleIndex)),...
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
            double(fad(CurrentChannel).channels(CurrentFrame).fits.z(CurrentParticleIndex)))
        title('Z profile')

    end
       
       
    
    figure(TraceFig)
    if ~strcmp(lower(ExperimentType),'inputoutput')
        %Only update the trace information if we have switched particles
        if (CurrentParticle~=PreviousParticle)|~exist('Amp')|(CurrentChannel~=PreviousChannel)
            PreviousParticle=CurrentParticle;
            [Frames,Amp]=PlotParticleTrace(CurrentParticle,Particles{CurrentChannel},fad(CurrentChannel),FISHPath,FilePrefix);
        end
        plot(Frames(Particles{CurrentChannel}(CurrentParticle).FrameApproved),Amp(Particles{CurrentChannel}(CurrentParticle).FrameApproved),'.-k')
        hold on
        plot(Frames(~Particles{CurrentChannel}(CurrentParticle).FrameApproved),Amp(~Particles{CurrentChannel}(CurrentParticle).FrameApproved),'.r')
        plot(Frames(Frames==CurrentFrame),Amp(Frames==CurrentFrame),'ob')
        hold off
        try
            xlim([min(Frames)-1,max(Frames)+1]);
        catch
        end
        xlabel('Frames')
        h = get(gca, 'xtick');
        set(gca,'xticklabel',h)
        ylabel('Intensity (A.U.)')
    else
        %Only update the trace information if we have switched particles
        if (CurrentParticle~=PreviousParticle)|~exist('Amp')|(CurrentChannel~=PreviousChannel)
            clf
            PreviousParticle=CurrentParticle;
            [Frames,Amp]=PlotParticleTrace(CurrentParticle,Particles{CurrentChannel},fad(CurrentChannel),FISHPath,FilePrefix);
        end
        
        plot(Frames(Particles{CurrentChannel}(CurrentParticle).FrameApproved),Amp(Particles{CurrentChannel}(CurrentParticle).FrameApproved),'.-k')
        hold on
        plot(Frames(~Particles{CurrentChannel}(CurrentParticle).FrameApproved),Amp(~Particles{CurrentChannel}(CurrentParticle).FrameApproved),'.r')
        plot(Frames(Frames==CurrentFrame),Amp(Frames==CurrentFrame),'ob')
        try
            xlim([min(schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).frames),max(schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).frames)])
        catch
        end
        ax1=gca;
        set(ax1,'XAxisLocation','Bottom')
        set(ax1,'YAxisLocation','Left')
        set(ax1,'Box','off')
        set(ax1,'YColor','k','XColor','k')
        ax2= axes('Position',get(ax1,'Position'),...
                   'XAxisLocation','top',...
                   'YAxisLocation','right',...
                   'Color','none',...
                   'Box','off',...
                   'XColor','k','YColor','k');
        hold on
        plot(schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).frames,...
            max(schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).Fluo,[],2),'g.-')
        hold off
        try
            xlim([min(schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).frames),max(schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).frames)])
        catch
        end
        xlabel('Frames')
        h = get(gca, 'xtick');
        set(gca,'xticklabel',h)
        ylabel('Intensity (A.U.)')
        
    end

    
    FigureTitle=['Particle: ',num2str(CurrentParticle),'/',num2str(length(Particles{CurrentChannel})),...
        ', Frame: ',num2str(CurrentFrame),'/',num2str(TotalFrames),...
        ', Z: ',num2str(CurrentZ),'/',num2str(ZSlices),', Ch: ',num2str(CurrentChannel)];
    
    if HideApprovedFlag==1
        FigureTitle=[FigureTitle,', Showing non-flagged particles'];
    elseif HideApprovedFlag==2
        FigureTitle=[FigureTitle,', Showing disapproved particles'];
    end
    title(FigureTitle)
    
    
    %Define the windows
    figure(Overlay)
    set(gcf,'units', 'normalized', 'position',[0.01, .6, .33, .33]);
    if UseHistoneOverlay 
        figure(HisOverlayFig)
        set(gcf,'units', 'normalized', 'position',[0.01, .2, .33, .33]);
    end
    figure(TraceFig);
    set(gcf,'units', 'normalized', 'position',[0.35, .6, .2, .33]);
    figure(SnippetFig);
    set(gcf,'units', 'normalized', 'position',[0.35, .36, .2/2, .33/2]);
    figure(ZProfileFig);
    set(gcf,'units', 'normalized', 'position',[0.46, .36, .2/2, .33/2]);

    
    
    
    figure(Overlay)
    ct=waitforbuttonpress;
    cc=get(Overlay,'currentcharacter');
    cm=get(gca,'CurrentPoint');
    
    if (cc=='.')&(CurrentFrame<length(fad(CurrentChannel).channels))
        CurrentFrame=CurrentFrame+1;
        ManualZFlag=0;
        %DisplayRange=[];
    elseif (cc==',')&(CurrentFrame>1)
        CurrentFrame=CurrentFrame-1;
        ManualZFlag=0;
        %DisplayRange=[];
    elseif (cc=='>')&(CurrentFrame+5<length(fad(CurrentChannel).channels))
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
        if (floor(iJump)>0)&(iJump<length(fad(CurrentChannel).channels))
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
        if (floor(ParticleJump)>0)&(ParticleJump<=length(Particles{CurrentChannel}))
            CurrentParticle=ParticleJump;
            CurrentFrame=Particles{CurrentChannel}(CurrentParticle).Frame(1);
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
        for i=1:length(Particles{CurrentChannel})
            FirstFrame(i)=Particles{CurrentChannel}(i).Frame(1);
        end
        [Dummy,Permutations]=sort(FirstFrame);
        Particles{CurrentChannel}=Particles{CurrentChannel}(Permutations);
        
    elseif cc=='f'
        Answer=input('Are you sure you want to redo the tracking?  (y/n) ','s');
        Answer=lower(Answer);
        if Answer=='y'
            warning('HG: Not clear that this feature will work with the multiple channels')
            
            %We need to save the data
            save([DataFolder,filesep,'FrameInfo.mat'],'FrameInfo')
            if UseHistoneOverlay
                save([DataFolder,filesep,'Particles.mat'],'Particles','fad','fad2','Threshold1','Threshold2')
                save([DropboxFolder,filesep,FilePrefix(1:end-1),filesep,FilePrefix(1:end-1),'_lin.mat'],'schnitzcells')
            else
                save([DataFolder,filesep,'Particles.mat'],'Particles','fad','fad2','Threshold1','Threshold2')            
            end
        display('Particles saved.')
        if NChannels==1
            Particles=Particles{1};
        end
            
           [Particles,schnitzcells,fad,fad2]=TrackmRNADynamics(FilePrefix(1:end-1),...
               Threshold1,Threshold2); 
        if NChannels==1
            Particles={Particles};
        end
           %Check the FrameApproved field
            for i=1:length(Particles{CurrentChannel})
                if isempty(Particles{CurrentChannel}(i).FrameApproved)
                    Particles{CurrentChannel}(i).FrameApproved=logical(ones(size(Particles{CurrentChannel}(i).Frame)));
                end
            end
        end
    elseif cc=='c'
        PreviousParticle=0;
        if ~sum(Particles{CurrentChannel}(CurrentParticle).Frame==CurrentFrame)
            ConnectPosition=ginput(1);
            
            
            if ~isempty(ConnectPosition)
                [ParticleOutput,IndexOutput]=FindClickedParticle(ConnectPosition,CurrentFrame,fad(CurrentChannel),Particles{CurrentChannel});
               %Check that the clicked particle doesn't exist in a previous
                %frame, that there is no overlap of frames. If it does
                %exist in a previous frame we will have to disconnect it.
                if sum(Particles{CurrentChannel}(ParticleOutput).Frame<CurrentFrame)
                    %Disconnect the clicked particle
                    Particles{CurrentChannel}=SeparateParticleTraces(ParticleOutput,CurrentFrame,Particles{CurrentChannel});
                    ParticleOutput=ParticleOutput+1;
                end
                    
                Particles{CurrentChannel}=JoinParticleTraces(CurrentParticle,ParticleOutput,Particles{CurrentChannel});
                %Do this in case the clicked particle comes before the current
                %particle in the structure
                if ParticleOutput<CurrentParticle
                    CurrentParticle=ParticleOutput;
                end
                %Sort the frames within the particle. This is useful if we
                %connected to a particle that came before.
                [SortedFrame,Permutations]=sort(Particles{CurrentChannel}(CurrentParticle).Frame);
                Particles{CurrentChannel}(CurrentParticle).Frame=Particles{CurrentChannel}(CurrentParticle).Frame(Permutations);
                Particles{CurrentChannel}(CurrentParticle).Index=Particles{CurrentChannel}(CurrentParticle).Index(Permutations);
                Particles{CurrentChannel}(CurrentParticle).FrameApproved=Particles{CurrentChannel}(CurrentParticle).FrameApproved(Permutations);                

            end
            
        else
            ConnectPosition=ginput(1);
            [ParticleOutput,IndexOutput]=FindClickedParticle(ConnectPosition,CurrentFrame,fad(CurrentChannel),Particles{CurrentChannel});
            
            %If it's an independent particle swap it with the frame in the
            %current particle
            if (length(Particles{CurrentChannel}(ParticleOutput).Frame)==1)&...
                    (sum(Particles{CurrentChannel}(ParticleOutput).Frame==CurrentFrame)==1)
                
                ParticleTemp=Particles{CurrentChannel}(ParticleOutput);
                
                %Copy the particle out
                Particles{CurrentChannel}(ParticleOutput).Index=...
                    Particles{CurrentChannel}(CurrentParticle).Index(Particles{CurrentChannel}(CurrentParticle).Frame==CurrentFrame);

                %Copy the new particle in
                Particles{CurrentChannel}(CurrentParticle).Index(Particles{CurrentChannel}(CurrentParticle).Frame==CurrentFrame)=...
                    ParticleTemp.Index;
                Particles{CurrentChannel}(CurrentParticle).FrameApproved(Particles{CurrentChannel}(CurrentParticle).Frame==CurrentFrame)=1;
            else
                display('Cannnot connect to two particles!')
            end
            
            
        end
    elseif cc=='p' %Identify a particle. It will also tell you the particle associated with
                   %  the clicked nucleus.
        [ConnectPositionx,ConnectPositiony]=ginputc(1,'color', 'b', 'linewidth',1);
        ConnectPosition = [ConnectPositionx,ConnectPositiony];
        if ~isempty(ConnectPosition)
            %Find the closest particle
            [ParticleOutput,IndexOutput]=FindClickedParticle(ConnectPosition,CurrentFrame,fad(CurrentChannel),Particles{CurrentChannel});
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
                for i=1:length(Particles{CurrentChannel})
                    if ~isempty(Particles{CurrentChannel}(i).Nucleus)
                        AssignedNuclei(i)=Particles{CurrentChannel}(i).Nucleus;
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
        
    elseif cc=='\' %Moves to clicked particle.
        [ConnectPositionx,ConnectPositiony]=ginputc(1,'color', 'b', 'linewidth',1);
        ConnectPosition = [ConnectPositionx,ConnectPositiony];
        if ~isempty(ConnectPosition)
            %Find the closest particle
            [ParticleOutput,IndexOutput]=FindClickedParticle(ConnectPosition,CurrentFrame,fad,Particles{CurrentChannel});
            display(['Clicked particle: ',num2str(ParticleOutput)]);
            try
                ParticleJump=ParticleOutput;
            end
            if (floor(ParticleJump)>0)&(ParticleJump<=length(Particles{CurrentChannel}))
                CurrentParticle=ParticleJump;
                CurrentFrame=Particles{CurrentChannel}(CurrentParticle).Frame(1);
                ManualZFlag=0;
            end
            
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
                for i=1:length(Particles{CurrentChannel})
                    if ~isempty(Particles{CurrentChannel}(i).Nucleus)
                        AssignedNuclei(i)=Particles{CurrentChannel}(i).Nucleus;
                    else
                        AssignedNuclei(i)=nan;
                    end
                end
                AssociatedParticle=find(AssignedNuclei==ClickedSchnitz);

                if isempty(AssociatedParticle)
                    display(['Nucleus ',num2str(ClickedSchnitz),' does not have an associated particle'])
                else
                    display(['Particle ',num2str(AssociatedParticle),' is associated with nucleus ',...
                        num2str(ClickedSchnitz)])
                end
            end
            
        end    
%     elseif cc=='1'  %Give the nucleus number in the schnitzcell segmentation structure. This
%                     %only  works for troubleshooting and you need to be online and on the
%                     %Princeton network/VPN for now.
%         if exist('Z:\FISHDrosophila\Analysis\schnitzcells\Analysis')
%             ClickedPosition=round(ginput(1));
%             if ~isempty(ClickedPosition)
%                 Dashes=strfind(FilePrefix,'-');
%                 FileDate=FilePrefix(1:Dashes(3)-1);
%                 load(['Z:\FISHDrosophila\Analysis\schnitzcells\Analysis\',...
%                     FileDate,filesep,FilePrefix(1:end-1),'\segmentation\',...
%                     FilePrefix(1:end-1),'seg',iIndex(CurrentFrame,3),'.mat'],...
%                     '-mat','LN')
%                 %See which nucleus it corresponds to
%                 ClickedNucleus=LN(ClickedPosition(2),ClickedPosition(1));
%                 display(['Clicked nucleus: ',num2str(ClickedNucleus)])
%                 
%                 %Now check if there's any schnitzcell associated with this
%                 SchnitzFound=0;
%                 for i=1:length(schnitzcells)
%                     if sum((schnitzcells(i).frames)==CurrentFrame)
%                         IndexToSearch=find((schnitzcells(i).frames)==CurrentFrame);
%                         try 
%                             if schnitzcells(i).cellno(IndexToSearch)==ClickedNucleus
%                                 SchnitzFound=i;
%                             end
%                         catch
%                             display(['WARNING: Consistency error in schnitz ',num2str(i)])
%                         end
%                     end
%                 end
%                 
%                 if SchnitzFound==0
%                     display('No schnitz associated with it! Creating it...')
%                     
%                     Nschnitz=length(schnitzcells);
% 
%                     schnitzcells(Nschnitz+1).P=0;
%                     schnitzcells(Nschnitz+1).D=0;
%                     schnitzcells(Nschnitz+1).E=0;
%                     schnitzcells(Nschnitz+1).frames=CurrentFrame;
%                     schnitzcells(Nschnitz+1).cenx=Ellipses{CurrentFrame}(ClickedNucleus,1)+1;
%                     schnitzcells(Nschnitz+1).ceny=Ellipses{CurrentFrame}(ClickedNucleus,2)+1;
%                     %schnitzcells(Nschnitz+1).ang=Ellipses{CurrentFrame}(ClickedNucleus,5);
%                     schnitzcells(Nschnitz+1).len=max([Ellipses{CurrentFrame}(ClickedNucleus,3),...
%                         Ellipses{CurrentFrame}(ClickedNucleus,4)]);
%                     schnitzcells(Nschnitz+1).cellno=ClickedNucleus;
%                 else
%                     display(['It is associated with schnitz ',num2str(SchnitzFound)])
%                 end
%             end
%            
%         else
%             display('This feature only works on Hernan''s and if you''re connected to the VPN')
%         end
%                     
    elseif cc=='u'
         [x2,y2]=fad2xyzFit(CurrentFrame,fad2(CurrentChannel), 'addMargin'); 
         if ~isempty(x2)
            fad2Position=ginput(1);
            [fad(CurrentChannel),fad2(CurrentChannel),Particles{CurrentChannel}]=...
                Integratefad2Particle(fad(CurrentChannel),fad2(CurrentChannel),...
                fad2Position,Particles{CurrentChannel},CurrentFrame);
         end
    elseif cc=='i'
        PreviousParticle=0;
         [x2,y2]=fad2xyzFit(CurrentFrame,fad2(CurrentChannel), 'addMargin'); 
         if ~isempty(x2)
            fad2Position=ginput(1);
            if (~isempty(fad2Position))
                [fad(CurrentChannel),fad2(CurrentChannel),Particles{CurrentChannel}]=...
                    Integratefad2Particle(fad(CurrentChannel),fad2(CurrentChannel),fad2Position,Particles{CurrentChannel},CurrentFrame);
            end
         end

         
         if (~sum(Particles{CurrentChannel}(CurrentParticle).Frame==CurrentFrame))&(~isempty(fad2Position))
            ConnectPosition=fad2Position;
            [ParticleOutput,IndexOutput]=FindClickedParticle(ConnectPosition,CurrentFrame,fad(CurrentChannel),Particles{CurrentChannel});
            
            %Check that the clicked particle doesn't exist in a previous
            %frame, that there is no overlap of frams.  Maybe I can have
            %those in a different color.

           
            if sum(Particles{CurrentChannel}(ParticleOutput).Frame<CurrentFrame)
                display(['Target particle (',num2str(ParticleOutput),') is already in a previous frame!']);
            else
                Particles=JoinParticleTraces(CurrentParticle,ParticleOutput,Particles{CurrentChannel});
                %Do this in case the clicked particle comes before the current
                %particle in the structure
                if ParticleOutput<CurrentParticle
                    CurrentParticle=ParticleOutput;
                end
                %Sort the frames within the particle. This is useful if we
                %connected to a particle that came before.
                [SortedFrame,Permutations]=sort(Particles{CurrentChannel}(CurrentParticle).Frame);
                Particles{CurrentChannel}(CurrentParticle).Frame=Particles{CurrentChannel}(CurrentParticle).Frame(Permutations);
                Particles{CurrentChannel}(CurrentParticle).Index=Particles{CurrentChannel}(CurrentParticle).Index(Permutations);
                Particles{CurrentChannel}(CurrentParticle).FrameApproved=Particles{CurrentChannel}(CurrentParticle).FrameApproved(Permutations);
                
                if UseHistoneOverlay
                    %Check for consistency within schnitzcell 
                    [Particles{CurrentChannel},schnitzcells]=CheckSchnitzLineage(Particles{CurrentChannel},CurrentParticle,schnitzcells,CurrentFrame,...
                        Overlay);
                end
                
                
                
                
            end
            
        else
            display('Cannnot connect to two particles!')
        end
         
         
         
         
     elseif cc=='d'  %d Separate traces. The separated particle won't have a nucleus assigned!
         PreviousParticle=0;
        %Check that the particle does actually exist in this frame
        if ~(Particles{CurrentChannel}(CurrentParticle).Frame(1)==CurrentFrame)
            if sum(Particles{CurrentChannel}(CurrentParticle).Frame==CurrentFrame)
                Particles{CurrentChannel}=SeparateParticleTraces(CurrentParticle,CurrentFrame,Particles{CurrentChannel});
            end
        elseif length(Particles{CurrentChannel}(CurrentParticle).Frame)==1
            Particles{CurrentChannel}(CurrentParticle).Nucleus=[];
        else
            display('Cannot divide a trace at the first time point')
        end
    elseif cc=='q'      %Approve a trace
        if Particles{CurrentChannel}(CurrentParticle).Approved==1
            Particles{CurrentChannel}(CurrentParticle).Approved=2;
        elseif Particles{CurrentChannel}(CurrentParticle).Approved==0
            Particles{CurrentChannel}(CurrentParticle).Approved=1;
        elseif Particles{CurrentChannel}(CurrentParticle).Approved==2
            Particles{CurrentChannel}(CurrentParticle).Approved=0;
        end
    elseif cc=='w'      %Disapproove a trace
        if Particles{CurrentChannel}(CurrentParticle).Approved==-1
            Particles{CurrentChannel}(CurrentParticle).Approved=0;
        else
            Particles{CurrentChannel}(CurrentParticle).Approved=-1;
        end    
        
    elseif cc=='s'
        
        %If we only have one channel bring Particles back to the legacy
        %format without any cells
        if NChannels==1
            Particles=Particles{1};
        end
        
        save([DataFolder,filesep,'FrameInfo.mat'],'FrameInfo')
        if UseHistoneOverlay
            save([DataFolder,filesep,'Particles.mat'],'Particles','fad','fad2','Threshold1','Threshold2')
            save([DropboxFolder,filesep,FilePrefix(1:end-1),filesep,FilePrefix(1:end-1),'_lin.mat'],'schnitzcells')
        else
            save([DataFolder,filesep,'Particles.mat'],'Particles','fad','fad2','Threshold1','Threshold2')            
        end
        display('Particles saved.')
        if NChannels==1
            Particles={Particles};
        end
        
    elseif cc=='t'
        ShowThreshold2=~ShowThreshold2;
    elseif (cc=='y')&(~UseHistoneOverlay)
            FrameInfo=DetermineNC(fad,Particles{CurrentChannel},FrameInfo);
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
    elseif (cc=='m')&(CurrentParticle<length(Particles{CurrentChannel}))
        
        NextParticle=CurrentParticle+1;
        
        if NextParticle>length(Particles{CurrentChannel})
            NextParticle=length(Particles{CurrentChannel});
        end
        
        
        %Mode 1 - skip approved or flagged traces
        while (HideApprovedFlag)==1&(NextParticle<length(Particles{CurrentChannel}))&...
                ((Particles{CurrentChannel}(NextParticle).Approved==1)|(Particles{CurrentChannel}(NextParticle).Approved==-1)|...
                (Particles{CurrentChannel}(NextParticle).Approved==2))
                NextParticle=NextParticle+1;
        end
        
        %Mode 2 - skip approved traces
        while ((HideApprovedFlag)==2)&(NextParticle<length(Particles{CurrentChannel}))&...
                ((Particles{CurrentChannel}(NextParticle).Approved==1)|(Particles{CurrentChannel}(NextParticle).Approved==2))
                NextParticle=NextParticle+1;
        end
        
        
        CurrentParticle=NextParticle;
%         if ~isempty(Particles(CurrentParticle).Nucleus)
%             CurrentFrame=schnitzcells(Particles(CurrentParticle).Nucleus).frames(1)-1;
%         else
%             CurrentFrame=Particles(CurrentParticle).Frame(1);
%         end
        CurrentFrame=Particles{CurrentChannel}(CurrentParticle).Frame(1);

        ParticleToFollow=[];
        DisplayRange=[];
        
        
        display('Missing frames:')
        Particles{CurrentChannel}(CurrentParticle).Frame(find(diff(Particles{CurrentChannel}(CurrentParticle).Frame)>1))
        
    elseif (cc=='n')&(CurrentParticle>1)
        Approved=(find([Particles{CurrentChannel}.Approved]));
        %NotApproved=(find(~[Particles.Approved]));
        
        NextParticle=CurrentParticle-1;
        
        
        
        %Mode 1 - show non-flagged traces
        while (HideApprovedFlag)==1&(NextParticle>1)&...
                ((Particles{CurrentChannel}(NextParticle).Approved==1)|(Particles{CurrentChannel}(NextParticle).Approved==-1)|...
                (Particles{CurrentChannel}(NextParticle).Approved==2))
            NextParticle=NextParticle-1;
            if NextParticle<1
                NextParticle=1;
            end
        end
        
        
        %Mode 2 - show dissapproved traces
        while ((HideApprovedFlag)==2)&(NextParticle>1)&...
                ((Particles{CurrentChannel}(NextParticle).Approved==1)|(Particles{CurrentChannel}(NextParticle).Approved==2))
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
        CurrentFrame=Particles{CurrentChannel}(CurrentParticle).Frame(1);

        
        ParticleToFollow=[];
        
        DisplayRange=[];
        
    elseif cc=='e'
        Particles{CurrentChannel}(CurrentParticle).FrameApproved(Particles{CurrentChannel}(CurrentParticle).Frame==CurrentFrame)=...
            ~Particles{CurrentChannel}(CurrentParticle).FrameApproved(Particles{CurrentChannel}(CurrentParticle).Frame==CurrentFrame);
    
    
    
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
                if Particles{CurrentChannel}(CurrentParticle).Nucleus==ClickedSchnitz %Split the lineage
                    [Particles{CurrentChannel},schnitzcells]=SplitSchnitz(Particles{CurrentChannel},schnitzcells,...
                            CurrentFrame,...
                            CurrentParticle);
                else
                    try
                        [Particles{CurrentChannel},schnitzcells]=JoinSchnitz(Particles{CurrentChannel},schnitzcells,Particles(CurrentParticle).Nucleus,...  
                            ClickedSchnitz,CurrentFrame);
                    catch
                        display('Error in JoinSchnitz')
                    end
                end
            elseif length(ClickedSchnitz)==2
                if ClickedSchnitz(1)~=ClickedSchnitz(2)
                    [Particles{CurrentChannel},schnitzcells]=SplitSchnitzDaughters(Particles{CurrentChannel},schnitzcells,...
                            CurrentFrame,...
                            Particles{CurrentChannel}(CurrentParticle).Nucleus,ClickedSchnitz(1),ClickedSchnitz(2));
                else
                    [Particles{CurrentChannel},schnitzcells]=SplitSchnitzDaughters(Particles{CurrentChannel},schnitzcells,...
                        CurrentFrame,...
                        Particles{CurrentChannel}(CurrentParticle).Nucleus,ClickedSchnitz(1),0);
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
                schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).P=ClickedSchnitz;

            elseif length(ClickedSchnitz)==0
                schnitzcells(Particles{CurrentChannel}(CurrentParticle).Nucleus).P=0;
            end                
        end    
        
    elseif cc=='8'      %Switch channels
        PreviousChannel=CurrentChannel;
        CurrentChannel=CurrentChannel+1;
        if CurrentChannel>NChannels
            CurrentChannel=1;
        end
        
        %If a particle is associated with this same nucleus in the new
        %channel then change to it
        AssignedNucleusPreviousChannel=Particles{PreviousChannel}(CurrentParticle).Nucleus;
        
        %Now, find its associated particle
        AssignedNucleusNewChannel=[];
        for i=1:length(Particles{CurrentChannel})
            if ~isempty(Particles{CurrentChannel}(i).Nucleus)
                AssignedNucleusNewChannel(i)=Particles{CurrentChannel}(i).Nucleus;
            else
                AssignedNucleusNewChannel(i)=nan;
            end
        end
        
        if ~isempty(find(AssignedNucleusNewChannel==AssignedNucleusPreviousChannel))
            CurrentParticle=find(AssignedNucleusNewChannel==AssignedNucleusPreviousChannel);
        end
        
        
    elseif cc=='0'      %Debugging mode
        keyboard;
    end
        
end


save([DataFolder,filesep,'FrameInfo.mat'],'FrameInfo')


%If we only have one channel bring Particles back to the legacy
%format without any cells
if NChannels==1
    Particles=Particles{1};
end


if UseHistoneOverlay
    save([DataFolder,filesep,'Particles.mat'],'Particles','fad','fad2','Threshold1','Threshold2')
    save([DropboxFolder,filesep,FilePrefix(1:end-1),filesep,FilePrefix(1:end-1),'_lin.mat'],'schnitzcells')
else
    save([DataFolder,filesep,'Particles.mat'],'Particles','fad','fad2','Threshold1','Threshold2')            
end
close all
display('Particles saved.')
display(['(Left off at Particle #', num2str(CurrentParticle), ')'])





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
