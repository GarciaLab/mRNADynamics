function CheckNuclearActivity(varargin)

%Determines the frame to use to count active and inactive nuclei. The user
%needs to pick a frame and then select the nuclei that should not be
%counted


%Usage:

%Left click: Toggle nucleus for counting

%Frame specific:
%. Move a frame forward
%, Move a frame backwards
%> Move five frames forward
%< Move five frames backwards
%m Move to the next nc
%n Move to the previous nc
%a Move up in Z
%z Move down in z
%j Jump to a specified frame
%g,b Increase/decrease histone channel contrast
%q Set this channel for counting active and inactive nuclei



%General:
%s Save the current Particles structure
%0 Enter debug mode to fix things manually

close all


%% Information about about folders

    
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,SchnitzcellsFolder]=...
    DetermineLocalFolders(varargin{1});


%%

Prefix=varargin{1};
FilePrefix=[Prefix,'_'];
load([DropboxFolder,filesep,Prefix,filesep,'Particles.mat'])
load([DropboxFolder,filesep,Prefix,filesep,'CompiledParticles.mat'])
load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'])


%Create an mask for the AP bins that were rejected

%Get the information about the AP axis as well as the image shifts
%used for the stitching of the two halves of the embryo
load([DropboxFolder,filesep,Prefix,filesep,'\APDetection.mat'])

%Angle between the x-axis and the AP-axis
APAngle=atan((coordPZoom(2)-coordAZoom(2))/(coordPZoom(1)-coordAZoom(1)));
APLength=sqrt((coordPZoom(2)-coordAZoom(2))^2+(coordPZoom(1)-coordAZoom(1))^2);

Image=imread([FISHPath,'\Data\',FilePrefix(1:end-1),filesep,...
            FilePrefix,iIndex(1,3),'_z',iIndex(10,2),'.tif']);

APPosImage=zeros(size(Image));
[Rows,Columns]=size(Image);

for i=1:Rows
    for j=1:Columns
        Angle=atan((i-coordAZoom(2))./(j-coordAZoom(1)));
        Distance=sqrt((coordAZoom(2)-i).^2+(coordAZoom(1)-j).^2);
        APPosition=Distance.*cos(Angle-APAngle);
        APPosImage(i,j)=APPosition/APLength;
    end
end


%Bin the pixels along the AP axis
APPosBinImage=zeros(size(APPosImage));
for i=1:(length(APbinID)-1)
    if isnan(APbinArea(i))
        FilteredMask=(APbinID(i)<=APPosImage)&(APbinID(i+1)>APPosImage);

        APPosBinImage=APPosBinImage+FilteredMask;
    end
end






%Check if we have the histone channel and we have done the nuclear
%segmentation.
if exist([FISHPath,'\Data\',FilePrefix(1:end-1),filesep,...
        FilePrefix(1:end-1),'-His_',iIndex(1,3),'.tif'])
    load([DropboxFolder,filesep,FilePrefix(1:end-1),'\Ellipses.mat'])
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
    


%Load the information about the nc from the XLS file
[Num,Txt]=xlsread([DropboxFolder,'\HGMovieDatabase.xlsx']);
XLSHeaders=Txt(1,:);
Txt=Txt(2:end,:);

%Find the different columns.
DataFolderColumn=7;

%Convert the prefix into the string used in the XLS file
Dashes=findstr(FilePrefix,'-');

%Find the corresponding entry in the XLS file
XLSEntry=find(strcmp(Txt(:,DataFolderColumn),...
    ['D:\Hernan\LivemRNA\Analysis\MS2\MCPNoNLS+MS2\',...
    FilePrefix(1:Dashes(3)-1),filesep,FilePrefix(Dashes(3)+1:end-1)]));



if strcmp(Txt(XLSEntry,4),'Yes')

    nc9=Num(XLSEntry,5);
    nc10=Num(XLSEntry,6);
    nc11=Num(XLSEntry,7);
    nc12=Num(XLSEntry,8);
    nc13=Num(XLSEntry,9);
    nc14=Num(XLSEntry,10);
    CF=Num(XLSEntry,11);


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
    display('No nc information provided in the XLS file.')
end

%Check if we have already determined nc
if (~isfield(FrameInfo,'nc'))&(~UseHistoneOverlay)
    FrameInfo=DetermineNC(fad,Particles,FrameInfo);
elseif UseSchnitz
    load([DropboxFolder,filesep,FilePrefix(1:end-1),filesep,FilePrefix(1:end-1),'_lin.mat'])
end


TotalFrames=length(ElapsedTime);


if ~exist([DropboxFolder,filesep,FilePrefix(1:end-1),filesep,'CountEllipses.mat'])

    %Generate the ellipse filter for all frames
    EdgeWidth=10;
    for i=1:TotalFrames
        CurrentEllipses=Ellipses{i};

        Radius=max(CurrentEllipses(:,3:4)')';

        EllipseFilter{i}=(CurrentEllipses(:,1)-Radius-EdgeWidth>0)&...
            (CurrentEllipses(:,2)-Radius-EdgeWidth>0)&...
            (CurrentEllipses(:,1)+Radius+EdgeWidth<FrameInfo(1).PixelsPerLine)&...
            (CurrentEllipses(:,2)+Radius+EdgeWidth<FrameInfo(1).LinesPerFrame);
        
        %Remove those nuclei tagged with code 2 or -1. These should not be used
        %at all

        schnitzCellNo=[];
        for j=1:length(Particles)
            if (Particles(j).Approved==2)|(Particles(j).Approved==-1)
                if ~isempty(Particles(j).Nucleus)
                    schnitzIndex=find((schnitzcells(Particles(j).Nucleus).frames-1)==i);
                    schnitzCellNo=[schnitzCellNo,schnitzcells(Particles(j).Nucleus).cellno(schnitzIndex)];
                end
            end
        end
           
        EllipseFilter{i}(schnitzCellNo)=0;
    end

    %Initialize the vector telling us which frame for each nc is used to do
    %this calculation
    CountEllipsesFrame=nan(14,1);
    for i=1:14
        if eval(['exist(''nc',num2str(i),''')'])
            CountEllipsesFrame(i)=eval(['nc',num2str(i)]);
        end
    end
else
    EdgeWidth=10;
    load([DropboxFolder,filesep,FilePrefix(1:end-1),filesep,'CountEllipses.mat'])
end


%Some flags and initial parameters
ZSlices=FrameInfo(1).NumberSlices+2;   %Note that the blank slices are included
CurrentZ=round(ZSlices/2);          
ManualZFlag=0;
CurrentNC=14;
CurrentFrame=CountEllipsesFrame(14);


DisplayRange=[];



Overlay=figure;

if UseHistoneOverlay
    HisOverlayFig=figure;
end


cc=1;

while (cc~=13)
    EllipseHandle=[];
    EllipseHandleYellow=[];
    EllipseHandleBlue=[];
    EllipseHandleWhite=[];
    EllipseHandleGreen=[];
    
    %Get the coordinates taking the margins into account
    [x,y]=fad2xyzFit(CurrentFrame,fad, 'addMargin'); 
    
    %Use the approved particles from CompiledParticles
    ApprovedParticles=[CompiledParticles.OriginalParticle];
    
   
    
    %These are the positions of all the approved and disapprooved particles
    %Find the particles in this frame
    IndexApprovedParticles=[];
    for i=1:length(ApprovedParticles)
        if sum(Particles(ApprovedParticles(i)).Frame==CurrentFrame)
            IndexApprovedParticles=[IndexApprovedParticles,...
                Particles(ApprovedParticles(i)).Index(Particles(ApprovedParticles(i)).Frame==CurrentFrame)];
        end
    end
    xApproved=x(IndexApprovedParticles);
    yApproved=y(IndexApprovedParticles);

    try
        Image=imread([FISHPath,'\Data\',FilePrefix(1:end-1),filesep,...
            FilePrefix,iIndex(CurrentFrame,3),'_z',iIndex(CurrentZ,2),'.tif']);
        
        %Create an overlay with the APbin image
        ImageBins=cat(3,mat2gray(Image),mat2gray(Image),mat2gray(Image)+APPosBinImage/3);
    catch
        display(['Warning: Could not load file: ',...
            FilePrefix,iIndex(CurrentFrame,3),'_z',iIndex(CurrentZ,2),'.tif'])
    end
    

    figure(Overlay)
    imshow(ImageBins,'Border','Tight')
    hold on
    plot(x,y,'or')
    plot(xApproved,yApproved,'ob')
    hold off
    if ~UseHistoneOverlay
        set(gcf,'Position',[10    89   677   599]);
        set(gcf,'MenuBar','none','ToolBar','none')
    else
        set(gcf,'Position',[10   414   677   351]);
        set(gcf,'MenuBar','none','ToolBar','none')
    end
    
    set(gcf,'Name',['Frame: ',num2str(CurrentFrame),'/',num2str(TotalFrames),...
        ', Z: ',num2str(CurrentZ),'/',num2str(ZSlices),' nc: ', num2str(CurrentNC),...
        '. ',num2str(round(ElapsedTime(CurrentFrame)-ElapsedTime(eval(['nc',num2str(CurrentNC)])))),...
        ' min']);
    
    
    
    
    if UseSchnitz
        %Show all the nuclei
        hold on       
        EllipseHandle=notEllipse(Ellipses{CurrentFrame}(:,3),...
            Ellipses{CurrentFrame}(:,4),...
            Ellipses{CurrentFrame}(:,5),...
            Ellipses{CurrentFrame}(:,1)+1,...
            Ellipses{CurrentFrame}(:,2)+1,'r',50);
        hold off
            

        %Show the ones that have been approved
        hold on
        schnitzCellNo=[];
        for i=1:length(CompiledParticles)
            schnitzIndex=find((schnitzcells(CompiledParticles(i).Nucleus).frames-1)==CurrentFrame);
            schnitzCellNo=[schnitzCellNo,schnitzcells(CompiledParticles(i).Nucleus).cellno(schnitzIndex)];
        end
        
        EllipseHandle=notEllipse(Ellipses{CurrentFrame}(schnitzCellNo,3),...
            Ellipses{CurrentFrame}(schnitzCellNo,4),...
            Ellipses{CurrentFrame}(schnitzCellNo,5),...
            Ellipses{CurrentFrame}(schnitzCellNo,1)+1,...
            Ellipses{CurrentFrame}(schnitzCellNo,2)+1,'b',50);
        hold off
     
        %Show the ones flagged with 2. These are not in CompileParticles
        %nor should they be used for anything.
        schnitzCellNo=[];
        for i=1:length(Particles)
            if Particles(i).Approved==2
                schnitzIndex=find((schnitzcells(Particles(i).Nucleus).frames-1)==CurrentFrame);
                schnitzCellNo=[schnitzCellNo,schnitzcells(Particles(i).Nucleus).cellno(schnitzIndex)];
            end
        end
            
        Ellipse2Handle=notEllipse(Ellipses{CurrentFrame}(schnitzCellNo,3),...
            Ellipses{CurrentFrame}(schnitzCellNo,4),...
            Ellipses{CurrentFrame}(schnitzCellNo,5),...
            Ellipses{CurrentFrame}(schnitzCellNo,1)+1,...
            Ellipses{CurrentFrame}(schnitzCellNo,2)+1,'c',50);    
         
    end
    

    %Add the rectangle we'll use for counting
    hold on
    rectangle('Position',[EdgeWidth,EdgeWidth,...
        size(Image,2)-EdgeWidth*2,size(Image,1)-EdgeWidth*2],...
        'EdgeColor','r','LineWidth',2);
    hold off
    
    %Show the ellipses that pass the filter

    
    hold on
    EllipseFillHandle=notEllipse(Ellipses{CurrentFrame}(EllipseFilter{CurrentFrame},3)*0.75,...
            Ellipses{CurrentFrame}(EllipseFilter{CurrentFrame},4)*0.75,...
            Ellipses{CurrentFrame}(EllipseFilter{CurrentFrame},5),...
            Ellipses{CurrentFrame}(EllipseFilter{CurrentFrame},1)+1,...
            Ellipses{CurrentFrame}(EllipseFilter{CurrentFrame},2)+1,'y',50);
    hold off
    
    
    %Show whether this is the frame we'll be using
    if CountEllipsesFrame(CurrentNC)==CurrentFrame
        set(gcf,'Color','g')
    else
        set(gcf,'Color','default')
    end
    

    if UseHistoneOverlay
        figure(HisOverlayFig)
        ImageHis=imread([FISHPath,'\Data\',FilePrefix(1:end-1),filesep,...
            FilePrefix(1:end-1),'-His_',iIndex(CurrentFrame,3),'.tif']);
        

        if isempty(DisplayRange)
            HisOverlayImage=cat(3,mat2gray(ImageHis)*2,mat2gray(Image),zeros(size(Image)));
        else
            HisOverlayImage=cat(3,mat2gray(ImageHis,double(DisplayRange))*2,mat2gray(Image),zeros(size(Image)));
        end
        imshow(HisOverlayImage,'Border','Tight')

        set(gcf,'Position',[-7    60   691   311])
        set(gcf,'MenuBar','none','ToolBar','none')
        
       
        hold on
        plot(x,y,'ow')
        plot(xApproved,yApproved,'ob')
        
        hold off
        

  
        
        if UseSchnitz
            
            copyobj(EllipseHandle,gca)
            copyobj(Ellipse2Handle,gca)
            copyobj(EllipseHandleBlue,gca)
            copyobj(EllipseHandleGreen,gca)
            copyobj(EllipseHandleWhite,gca)
            copyobj(EllipseHandleYellow,gca)

            
            

        end
    
        set(gcf,'Name',['Frame: ',num2str(CurrentFrame),'/',num2str(TotalFrames),...
            ', Z: ',num2str(CurrentZ),'/',num2str(ZSlices),' nc: ', num2str(CurrentNC),...
            '. ',num2str(round(ElapsedTime(CurrentFrame)-ElapsedTime(eval(['nc',num2str(CurrentNC)])))),...
            ' min'])
     
    end
    
    
   
    
    figure(Overlay)
    ct=waitforbuttonpress;
    cc=get(Overlay,'currentcharacter');
    cm=get(gca,'CurrentPoint');
    
    if (ct==0)&(strcmp(get(Overlay,'SelectionType'),'normal'))
        cc=1;
        if (cm(1,2)>0)&(cm(1,1)>0)&(cm(1,2)<=size(Image,1))&(cm(1,1)<=size(Image,2))
            %Find out which ellipses we clicked on so we can delete it
            
            %(x, y, a, b, theta, maxcontourvalue, time, particle_id)
            Distances=sqrt((Ellipses{CurrentFrame}(:,1)-cm(1,1)).^2+...
                (Ellipses{CurrentFrame}(:,2)-cm(1,2)).^2);
            [MinValue,MinIndex]=min(Distances);

            EllipseFilter{CurrentFrame}(MinIndex)=~EllipseFilter{CurrentFrame}(MinIndex);
            
        end
    elseif (cc=='.')&(CurrentFrame<TotalFrames)
        CurrentFrame=CurrentFrame+1;
        ManualZFlag=0;
        %DisplayRange=[];
    elseif (cc==',')&(CurrentFrame>1)
        CurrentFrame=CurrentFrame-1;
        ManualZFlag=0;
        %DisplayRange=[];
    elseif (cc=='>')&((CurrentFrame+5)<TotalFrames)
        CurrentFrame=CurrentFrame+5;
        ManualZFlag=0;
        %DisplayRange=[];
    elseif (cc=='<')&((CurrentFrame-5)>0)
        CurrentFrame=CurrentFrame-5;
        ManualZFlag=0;
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
    elseif cc=='g'      %Increase histone channel contrast
        if isempty(DisplayRange)
            DisplayRange=[min(min(ImageHis)),max(max(ImageHis))/1.5];
        else
            DisplayRange=[DisplayRange(1),DisplayRange(2)/1.5];
        end
        
    elseif cc=='b'      %Decrease histone channel contrast
        DisplayRange=[min(min(ImageHis)),max(max(ImageHis))*1.5];
 
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
                            schnitzcells(j).cenx(find((schnitzcells(j).frames-1)==CurrentFrame))];
                        yPosSuspect=[yPosSuspect,...
                            schnitzcells(j).ceny(find((schnitzcells(j).frames-1)==CurrentFrame))];
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
  
    elseif cc=='s'
        save([DropboxFolder,filesep,FilePrefix(1:end-1),filesep,'CountEllipses.mat'],...
            'CountEllipsesFrame','EllipseFilter')
        display('CountEllipses.mat saved')
    elseif (cc=='m')&(CurrentNC<14)
        
        CurrentNC=CurrentNC+1;
        CurrentFrame=CountEllipsesFrame(CurrentNC);
        
        
        DisplayRange=[];
     
        
        
    elseif (cc=='n')&CurrentNC>8
        if eval(['nc',num2str(CurrentNC-1)])>0
        
            CurrentNC=CurrentNC-1;
            CurrentFrame=CountEllipsesFrame(CurrentNC);
     
            DisplayRange=[];
        end

            
    elseif cc=='0'      %Debugging mode
        keyboard;
   
    elseif cc=='q'
        CountEllipsesFrame(CurrentNC)=CurrentFrame;
    end
end
        


save([DropboxFolder,filesep,FilePrefix(1:end-1),filesep,'CountEllipses.mat'],...
    'CountEllipsesFrame','EllipseFilter')
display('CountEllipses.mat saved')




