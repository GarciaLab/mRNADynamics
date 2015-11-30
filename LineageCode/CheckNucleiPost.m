function  schnitzcells=CheckNuclei(Prefix)
% GUI to select the daughter nuclei and manually curate the lineage


%% Extract all data
%Get the folders, including the default Dropbox one
[SourcePath,FISHPath,DefaultDropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders;
%Now get the actual DropboxFolder
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders(Prefix);



%Determine division times
%Load the information about the nc from the XLS file
[Num,Txt,XLSRaw]=xlsread([DefaultDropboxFolder,filesep,'MovieDatabase.xlsx']);
XLSHeaders=Txt(1,:);
Txt=Txt(2:end,:);

ExperimentTypeColumn=find(strcmp(XLSRaw(1,:),'ExperimentType'));
ExperimentAxisColumn=find(strcmp(XLSRaw(1,:),'ExperimentAxis'));

DataFolderColumn=find(strcmp(XLSRaw(1,:),'DataFolder'));
Dashes=findstr(Prefix,'-');

PrefixRow=find(strcmp(XLSRaw(:,DataFolderColumn),[Prefix(1:Dashes(3)-1),'\',Prefix(Dashes(3)+1:end)]));
if isempty(PrefixRow)
    PrefixRow=find(strcmp(XLSRaw(:,DataFolderColumn),[Prefix(1:Dashes(3)-1),'/',Prefix(Dashes(3)+1:end)]));
    if isempty(PrefixRow)
        error('Could not find data set in MovieDatabase.XLSX. Check if it is defined there.')
    end
end

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
Channel1Column=find(strcmp(XLSRaw(1,:),'Channel1'));
Channel2Column=find(strcmp(XLSRaw(1,:),'Channel2'));


%Find the corresponding entry in the XLS file
if (~isempty(findstr(Prefix,'Bcd')))&&(isempty(findstr(Prefix,'BcdE1')))&...
        (isempty(findstr(Prefix,'NoBcd')))&(isempty(findstr(Prefix,'Bcd1x')))&(isempty(findstr(Prefix,'Bcd4x')))
    warning('This step in CheckParticleTracking will most likely have to be modified to work')
    XLSEntry=find(strcmp(XLSRaw(:,DataFolderColumn),...
        [Date,'\BcdGFP-HisRFP']));
else
    XLSEntry=find(strcmp(XLSRaw(:,DataFolderColumn),...
        [Prefix(1:Dashes(3)-1),'\',Prefix(Dashes(3)+1:end)]));
    
    if isempty(XLSEntry)
        XLSEntry=find(strcmp(XLSRaw(:,DataFolderColumn),...
            [Prefix(1:Dashes(3)-1),'/',Prefix(Dashes(3)+1:end)]));
        if isempty(XLSEntry)
            display('%%%%%%%%%%%%%%%%%%%%%')
            error('Dateset could not be found. Check MovieDatabase.xlsx')
            display('%%%%%%%%%%%%%%%%%%%%%')
        end
    end
end


if (strcmp(XLSRaw(XLSEntry,Channel2Column),'His-RFP'))||...
        (strcmp(XLSRaw(XLSEntry,Channel1Column),'His-RFP'))||...
        (strcmp(XLSRaw(XLSEntry,Channel2Column),'MCP-TagRFP(1)'))||...
        (strcmp(XLSRaw(XLSEntry,Channel1Column),'MCP-mCherry(3)'))||...
        (strcmp(XLSRaw(XLSEntry,Channel2Column),'MCP-mCherry(3)'))
    nc9=XLSRaw{XLSEntry,nc9Column};
    nc10=XLSRaw{XLSEntry,nc10Column};
    nc11=XLSRaw{XLSEntry,nc11Column};
    nc12=XLSRaw{XLSEntry,nc12Column};
    nc13=XLSRaw{XLSEntry,nc13Column};
    nc14=XLSRaw{XLSEntry,nc14Column};
    CF=XLSRaw{XLSEntry,CFColumn};
end

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


%Convert the prefix into the string used in the XLS file
Dashes=findstr(Prefix,'-');

%Find the corresponding entry in the XLS file
if (~isempty(findstr(Prefix,'Bcd')))&&(isempty(findstr(Prefix,'BcdE1')))&&...
        (isempty(findstr(Prefix,'NoBcd')))&&(isempty(findstr(Prefix,'Bcd1')))&&(isempty(findstr(Prefix,'Bcd4x')))
    XLSEntry=find(strcmp(Txt(:,DataFolderColumn),...
        [Date,'\BcdGFP-HisRFP']));
else
    XLSEntry=find(strcmp(Txt(:,DataFolderColumn),...
        [Prefix(1:Dashes(3)-1),filesep,Prefix(Dashes(3)+1:end)]));
end
ncs=[nc9,nc10,nc11,nc12,nc13,nc14];

if (length(find(isnan(ncs)))==length(ncs))||(length(ncs)<6)
    error('Have the ncs been defined in MovieDatabase.XLSX?')
end

%Now do the nuclear segmentation and lineage tracking. This should be put
%into an independent function.

%Create the cell array with the names.
D=dir([PreProcPath,filesep,Prefix,filesep,'*His*.tif']);
for i=1:length(D)
    names{i}=[PreProcPath,filesep,Prefix,filesep,D(i).name];
end

%% Extract images and structures
DataFolder=[DropboxFolder,filesep,Prefix];
FilePrefix=[DataFolder(length(DropboxFolder)+2:end),'_'];

try
    dummy0=load([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat']);
    linbackup=[DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin_backup.mat'];
    schnitzcells=dummy0.schnitzcells;
    save(linbackup,'schnitzcells');
    
catch
    error('Unable to load the lineage file. Cannot Proceed');
end

% Create empty structures for Particles and FrameInfo
Particles=struct;
FrameInfo=struct;
cpexist=true; % Does Compiled Particles Exist?
try
    dummy1=load([DropboxFolder,filesep,Prefix,filesep,'Particles.mat']);
    Particles=dummy1.Particles;
catch
    warning('No Particles.mat found. No changes will be made to the Particles structure.');
    cpexist=false;
end
try
    dummy2=load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat']);
    FrameInfo=dummy2.FrameInfo;
catch
    error('No Frame Info found!');
end


dummy=FrameInfo(1);
dims(1)=dummy.LinesPerFrame;
dims(2)=dummy.PixelsPerLine;


%% Determine Order of Nuclei - Optional
% This functionality introduces a different way of going through the
% nuclei. I haven't used it, and it may be very inefficient
% Nuclei=zeros(6,length(schnitzcells));
% count=zeros(6,1);
% for i=1:length(schnitzcells)
%     for ii=1:ncs
%         if abs(max(schnitzcells(i).frames)-ncs(ii))<3
%             Nuclei(ii,count(ii))=i;
%             count(ii)=count(ii)+1;
%         end
%     end
% end
% for ii=1:ncs
%     keep=Nuclei(ii,:)~=0;
%     Nuclei(ii,:)=Nuclei(ii,keep);
% end
%% PreProcessing
% Optional step to include pre-processing
% schnitzcells=FixLineage(Prefix, FrameInfo);

% Pre-processing to create families. This will be updated at each frame
% The constant updating is inefficient, should place differently later.
len=size(schnitzcells,2);
count=0;
families=[];
for i=1:len
    if (schnitzcells(1,i).D>0)
        count=count+1;
        if ~isempty(schnitzcells(1,i).E) && (schnitzcells(1,i).E~=0)
            families(count,:)=[i schnitzcells(1,i).E schnitzcells(1,i).D]';
        else
            families(count,:)=[i 0 schnitzcells(1,i).D]';
        end
    elseif (schnitzcells(1,i).E>0)
        count=count+1;
        families(count,:)=[i schnitzcells(1,i).E 0]';
    end
end
%% GUI

keeptracking=true;
ncexist=find(ncs~=0 & ~isnan(ncs)); % Logical array to tell if the Nuclear cycle exists
Current_Nuclear_Cycle=1;

% The user can choose to cycle between Nuclear Cycles or Frames, so these
% variables keep track of the Current Nuclear Cycle and the Current Frame
CurrentNuc=1;
CurrentFrame=ncs(ncexist(Current_Nuclear_Cycle))-1;

% This is an array that stores a random colour for each family
colours=zeros(length(families),3);

while keeptracking
    
    %%% Family Checking
    % The code below finds the families once again. This is inefficient,
    % consider changing
    len=size(schnitzcells,2);
    count=0;
    families=[];
    for i=1:len
        if (schnitzcells(1,i).D>0)
            count=count+1;
            if ~isempty(schnitzcells(1,i).E) && (schnitzcells(1,i).E~=0)
                families(count,:)=[i schnitzcells(1,i).E schnitzcells(1,i).D]';
            else
                families(count,:)=[i 0 schnitzcells(1,i).D]';
            end
        elseif (schnitzcells(1,i).E>0)
            count=count+1;
            families(count,:)=[i schnitzcells(1,i).E 0]';
        end
    end
    
    % displayFrame is the function that generates the display
    displayFrame(CurrentFrame,[]);
    
    % Current_NC is only for display. It will not be used computationally
    Current_NC=Current_Nuclear_Cycle+7+min(ncexist);
    display(Current_NC);
    
    % Instructions
    display('x to save , to move frame backward . to move frame forward mouseclick to start');
    
    % This stores the variable to see if the current frame is the default
    % frame - it is used to toggle whether the user is moving between
    % frames
    frameReset=true;
    
    % cc, ct and cm gather the current user input (cm is not used yet)
    ct=waitforbuttonpress;
    cc=get(gcf,'currentcharacter')
    cm=get(gcf,'CurrentPoint');
    display('Click to start Lineaging (Mitosis Selection). ');
    
    % Delete mode
    if cc=='d'
        cc3='';
        while cc3~='n'
            ct3=waitforbuttonpress
            cc3=get(gcf,'currentcharacter');
            cm3=get(gcf,'CurrentPoint');
            toDelete=identify(cm3,CurrentFrame);
            childD=schnitzcells(toDelete).D;
            childE=schnitzcells(toDelete).E;
            if ~isempty(childD)
                schnitzcells(childD).P=[];
            end
            if ~isempty(childE)
                schnitzcells(childE).P=[];
            end
        end
    end
    
    % Save and quit
    if cc=='x'
        save([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'],'schnitzcells');
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
    elseif cc=='w'
        CurrentFrame=CurrentFrame+1;
        frameReset=false;
    elseif cc=='W'
        CurrentFrame=CurrentFrame+5;
        frameReset=false;
    elseif cc==','
        CurrentFrame=CurrentFrame-1;
        frameReset=false;
    elseif cc=='q'
        CurrentFrame=CurrentFrame-1;
        frameReset=false;
    elseif cc=='Q'
        CurrentFrame=CurrentFrame-5;
        frameReset=false;    
    % Save
    elseif cc=='s'
        save([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'],'schnitzcells');
        
    % Moving between Nuclear Cycles
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
    
    % Switch to lineaging mode
    if ct==0
        NewNuclei=ginput_red(1); % Input - requires oddball dependency
        Parent=identify(NewNuclei,CurrentFrame); % Calls function to identify current nuclei whose daughters will be selected
        CurrentFrame=ncs(ncexist(Current_Nuclear_Cycle))+1;
        cc2=''; % Stores current imput
        keepDaughtering=true;
        
        while(keepDaughtering)
            display('Select Daughter Nucleus (. to move a frame forward and , to move a frame backward)');
            display('Again, mouseclick to start selecting daughters');
            display('Click on the same nucleus twice if only child, n if no children');
            displayFrame(CurrentFrame,Parent,1); % Displays current frame with parent highlighted
            
            % Waits for user input again
            ct2=waitforbuttonpress;
            cc2=get(gcf,'currentcharacter');
            cm2=get(gcf,'CurrentPoint');
            if cc2=='.'
                CurrentFrame=CurrentFrame+1;
            elseif cc2=='w'
                CurrentFrame=CurrentFrame+1;
            elseif cc2=='W'
                CurrentFrame=CurrentFrame+5;
            elseif cc2==','
                CurrentFrame=CurrentFrame-1;
            elseif cc2=='q'
                CurrentFrame=CurrentFrame-1;
            elseif cc2=='Q'
                CurrentFrame=CurrentFrame-5;
            elseif cc2=='n'
                keepDaughtering=false;
            end
            if ct2==0
                Daughters=zeros(2,1);
                NewNuclei=zeros(2);
                NewNuclei(1,:)=ginput_red(1);
                Daughters(1)=identify(NewNuclei(1,:),CurrentFrame);
                NewNuclei(2,:)=ginput_red(1);
                Daughters(2)=identify(NewNuclei(2,:),CurrentFrame);
                if ~any(Daughters==Parent)
                    schnitzcells(Parent).E=Daughters(1);
                    schnitzcells(Daughters(1)).P=Parent;
                    if(Daughters(1)~=Daughters(2))
                        schnitzcells(Parent).D=Daughters(2);
                        schnitzcells(Daughters(2)).P=Parent;
                    end
                    schnitzcells(Daughters(1)).D=[];
                    schnitzcells(Daughters(2)).D=[];
                    schnitzcells(Daughters(1)).E=[];
                    schnitzcells(Daughters(2)).E=[];
                else
                    display('Parent Splitting');
                    FrameIndex=find(schnitzcells(Parent).frames==CurrentFrame);
                    newIndex=size(schnitzcells,1)+1;
                    schnitzcells(newIndex).frames=schnitzcells(Parent).frames(FrameIndex:end);
                    schnitzcells(newIndex).cenx=schnitzcells(Parent).cenx(FrameIndex:end);
                    schnitzcells(newIndex).ceny=schnitzcells(Parent).ceny(FrameIndex:end);
                    schnitzcells(newIndex).cellno=schnitzcells(Parent).cellno(FrameIndex:end);
                    schnitzcells(newIndex).len=schnitzcells(Parent).len(FrameIndex:end);
                    
                    schnitzcells(Parent).frames=schnitzcells(Parent).frames(1:FrameIndex-1);
                    schnitzcells(Parent).cenx=schnitzcells(Parent).cenx(1:FrameIndex-1);
                    schnitzcells(Parent).ceny=schnitzcells(Parent).ceny(1:FrameIndex-1);
                    schnitzcells(Parent).cellno=schnitzcells(Parent).cellno(1:FrameIndex-1);
                    schnitzcells(Parent).len=schnitzcells(Parent).len(1:FrameIndex-1);
                    if (Daughters(1)==Parent)
                        schnitzcells(Parent).E=newIndex;
                        schnitzcells(Parent).D=Daughters(2);
                        
                    else
                        schnitzcells(Parent).D=newIndex;
                        schnitzcells(Parent).E=Daughters(1);
                    end
                end
                keepDaughtering=false;
            end
        end
    end
    if keeptracking && frameReset
        CurrentFrame=ncs(ncexist(Current_Nuclear_Cycle))-1;
    end
    CurrentFrame=max(1,CurrentFrame);
end
close all;


%% Display
    function displayFrame(CurrentFrame,select,varargin)
        daughtering=false;
        
        % If daughtering, the display is slightly different: there are no
        % families being displayed
        if ~isempty(varargin)
            if varargin{1}==1
                daughtering=true;
            end
        end
        % Displays current frame
        %         display('Image Current Frame')
        %         display(CurrentFrame);
        ImageHis=imread([PreProcPath,filesep,FilePrefix(1:end-1),filesep,...
            FilePrefix(1:end-1),'-His_',iIndex(CurrentFrame,3),'.tif']);
        
        for j=1:size(schnitzcells,2)
            
            if any(schnitzcells(j).frames==CurrentFrame)
                index=find(schnitzcells(j).frames==CurrentFrame);
                circles=[schnitzcells(j).cenx(index) schnitzcells(j).ceny(index) 7];
                colour='w';
                if ~isempty(schnitzcells(j).D) || ~isempty(schnitzcells(j).E)
                    colour='g';
                end
                if select==j
                    colour='y';
                end
                if daughtering % All the nuclei are shown without families only if in daughtering mode
                    ImageHis=insertShape(ImageHis, 'filledcircle', circles,'Color',colour);
                else
                    % ImageHis=insertText(ImageHis, [schnitzcells(j).cenx(index) schnitzcells(j).ceny(index)], j, 'FontSize',9, 'BoxOpacity',0, 'TextColor','r');
                    
                    % If not daughtering, all the families are shown
                    for jj=1:size(families,1)
                        if any(families(jj,:)==j)
                            gen=find(families(jj,:)==j);
                            if length(gen)>1
                                error('ABORT! You made an incorrect click. For everybodys welfare, this will be aborted');
                            end
                            if gen>1 && schnitzcells(families(jj,1)).frames(end)>CurrentFrame-10
                                if jj>size(colours,1)
                                    colours(jj,:)=[rand rand rand];
                                end
                                if all(colours(jj,:)==0)
                                    colours(jj,:)=[rand rand rand];
                                end
                                ImageHis=insertShape(ImageHis, 'filledcircle', circles,'Color',65536*colours(jj,:));
                                parent=families(jj,1);
                                circles2=[schnitzcells(parent).cenx(end) schnitzcells(parent).ceny(end) 7];
                                ImageHis=insertShape(ImageHis, 'filledcircle', circles2,'Color',65536*colours(jj,:)/2);
                                ImageHis=insertShape(ImageHis,'line',[schnitzcells(parent).cenx(end),schnitzcells(parent).ceny(end),...
                                    schnitzcells(j).cenx(index) schnitzcells(j).ceny(index)]);
                            else
                                circles=[schnitzcells(j).cenx(index) schnitzcells(j).ceny(index) 12];
                                ImageHis=insertShape(ImageHis, 'circle', circles,'Color',colour);
                            end
                        end
                    end
                end
            end
        end
        
        % Shows all the nuclei from the previous frame if picking daughters
        if CurrentFrame~=1 && daughtering
            for j=1:size(schnitzcells,2)
                
                if any(schnitzcells(j).frames==CurrentFrame-1)
                    index=find(schnitzcells(j).frames==CurrentFrame-1);
                    circles=[schnitzcells(j).cenx(index) schnitzcells(j).ceny(index) 15];
                    colour=[0 65536/4 65536/4];
                    ImageHis=insertShape(ImageHis, 'circle', circles,'Color',colour);
                    
                end
            end
        end
        
        title(['Current Frame: ',CurrentFrame]);
        imshow(ImageHis,'Border','Tight')
        %    set(gcf,'units', 'normalized', 'position',[.1   .55   .4   .35])
        set(gcf,'MenuBar','none','ToolBar','none')
        
        
    end

%% Identify Clicked Nuclei
% Subfunction to identify clicked nuclei
    function ClickedSchnitz=identify(NewNuclei,CurrentFrame)
        
        
        ClickedSchnitz=[];
        if isempty(NewNuclei)
            warning('Please select the correct nucleus')
        else
            [NClicks,Dummy]=size(NewNuclei);
            
            %Find the nuclei/schnitz that we clicked on
            for i=1:NClicks
                %Find which schnitz this corresponds to by finding the
                %minimum distance among all nuclei to the point that was
                %clicked
                ClickedSchnitz=[];
                xPosSuspect=[];
                yPosSuspect=[];
                min_Distance=Inf;
                for j=1:size(schnitzcells,2)
                    if any(schnitzcells(j).frames==CurrentFrame)
                        indx=find(schnitzcells(j).frames==CurrentFrame);
                        dist=(schnitzcells(j).cenx(indx)-NewNuclei(1))^2+(schnitzcells(j).ceny(indx)-NewNuclei(2))^2;
                        if dist<min_Distance
                            min_Distance=dist;
                            ClickedSchnitz=j;
                        end
                    end
                end
            end
        end
    end
end