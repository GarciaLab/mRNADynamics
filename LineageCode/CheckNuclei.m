function  schnitzcells=CheckNuclei(Prefix,varargin)
% GUI to select the daughter nuclei and manually curate the lineage
SpeedMode=false;
if strcmpi(varargin, 'speedmode')
    SpeedMode = true;
end

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
if (~isempty(findstr(Prefix,'Bcd')))&&(isempty(findstr(Prefix,'BcdE1')))&&...
        (isempty(findstr(Prefix,'NoBcd')))&&(isempty(findstr(Prefix,'Bcd1x')))&&...
        (isempty(findstr(Prefix,'Bcd4x')))
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
pexist=true; % Does Particles.mat Exist?
try
    dummy1=load([DropboxFolder,filesep,Prefix,filesep,'Particles.mat']);
    Particles=dummy1.Particles;
catch
    warning('No Particles.mat found. No changes will be made to the Particles structure.');
    pexist=false;
end
CompiledParticles=struct;
cpexist=true; % Does Compiled Particles Exist?
try
    dummy2=load([DropboxFolder,filesep,Prefix,filesep,'CompiledParticles.mat']);
    CompiledParticles=dummy2.CompiledParticles;
catch
    warning('No CompiledParticles.mat found. No changes will be made to the CompiledParticles structure.');
    cpexist=false;
end

FrameInfo=struct;
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
for n=1:len
    if (schnitzcells(1,n).D>0)
        count=count+1;
        if ~isempty(schnitzcells(1,n).E) && (schnitzcells(1,n).E~=0)
            families(count,:)=[n schnitzcells(1,n).E schnitzcells(1,n).D]';
        else
            families(count,:)=[n 0 schnitzcells(1,n).D]';
        end
    elseif (schnitzcells(1,n).E>0)
        count=count+1;
        families(count,:)=[n schnitzcells(1,n).E 0]';
    end
end
%% GUI

keeptracking=true;
ncexist=find(ncs~=0 & ~isnan(ncs)); % Logical array to tell if the Nuclear cycle exists
Current_Nuclear_Cycle=1;

% The user can choose to cycle between Nuclear Cycles or Frames, so these
% variables keep track of the Current Nuclear Cycle and the Current Frame
CurrentFrame=ncs(ncexist(Current_Nuclear_Cycle))-1;

% This is an array that stores a random colour for each family
colours=zeros(length(families),3);
HistoneMode=false;
while keeptracking
    %% Fix Corner Cases
    for n=1:len
        if schnitzcells(n).P==0
            schnitzcells(n).P=[];
        end
        if schnitzcells(n).D==0
            schnitzcells(n).D=[];
        end
        if schnitzcells(n).E==0
            schnitzcells(n).E=[];
        end
        Parent=schnitzcells(n).P;
        ChildD=schnitzcells(n).D;
        ChildE=schnitzcells(n).E;
        % Make sure all families are unique and two-sided
        if ~isempty(ChildE)
            if schnitzcells(ChildE).P~=n
                schnitzcells(n).E=[];
            end
        end
        if ~isempty(ChildD)
            if schnitzcells(ChildD).P~=n
                schnitzcells(n).D=[];
            end
        end
        if ~isempty(Parent)
            if ~isempty(schnitzcells(Parent).D) && ~isempty(schnitzcells(Parent).E)
                if schnitzcells(Parent).D~=n & schnitzcells(Parent).E~=n
                    schnitzcells(n).P=[];
                end
            end
        end
        
        % Daughters are older than parents
        if ~isempty(ChildD)
            if schnitzcells(ChildD).frames(1)<schnitzcells(n).frames(end)
                schnitzcells(n).D=[];
                schnitzcells(ChildD).P=[];
            end
        end
        if ~isempty(ChildE)
            if schnitzcells(ChildE).frames(1)<schnitzcells(n).frames(end)
                schnitzcells(n).E=[];
                schnitzcells(ChildE).P=[];
            end
        end
        
        % Sister is a mother
        
        % Two families in same frame
    end
    %% Family Checking
    % The code below finds the families once again. This is inefficient,
    % consider changing
    len=size(schnitzcells,2);
    count=0;
    families=[];
    for n=1:len
        if (schnitzcells(1,n).D>0)
            count=count+1;
            if ~isempty(schnitzcells(1,n).E) && (schnitzcells(1,n).E~=0)
                families(count,:)=[n schnitzcells(1,n).E schnitzcells(1,n).D]';
            else
                families(count,:)=[n 0 schnitzcells(1,n).D]';
            end
        elseif (schnitzcells(1,n).E>0)
            count=count+1;
            families(count,:)=[n schnitzcells(1,n).E 0]';
        end
    end
    
    % displayFrame is the function that generates the display
    if HistoneMode
        displayFrame(CurrentFrame,[],2);
    else
        displayFrame(CurrentFrame,[]);
    end
    
    % Current_NC is only for display. It will not be used computationally
    [~,Current_Nuclear_Cycle]=min(abs(CurrentFrame-ncs(ncexist)));
    Current_NC=Current_Nuclear_Cycle+7+min(ncexist);
    display(Current_NC);
    
    % Instructions
    if ~SpeedMode
    display('x to save , to move frame backward . to move frame forward mouseclick to start');
    end
    % This stores the variable to see if the current frame is the default
    % frame - it is used to toggle whether the user is moving between
    % frames
    frameReset=true;
    
    % cc, ct and cm gather the current user input (cm is not used yet)
    if ~SpeedMode
    display('Click to start Lineaging (Mitosis Selection). ');
    end
    ct=waitforbuttonpress;
    cc=get(gcf,'currentcharacter');
    cm=get(gcf,'CurrentPoint');
    
    
    % Delete mode
    %     if cc=='d'
    %
    %         cc3='';
    %         while strcmp(cc3,'n')==0
    %             display('EnteringDeleteMode');
    %             ct3=waitforbuttonpress;
    %             cc3=get(gcf,'currentcharacter');
    %             cm3=get(gcf,'CurrentPoint');
    %             toDelete=identify(cm3,CurrentFrame);
    %             childD=schnitzcells(toDelete).D;
    %             childE=schnitzcells(toDelete).E;
    %             paren=schnitzcells(toDelete).P;
    %             if ~isempty(childD)
    %                 schnitzcells(childD).P=[];
    %             end
    %             if ~isempty(childE)
    %                 schnitzcells(childE).P=[];
    %             end
    %             if ~isempty(paren)
    %                 if schnitzcells(paren).D==toDelete
    %                     schnitzcells(paren).D=[];
    %                 elseif schnitzcells(paren).D==toDelete
    %                     schnitzcells(paren).D=[];
    %                 end
    %             end
    %             schnitzcells(toDelete)=[];
    %         end
    %     end
    
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
    elseif cc=='h'
        HistoneMode=~HistoneMode;
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
    %% Identity Mode
    elseif cc=='i'
       % display('Reached0');
        cc7='a';
        while cc7~='i'
            %display('Reached1');
            displayFrame(CurrentFrame,[],3);
            ct7=waitforbuttonpress;
            cc7=get(gcf,'currentcharacter');
            if cc7=='.'
                CurrentFrame=CurrentFrame+1;
            elseif cc7=='w'
                CurrentFrame=CurrentFrame+1;
            elseif cc7=='W'
                CurrentFrame=CurrentFrame+5;
            elseif cc7=='>'
                CurrentFrame=CurrentFrame+5;
            elseif cc7==','
                CurrentFrame=CurrentFrame-1;
            elseif cc7=='q'
                CurrentFrame=CurrentFrame-1;
            elseif cc7=='Q'
                CurrentFrame=CurrentFrame-5;
            elseif cc7=='j'
                jump=input('Which Frame do you want to jump to?');
                CurrentFrame=jump;
                frameReset=false;
            elseif cc7=='<'
                CurrentFrame=CurrentFrame-5;
            end
            if ct7==0
                NewNuclei=ginput_red(1); % Input - requires oddball dependency
                Parent=identify(NewNuclei,CurrentFrame); % Calls function to identify current nuclei whose daughters will be selected
                
                frameReset=false;
                num='a';
                while ~isnumeric(num)
                    num=input('Input the number     ');
                end
                
                % Fix the current schnitz
                display('Current schnitz fixed');
                frameIndex=find(schnitzcells(Parent).frames==CurrentFrame);
                posX=schnitzcells(Parent).cenx(frameIndex);
                posY=schnitzcells(Parent).ceny(frameIndex);
                cellF=schnitzcells(Parent).cellno(frameIndex);
                lenF=schnitzcells(Parent).len(frameIndex);
                if frameIndex~=1
                    beforeFrames=schnitzcells(Parent).frames(1:frameIndex-1);
                    beforeCenx=schnitzcells(Parent).cenx(1:frameIndex-1);
                    beforeCeny=schnitzcells(Parent).ceny(1:frameIndex-1);
                    beforeCell=schnitzcells(Parent).cellno(1:frameIndex-1);
                    beforeLen=schnitzcells(Parent).len(1:frameIndex-1);
                    
                    % and the others
                else
                    beforeFrames=[];
                    beforeCenx=[];
                    beforeCeny=[];
                    beforeCell=[];
                    beforeLen=[];
                    % and the others
                end
                if frameIndex~=length(schnitzcells(Parent).frames)
                    afterFrames=schnitzcells(Parent).frames(frameIndex+1:end);
                    afterCenx=schnitzcells(Parent).cenx(frameIndex+1:end);
                    afterCeny=schnitzcells(Parent).ceny(frameIndex+1:end);
                    afterCell=schnitzcells(Parent).cellno(frameIndex+1:end);
                    afterLen=schnitzcells(Parent).len(frameIndex+1:end);
                else
                    afterFrames=[];
                    afterCenx=[];
                    afterCeny=[];
                    afterCell=[];
                    afterLen=[];
                end
                if length(schnitzcells(num).frames)==1
                    beforeFrames=Inf;
                end
                schnitzcells(Parent).frames=[beforeFrames; afterFrames];
                schnitzcells(Parent).cenx=[beforeCenx afterCenx];
                schnitzcells(Parent).ceny=[beforeCeny afterCeny];
                schnitzcells(Parent).cellno=[beforeCell afterCell];
                schnitzcells(Parent).len=[beforeLen afterLen];
                
                % If match in frame, put conflict at end of structure
                if any(schnitzcells(num).frames==CurrentFrame)
                    %                 display('Conflict');
                    frameIndex2=find(schnitzcells(num).frames==CurrentFrame);
                    newSpot=size(schnitzcells,2)+1;
                    schnitzcells(newSpot).frames=schnitzcells(num).frames(frameIndex2);
                    schnitzcells(newSpot).ceny=schnitzcells(num).ceny(frameIndex2);
                    schnitzcells(newSpot).cenx=schnitzcells(num).cenx(frameIndex2);
                    schnitzcells(newSpot).cellno=schnitzcells(num).cellno(frameIndex2);
                    schnitzcells(newSpot).len=schnitzcells(num).len(frameIndex2);
                    
                    schnitzcells(num).frames(frameIndex2)=CurrentFrame;
                    schnitzcells(num).ceny(frameIndex2)=posX;
                    schnitzcells(num).cenx(frameIndex2)=posY;
                    schnitzcells(num).cellno(frameIndex2)=cellF;
                    schnitzcells(num).len(frameIndex2)=lenF;
                else
                    %                 display('Attachment done');
                    % Otherwise attach to target
                    frameIndex=1;
                    for yy=1:length(schnitzcells(num).frames)
                        if schnitzcells(num).frames(yy)>CurrentFrame
                            frameIndex=yy;
                        end
                    end
                    if frameIndex~=1
                        beforeFrames=schnitzcells(num).frames(1:frameIndex-1);
                        beforeCenx=schnitzcells(num).cenx(1:frameIndex-1);
                        beforeCeny=schnitzcells(num).ceny(1:frameIndex-1);
                        beforeCell=schnitzcells(num).cellno(1:frameIndex-1);
                        beforeLen=schnitzcells(num).len(1:frameIndex-1);
                        
                        % and the others
                    else
                        beforeFrames=[];
                        beforeCenx=[];
                        beforeCeny=[];
                        beforeCell=[];
                        beforeLen=[];
                        % and the others
                    end
                    if frameIndex~=length(schnitzcells(num).frames)
                        afterFrames=schnitzcells(num).frames(frameIndex+1:end);
                        afterCenx=schnitzcells(num).cenx(frameIndex+1:end);
                        afterCeny=schnitzcells(num).ceny(frameIndex+1:end);
                        afterCell=schnitzcells(num).cellno(frameIndex+1:end);
                        afterLen=schnitzcells(num).len(frameIndex+1:end);
                    else
                        afterFrames=[];
                        afterCenx=[];
                        afterCeny=[];
                        afterCell=[];
                        afterLen=[];
                    end
                    if length(schnitzcells(num).frames)==1
                        beforeFrames=Inf;
                    end
                    schnitzcells(num).frames=[beforeFrames; CurrentFrame; afterFrames];
                    schnitzcells(num).cenx=[beforeCenx posX afterCenx];
                    schnitzcells(num).ceny=[beforeCeny posY afterCeny];
                    schnitzcells(num).cellno=[beforeCell cellF afterCell];
                    schnitzcells(num).len=[beforeLen lenF afterLen];
                end
            end
        end
    end
    
    % Switch to lineaging mode
    if ct==0
        NewNuclei=ginput_red(1); % Input - requires oddball dependency
        Parent=identify(NewNuclei,CurrentFrame); % Calls function to identify current nuclei whose daughters will be selected
        CurrentFrame=max(CurrentFrame,ncs(ncexist(Current_Nuclear_Cycle))+1);
        cc2=''; % Stores current imput
        keepDaughtering=true;
        
        while(keepDaughtering)
            if ~SpeedMode
            display('Select Daughter Nucleus (. to move a frame forward and , to move a frame backward)');
            display('Again, mouseclick to start selecting daughters');
            display('Click on the same nucleus twice if only child, n if no children');
            end
            if HistoneMode
                displayFrame(CurrentFrame,Parent,2);
            else
                displayFrame(CurrentFrame,Parent,1); % Displays current frame with parent highlighted
            end
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
            elseif cc2=='>'
                CurrentFrame=CurrentFrame+5;
            elseif cc2==','
                CurrentFrame=CurrentFrame-1;
            elseif cc2=='q'
                CurrentFrame=CurrentFrame-1;
            elseif cc2=='Q'
                CurrentFrame=CurrentFrame-5;
            elseif cc2=='j'
                jump=input('Which Frame do you want to jump to?');
                CurrentFrame=jump;
                frameReset=false;
            elseif cc2=='<'
                CurrentFrame=CurrentFrame-5;
            elseif cc2=='n'
                keepDaughtering=false;
            elseif cc2=='h'
                HistoneMode=~HistoneMode;
                
            elseif cc2=='d'
                keepDaughtering=false;
                if ~isempty(schnitzcells(Parent).E)
                    schnitzcells(schnitzcells(Parent).E).P=[];
                end
                if ~isempty(schnitzcells(Parent).D)
                    schnitzcells(schnitzcells(Parent).D).P=[];
                end
                schnitzcells(Parent).E=[];
                schnitzcells(Parent).D=[];
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
                    if FrameIndex==1
                        error('Cannot divide in the first frame that it exists');
                    end
                    newIndex=size(schnitzcells,2)+1;
                    display(Parent);
                    display(FrameIndex);
                    schnitzcells(newIndex).frames=schnitzcells(Parent).frames(FrameIndex:end);
                    schnitzcells(newIndex).cenx=schnitzcells(Parent).cenx(FrameIndex:end);
                    schnitzcells(newIndex).ceny=schnitzcells(Parent).ceny(FrameIndex:end);
                    schnitzcells(newIndex).cellno=schnitzcells(Parent).cellno(FrameIndex:end);
                    schnitzcells(newIndex).len=schnitzcells(Parent).len(FrameIndex:end);
                    schnitzcells(newIndex).P=Parent;
                    schnitzcells(newIndex).D=[];
                    schnitzcells(newIndex).E=[];
                    
                    schnitzcells(Parent).frames=schnitzcells(Parent).frames(1:FrameIndex-1);
                    schnitzcells(Parent).cenx=schnitzcells(Parent).cenx(1:FrameIndex-1);
                    schnitzcells(Parent).ceny=schnitzcells(Parent).ceny(1:FrameIndex-1);
                    schnitzcells(Parent).cellno=schnitzcells(Parent).cellno(1:FrameIndex-1);
                    schnitzcells(Parent).len=schnitzcells(Parent).len(1:FrameIndex-1);
                    if (Daughters(1)==Parent)
                        schnitzcells(Parent).E=newIndex;
                        schnitzcells(Parent).D=Daughters(2);
                        schnitzcells(Daughters(2)).P=Parent;
                    else
                        schnitzcells(Parent).D=newIndex;
                        schnitzcells(Parent).E=Daughters(1);
                        schnitzcells(Daughters(1)).P=Parent;
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
if pexist
    if cpexist
        for sc=1:size(schnitzcells,2)
            schnitzcells(sc).AssociatedParticle=[];
        end
        for parts=1:size(CompiledParticles,2)
            if CompiledParticles(parts).Approved~=-1
                NumFrames=length(CompiledParticles(parts).Frame);
                Mid=ceil(NumFrames/2);
                FrameNo=CompiledParticles(parts).Frame(Mid);
                Nuc=identify([CompiledParticles(parts).xPos(Mid),CompiledParticles(parts).yPos(Mid)],FrameNo);
                CompiledParticles(parts).Nucleus=Nuc;
                try
                    schnitzcells(Nuc).AssociatedParticle=[schnitzcells(Nuc).AssociatedParticle,parts];
                catch
                    schnitzcells(Nuc).AssociatedParticle=parts;
                end
            end
        end
        save([DropboxFolder,filesep,Prefix,filesep,'CompiledParticles.mat'],'CompiledParticles');
    end
    
    for sc=1:size(schnitzcells,2)
        schnitzcells(sc).AssociatedParticle=[];
    end
    fad=dummy1.fad;
    fad2=dummy1.fad2;
    Threshold1=dummy1.Threshold1;
    Threshold2=dummy1.Threshold2;
    if ~isfield(Particles,'xPos')
        for i=1:length(Particles)
            for j=1:length(Particles(i).Frame)
                [x,y]=fad2xyzFit(Particles(i).Frame(j),fad, 'addMargin');
                Particles(i).xPos(j)=x(Particles(i).Index(j));
                Particles(i).yPos(j)=y(Particles(i).Index(j));
            end
        end
    end
    for parts=1:size(Particles,2)
        if Particles(parts).Approved~=-1
            NumFrames=length(Particles(parts).Frame);
            Mid=ceil(NumFrames/2);
            FrameNo=Particles(parts).Frame(Mid);
            Nuc=identify([Particles(parts).xPos(Mid),Particles(parts).yPos(Mid)],FrameNo);
            Particles(parts).Nucleus=Nuc;
            try
                schnitzcells(Nuc).AssociatedParticle=[schnitzcells(Nuc).AssociatedParticle,parts];
            catch
                schnitzcells(Nuc).AssociatedParticle=parts;
            end
        end
    end
    save([DataFolder,filesep,'Particles.mat'],'Particles','fad','fad2','Threshold1','Threshold2')
    save([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'],'schnitzcells');
end

%% Display
    function displayFrame(CurrentFrame,select,varargin)
        daughtering=false;
        
        % If daughtering, the display is slightly different: there are no
        % families being displayed
        Histone=false;
        Identity=false;
        if ~isempty(varargin)
            if varargin{1}==1
                daughtering=true;
            end
            if varargin{1}==2
                Histone=true;
            end
            if varargin{1}==3
                Identity=true;
            end
        end
        % Displays current frame
        %         display('Image Current Frame')
        %         display(CurrentFrame);
        ImageHis=imread([PreProcPath,filesep,FilePrefix(1:end-1),filesep,...
            FilePrefix(1:end-1),'-His_',iIndex(CurrentFrame,3),'.tif']);
        if Identity
            if CurrentFrame~=1
                for j=1:size(schnitzcells,2)
                    
                    if any(schnitzcells(j).frames==CurrentFrame-1)
                        index=find(schnitzcells(j).frames==CurrentFrame-1);
                        circles=[schnitzcells(j).cenx(index) schnitzcells(j).ceny(index) 15];
                        colour=[0 65536/4 65536/4];
                        ImageHis=insertShape(ImageHis, 'circle', circles,'Color',colour);
                        
                    end
                end
            end
            for j=1:size(schnitzcells,2)
                
                if any(schnitzcells(j).frames==CurrentFrame)
                    index=find(schnitzcells(j).frames==CurrentFrame);
                    circles=[schnitzcells(j).cenx(index) schnitzcells(j).ceny(index) 15];
                    colour='w';
                    if ~isempty(schnitzcells(j).D) || ~isempty(schnitzcells(j).E)
                        colour='g';
                    end
                    if select==j
                        colour='y';
                    end
                    ImageHis=insertShape(ImageHis, 'circle', circles,'Color',colour);
                    ImageHis=insertText(ImageHis, [schnitzcells(j).cenx(index) schnitzcells(j).ceny(index)], j, 'FontSize',9, 'BoxOpacity',0, 'TextColor','r');
                    shifts=0;
                    for jj=1:size(families,1)
                        if any(families(jj,:)==j)
                            gen=find(families(jj,:)==j);
                            if length(gen)>1
                                warning('Major Incest');
                            end
                            if gen==1
                                col='y';
                            else
                                col='w';
                            end
                            ImageHis=insertText(ImageHis, ...
                                [schnitzcells(j).cenx(index)-15+15*shifts schnitzcells(j).ceny(index)-15], jj,...
                                'FontSize',9, 'BoxOpacity',0, 'TextColor',col);
                            shifts=shifts+1;
                        end
                    end
                end
            end
        elseif ~Histone
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
                        PaintCircles=true;
                        for jj=1:size(families,1)
                            if any(families(jj,:)==j)
                                gen=find(families(jj,:)==j);
                                if length(gen)>1
                                    error('ABORT! You made an incorrect click. For everybodys welfare, this will be aborted');
                                end
                                % Nucleus belonging to the next nuclear cycle
                                if gen>1 && schnitzcells(families(jj,1)).frames(end)>CurrentFrame-10
                                    if jj>size(colours,1)
                                        colours(jj,:)=[rand rand rand];
                                    end
                                    if all(colours(jj,:)==0)
                                        colours(jj,:)=[rand rand rand];
                                    end
                                    ImageHis=insertShape(ImageHis, 'filledcircle', circles,'Color',65536*colours(jj,:),'Opacity',0.2);
                                    parent=families(jj,1);
                                    circles2=[schnitzcells(parent).cenx(end) schnitzcells(parent).ceny(end) 7];
                                    ImageHis=insertShape(ImageHis, 'filledcircle', circles2,'Color',65536*colours(jj,:)/2,'Opacity',0.2);
                                    ImageHis=insertShape(ImageHis,'line',[schnitzcells(parent).cenx(end),schnitzcells(parent).ceny(end),...
                                        schnitzcells(j).cenx(index) schnitzcells(j).ceny(index)]);
                                    continue
                                    % Nucleus belonging to previous cycle
                                    
                                    PaintCircles=false;
                                    
                                end
                                
                            end
                        end
                        if PaintCircles
                            circles=[schnitzcells(j).cenx(index) schnitzcells(j).ceny(index) 12];
                            ImageHis=insertShape(ImageHis, 'circle', circles,'Color',colour);
                        end
                    end
                end
            end
            
            % Shows all the nuclei from the previous frame if picking daughters
            if CurrentFrame~=1 && daughtering && ~SpeedMode
                for j=1:size(schnitzcells,2)
                    
                    if any(schnitzcells(j).frames==CurrentFrame-1)
                        index=find(schnitzcells(j).frames==CurrentFrame-1);
                        circles=[schnitzcells(j).cenx(index) schnitzcells(j).ceny(index) 15];
                        colour=[0 65536/4 65536/4];
                        ImageHis=insertShape(ImageHis, 'circle', circles,'Color',colour);
                    end
                end
            end
        end
        
        
            
        
        set(gcf,'Name',['Current Frame: ',num2str(CurrentFrame)]);
        imshow(ImageHis*2,'Border','Tight')
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
            for nn=1:NClicks
                %Find which schnitz this corresponds to by finding the
                %minimum distance among all nuclei to the point that was
                %clicked
                ClickedSchnitz=[];
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