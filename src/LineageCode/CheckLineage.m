function  schnitzcells=CheckLineage(Prefix,FrameInfo,Ellipses)
% Finds the times of cell division (by using nuclei radius,
% and orphan the nuclei if they survive mitosis)


%% Extract all data
%Get the folders, including the default Dropbox one
[SourcePath, FISHPath, DefaultDropboxFolder, DropboxFolder, MS2CodePath, PreProcPath,...
configValues, movieDatabasePath] = DetermineAllLocalFolders(Prefix);



%Determine division times
%Load the information about the nc from moviedatabase file
[Date, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution,...
Channel1, Channel2, Objective, Power, DataFolder, DropboxFolderName, Comments,...
nc9, nc10, nc11, nc12, nc13, nc14, CF] = getExperimentDataFromMovieDatabase(Prefix, DefaultDropboxFolder)

if (~isempty(findstr(Prefix,'Bcd')))&(isempty(findstr(Prefix,'BcdE1')))&...
        (isempty(findstr(Prefix,'NoBcd')))&(isempty(findstr(Prefix,'Bcd1x')))&(isempty(findstr(Prefix,'Bcd4x')))
    warning('This step in CheckParticleTracking will most likely have to be modified to work')
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


%Convert the prefix into the string used in moviedatabase file
Dashes=findstr(Prefix,'-');

ncs=[nc9,nc10,nc11,nc12,nc13,nc14];

if (length(find(isnan(ncs)))==length(ncs))|(length(ncs)<6)
    error('Have the ncs been defined in MovieDatabase?')
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

load([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat']);
linbackup=[DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin_backup.mat'];
save(linbackup,'schnitzcells');


cpexist=true; % Does Compiled Particles Exist?
try
    load([DropboxFolder,filesep,Prefix,filesep,'Particles.mat']);
catch
    warning('No Compiled Particles found!');
    cpexist=false;
end

try
    load([DropboxFolder,filesep,FilePrefix(1:end-1),filesep,'Ellipses.mat'])
    UseHistoneOverlay=1;
catch
    UseHistoneOverlay=0;
    warning('No Ellipses found!');
end
try
    load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat']);
catch
    warning('No Frame Info found!');
end

dummy=FrameInfo(1);
dims(1)=dummy.LinesPerFrame;
dims(2)=dummy.PixelsPerLine;


%% Determine Order of Nuclei
Nuclei=zeros(6,length(schnitzcells));
count=zeros(6,1);
for i=1:length(schnitzcells)
    for ii=1:ncs
        if abs(max(schnitzcells(i).frames)-ncs(ii))<3
            Nuclei(ii,count(ii))=i;
            count(ii)=count(ii)+1;
        end
    end
end
for ii=1:ncs
    keep=Nuclei(ii,:)~=0;
    Nuclei(ii,:)=Nuclei(ii,keep);
end
%% GUI
schnitzcells=FixLineage(Prefix,FrameInfo)

keeptracking=true;
ncexist=find(ncs~=0);
Current_Nuclear_Cycle=1;

CurrentNuc=1;
CurrentFrame=ncs(ncexist(Current_Nuclear_Cycle))-1;
while keeptracking
    disp(CurrentFrame,[]);
    display(CurrentFrame);
    Current_NC=Current_Nuclear_Cycle+7+min(ncexist);
    display(Current_NC);
    display('x to save , to move frame backward . to move frame forward mouseclick to start');
    ct=waitforbuttonpress;
    cc=get(gcf,'currentcharacter');
    cm=get(gcf,'CurrentPoint');
    display('Click to start Lineaging');
    frameReset=true;
    
    if cc=='x'        
        save([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'],'schnitzcells');
        keeptracking=false;
    elseif cc=='.'
        CurrentFrame=CurrentFrame+1;
        frameReset=false;
    elseif cc=='w'
        CurrentFrame=CurrentFrame+1;
        frameReset=false;
    elseif cc==','
        CurrentFrame=CurrentFrame-1;
        frameReset=false;
    elseif cc=='q'
        CurrentFrame=CurrentFrame-1;
        frameReset=false;
    elseif cc=='s'      
        save([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'],'schnitzcells');
    elseif cc=='t'
        Current_Nuclear_Cycle=Current_Nuclear_Cycle+1;
        if Current_Nuclear_Cycle>length(ncexist)
            keeptracking=false;
        end
    end
    if ct==0
        display('Now Select Parent Nucleus');
        NewNuclei=zeros(1,1);
        NewNuclei=ginput_red(1);
        Parent=identify(NewNuclei)
        CurrentFrame=ncs(ncexist(Current_Nuclear_Cycle))+1;
        cc2='';
        keepDaughtering=true;
        display('Select Daughter Nucleus (. to move a frame forward and , to move a frame backward)');
        display('Again, mouseclick to start selecting daughters');
        display('Click on the same nucleus twice if only child, n if no children');
        while(keepDaughtering)
            disp(CurrentFrame,Parent);
            
            ct2=waitforbuttonpress;
            cc2=get(gcf,'currentcharacter');
            cm2=get(gcf,'CurrentPoint');
            if cc2=='.' 
                CurrentFrame=CurrentFrame+1;
            elseif cc2=='w'
                CurrentFrame=CurrentFrame+1;
            elseif cc2==',' 
                CurrentFrame=CurrentFrame-1;
            elseif cc2=='q'
                CurrentFrame=CurrentFrame-1;
            elseif cc2=='n'
                keepDaughtering=false;
            end
            if ct2==0
                NewNuclei=ginput_red(2);
                
                Daughters=identify(NewNuclei);
                if ~any(Daughters==Parent)
                    schnitzcells(Parent).E=Daughters(1);
                    schnitzcells(Daughters(1)).P=Parent
                    if(Daughters(1)~=Daughters(2))
                        schnitzcells(Parent).D=Daughters(2);
                        schnitzcells(Daughters(2)).P=Parent
                    end
                    schnitzcells(Daughters(1)).D=[];
                    schnitzcells(Daughters(2)).D=[];
                    schnitzcells(Daughters(1)).E=[];
                    schnitzcells(Daughters(2)).E=[];
                elseif (Daughters(1)==Parent)
                    if ~isempty(schnitzcells(Daughters(1)).D)
                        schnitzcells(Parent).E=schnitzcells(Daughters(1)).D;
                    elseif ~isempty(schnitzcells(Daughters(1)).E)
                        schnitzcells(Parent).E=schnitzcells(Daughters(1)).E;
                    end
                    schnitzcells(Parent).D=Daughters(2);
                else
                    if ~isempty(schnitzcells(Daughters(2)).D)
                        schnitzcells(Parent).E=schnitzcells(Daughters(2)).D;
                    elseif ~isempty(schnitzcells(Daughters(2)).E)
                        schnitzcells(Parent).E=schnitzcells(Daughters(2)).E;
                    end
                    schnitzcells(Parent).D=Daughters(1);
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
    function disp(CurrentFrame,select)
        close all;
        % Displays current frame
        if UseHistoneOverlay
            figure;
            
            try
                ImageHis=imread([PreProcPath,filesep,FilePrefix(1:end-1),filesep,...
                    FilePrefix(1:end-1),'-His_',iIndex(CurrentFrame,3),'.tif']);
            catch %Had to do this for KITP
                ImageHis=imread([PreProcPath,filesep,FilePrefix(1:end-1),filesep,...
                    FilePrefix(1:end-1),'_His_',iIndex(CurrentFrame,3),'.tif']);
            end
        else
            % ImageHis=ones(512,512);
            ImageHis=ones(dims(1),dims(2));
            x=Ellipses{CurrentFrame,1}(:,1); % x distance array for the frame
            y=Ellipses{CurrentFrame,1}(:,2); % y distance array for the frame
            
            % Choose a nucleus
            for j=1:length(y)
                I=imread(strcat(PreProcPath,filesept,Prefix,filesep,Prefix,'-His_',num2str(CurrentFrame),'.tif'));
                circles=[Ellipses{CurrentFrame,1}(j,1) Ellipses{CurrentFrame,1}(j,2) 15];
                I=insertShape(I, 'circle', circles,'Color','w');
            end
            for j=1:length(Nuclei)
                frame=find(schnitzcells(1,Nuclei(j)).frames==CurrentFrame);
                if(isempty(frame))
                    frame=length(schnitzcells(1,Nuclei(j)).frames);
                    %display('Couldnt find frame');
                    failed=failed+1;
                end
                circles=[schnitzcells(1,Nuclei(j)).cenx(frame) schnitzcells(1,1,Nuclei(j)).ceny(frame) 15];
                I=insertShape(I, 'circle', circles,'Color','g');
                
            end
            ImageHis=I;
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
            end
        end
        


%         failed=0;
%         for j=1:size(Ellipses{CurrentFrame},1)
%             I=imread(strcat(PreProcPath,filesept,Prefix,filesep,Prefix,'-His_',num2str(i),'.tif'));
%             circles=[Ellipses{CurrentFrame,1}(j,1) Ellipses{CurrentFrame,1}(j,2) 15];
%             colour='w';
%             if 
%             I=insertShape(I, 'circle', circles,'Color',colour);
%         end
%             for j=1:length(Nuclei)
%                 frame=find(schnitzcells(1,Nuclei(j)).frames==CurrentFrame);
%                 if(isempty(frame))
%                     frame=length(schnitzcells(1,Nuclei(j)).frames);
%                     %display('Couldnt find frame');
%                     failed=failed+1;
%                 end
%                 circles=[schnitzcells(1,Nuclei(j)).cenx(frame) schnitzcells(1,1,Nuclei(j)).ceny(frame) 15];
%                 I=insertShape(I, 'circle', circles,'Color','g');
%                 
%             end
        imshow(ImageHis,'Border','Tight')
        %    set(gcf,'units', 'normalized', 'position',[.1   .55   .4   .35])
        set(gcf,'MenuBar','none','ToolBar','none')
        
        
        
        %         if ZoomMode
        %             %Find the closest frame
        %             [Dummy,MinIndex]=min((Particles{CurrentChannel}(CurrentParticle).Frame-CurrentFrame).^2);
        %             if length(MinIndex)>1
        %                 MinIndex=MinIndex(1);
        %             end
        %             [xForZoom,yForZoom]=fad2xyzFit(Particles{CurrentChannel}(CurrentParticle).Frame(MinIndex),fad(CurrentChannel), 'addMargin');
        %             xForZoom=xForZoom(Particles{CurrentChannel}(CurrentParticle).Index(MinIndex));
        %             yForZoom=yForZoom(Particles{CurrentChannel}(CurrentParticle).Index(MinIndex));
        %
        %
        %             xlim([xForZoom-ZoomRange,xForZoom+ZoomRange])
        %             ylim([yForZoom-ZoomRange/2,yForZoom+ZoomRange/2])
        %         end
    end

%% Identify Clicked Nuclei
    function ClickedSchnitz=identify(NewNuclei)
        
        
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
        end
    end
end