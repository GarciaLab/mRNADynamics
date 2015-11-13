function LineageFixer(varargin)
%Author: Joe Robinson
%Description: Run this code after TrackNuclei and before TrackmRNADynamics
%This is a stopgap code designed to make certain that no
%nuclei persist between nuclear cycles. It finds nuclei which pass through
%the cycle, and then splits them into a mother and a single daughter at
%anaphase
%The primary utility of this code is for when you need to be certain your
%particles are from only one nuclear cycle, but the lineages themselves are
%not very important.
%Hopefully, Hernan will fix track nuclei so this code is irrelevant. 

%Args
%Prefix=The file prefix
%Window= The window around anaphase. Default is 4

%Did you run this without parameters?
if isEmpty(varargin)
    error('Forgot Parameter. Run LineageFixer(Prefix) or LineageFixer(Prefix,window)')
end
%Prefix is first parameter
Prefix=varargin{1};
%Optional Window
if length(varargin)>1
    window=varargin{2};
else
    window=4;
end

%Get the folders, including the default Dropbox one
[SourcePath,FISHPath,DefaultDropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders;
%Now get the actual DropboxFolder
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders(Prefix);

%Get the nuclear cycle data
[Num,Txt,XLSRaw]=xlsread([DefaultDropboxFolder,filesep,'MovieDatabase.xlsx']);
XLSHeaders=Txt(1,:);
Txt=Txt(2:end,:);

ExperimentTypeColumn=find(strcmp(XLSRaw(1,:),'ExperimentType'));
ExperimentAxisColumn=find(strcmp(XLSRaw(1,:),'ExperimentAxis'));
Channel1Column=find(strcmp(XLSRaw(1,:),'Channel1'));
Channel2Column=find(strcmp(XLSRaw(1,:),'Channel2'));

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
Channel1=XLSRaw(PrefixRow,Channel1Column);
Channel2=XLSRaw(PrefixRow,Channel2Column);

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
XLSEntry=find(strcmp(XLSRaw(:,DataFolderColumn),...
    [Prefix(1:Dashes(3)-1),'\',Prefix(Dashes(3)+1:end)]));

if isempty(XLSEntry)
    XLSEntry=find(strcmp(XLSRaw(:,DataFolderColumn),...
        [Prefix(1:Dashes(3)-1),'/',Prefix(Dashes(3)+1:end)]));
    if isempty(XLSEntry)
        disp('%%%%%%%%%%%%%%%%%%%%%')
        error('Dateset could not be found. Check MovieDatabase.xlsx')
        disp('%%%%%%%%%%%%%%%%%%%%%')
    end
end


%I like this data structure, so I modified this section heavily from the
%original nuclear cycle retrieving code

nuclearCycles={0 0}; %Initializing, removed later

nuclearCycles=[nuclearCycles;{9 XLSRaw{XLSEntry,nc9Column}}];
nuclearCycles=[nuclearCycles;{10 XLSRaw{XLSEntry,nc10Column}}];
nuclearCycles=[nuclearCycles;{11 XLSRaw{XLSEntry,nc11Column}}];
nuclearCycles=[nuclearCycles;{12 XLSRaw{XLSEntry,nc12Column}}];
nuclearCycles=[nuclearCycles;{13 XLSRaw{XLSEntry,nc13Column}}];
nuclearCycles=[nuclearCycles;{14 XLSRaw{XLSEntry,nc14Column}}];

%this trims the data structure so it only has the relevant nuclear
%cycles
firstCycle=0;
for i=2:length(nuclearCycles)
    if(nuclearCycles{i,2}==0)
        firstCycle=firstCycle+1;
    end
end
nuclearCycles(1:firstCycle,:)=[];

%Load in the schnitzcells
load([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'])

%Because this is a hack code, I save the original data so that it
%doesn't get screwed up
save([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin1.mat'],'schnitzcells')

%This structure holds everything about what is being fixed
%nc= nuclear cycle
%start=The nuclear cycle start
%nuclei=indicies of nuclei which start in this cycle
%toFix=indicies of nuclei which start in this cycle and end out of it
ncData= struct('nc','start','nuclei','toFix');
%This loop makes it continue fixing until there is nothing left to fix
changes=1;
while changes>0
    changes=0;
    %Initialize ncData
    for i=1:length(nuclearCycles)
        ncData(i).nc=nuclearCycles(i,1);
        ncData(i).start=nuclearCycles(i,2);
        ncData(i).nuclei=[];
        ncData(i).toFix=[];
    end
    %Cycle through nuclei
    for j=1:length(schnitzcells)
        %done=while loop control boolean
        done=0;
        %nci=nuclear cycle index. iterates through ncData
        nci=1;
        while ~done
            if(nci<length(ncData))
                %cEnd=end of this cycle
                cEnd=cell2mat(ncData(nci+1).start);
                %If the first frame for this nucleus is before the end of
                %this cycle...
                if min(schnitzcells(j).frames)<cEnd-window
                    %Add to the nuclei for this cycle
                    ncData(nci).nuclei(end+1)=(j);
                    %If the last frame of this nucleus is after the end of
                    %this cycle...
                    if(max(schnitzcells(j).frames)>cEnd+window)
                        %Add to ToFix for this cycle
                        ncData(nci).toFix(end+1)=(j);
                    end
                    %done, get the next nucleus
                    done=1;
                end
                %otherwise check if it starts in the next nuclear cycle
                nci=nci+1;
            else
                %This is in the last nuclear cycle
                ncData(nci).nuclei(end+1)=(j);
                done=1;
            end
        end
    end
    %Iterate through the nuclear cycles
    for nci=1:length(ncData)-1
        %pull out the toFix array
        fix=ncData(nci).toFix;
        %pull out the anaphase between this cycle and the next
        anaphase=cell2mat(ncData(nci+1).start);
        %Cycle through the nuclei to fix
        for n=1:length(fix)
            %Get the relevant data for the nucleus
            nucID=fix(n);
            nuc=schnitzcells(nucID);
            nFrames=nuc.frames;
            %Gives you the frame index for anaphase
            anaFrame=anaphase-nFrames(1);
            deltas=[];
            %pulls out the distance the nucleus traveled between frames for the initial cycle
            for j=1:anaFrame-window;
                points=[nuc.cenx(j),nuc.ceny(j);nuc.cenx(j+1),nuc.ceny(j+1)];
                deltas(end+1)=pdist(points);
            end
            %get the average and standard deviation of the distances
            %travelled 
            average=mean(deltas);
            deviation=std(deltas);
            %Try to split the nucleus in a sensible fashion
            splitFrame=0;
            frame=anaFrame-window;
            done=0;
            while ~done
                if(frame>anaFrame+window)
                    done=1;
                end
                %If the distance traveled between these frames is greater
                %than one sigma from the average before anaphase, this is
                %the split
                points=[nuc.cenx(frame),nuc.ceny(frame);nuc.cenx(frame+1),nuc.ceny(frame+1)];
                if pdist(points)>average+deviation
                    splitFrame=frame+1;
                    done=1;
                end
                frame=frame+1;
            end
            %If we have a split frame
            if splitFrame>0
                nucA=schnitzcells(nucID);
                nucB=schnitzcells(nucID);
                nucA.frames= nucA.frames(1:splitFrame-1);
                nucB.frames= nucB.frames(splitFrame:end);
                
                nucA.cenx= nucA.cenx(1:splitFrame-1);
                nucB.cenx= nucB.cenx(splitFrame:end);
                
                nucA.ceny= nucA.ceny(1:splitFrame-1);
                nucB.ceny= nucB.ceny(splitFrame:end);
                
                nucA.len= nucA.len(1:splitFrame-1);
                nucB.len= nucB.len(splitFrame:end);
                
                nucA.cellno= nucA.cellno(1:splitFrame-1);
                nucB.cellno= nucB.cellno(splitFrame:end);
                %might need a test case to handle daughters that would only
                %be relevant if it "buds" so that it is maintained and
                %has a daughter
                %if(isEmpty(nucA.E))
                
                nucA.E=length(schnitzcells)+1;
                nucB.P=nucID;
                
                %end
                
                schnitzcells(end+1)=(nucB);
                schnitzcells(nucID)=nucA;
                ncData(nci).toFix(n)=0;
                changes=changes+1;
            end
            %In theory, there should be an else here. I didn't write it,
            %and it may not be necessary because anaphase really does
            %create that much movement. The only time this doesn't happen
            %it seems is if cell division fails, and then that nuclei
            %will be and dropped from the data set
            %However, you will notice it remains in toFix.
            
        end
    end
    %Displays the result of this loop. 
    disp('Split nuclei:')
    disp(changes)
    disp('Remaining toFix:')
    for nci=1:length(ncData)-1
        fix=ncData(nci).toFix;
        disp(fix)
    end
end

%Save the schnitzcells variable to the _lin file
save([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'],'schnitzcells')


