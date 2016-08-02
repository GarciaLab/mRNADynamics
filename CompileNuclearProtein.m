function CompileNuclearProtein(varargin)
%THis code gives me the input
%varargin Variable length input argument list.
%allows any number of arguments to a function.  The variable
%varargin is a cell array containing the optional arguments to the
%function.  varargin must be declared as the last input argument
%and collects all the inputs from that point onwards.


%This function will add fluorescence information to each schnitz.

close all

%gives the location of these 5 quantities.
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders;

%Look at the input parameter and use defaults if missing
if isempty(varargin)% returns 1 if it is an empty array and 0 otherwise.an empty array has prod(size(X))==0
    FolderTemp=uigetdir(DropboxFolder,'Select folder with data to analyze');
    Dashes=strfind(FolderTemp,'\');
    Prefix=FolderTemp((Dashes(end)+1):end);
else
    Prefix=varargin{1};
end
FilePrefix=[Prefix,'_'];

[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders(Prefix);




%Load all the information
load([DropboxFolder,filesep,Prefix,'\Ellipses.mat'])
load([DropboxFolder,filesep,Prefix,'\FrameInfo.mat'])
load([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'])


%What type of experiment are we dealing with? Get this out of
%MovieDatabase.xlsx

[XLSNum,XLSTxt,XLSRaw]=xlsread([DropboxFolder,filesep,'MovieDatabase.xlsx']);
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
ExperimentTypeColumn=find(strcmp(XLSRaw(1,:),'ExperimentType'));
ExperimentAxisColumn=find(strcmp(XLSRaw(1,:),'ExperimentAxis'));
APResolutionColumn = find(strcmp(XLSRaw(1,:),'APResolution'));

DataFolderColumn=find(strcmp(XLSRaw(1,:),'DataFolder'));
Dashes=findstr(Prefix,'-');
PrefixRow=find(strcmp(XLSRaw(:,DataFolderColumn),[Prefix(1:Dashes(3)-1),'\',Prefix(Dashes(3)+1:end)]));
    if isempty(PrefixRow)
        PrefixRow=find(strcmp(XLSRaw(:,DataFolderColumn),[Prefix(1:Dashes(3)-1),'/',Prefix(Dashes(3)+1:end)]));
        if isempty(PrefixRow)
            error('Could not find data set in MovieDatabase.XLSX. Check if it is defined there.')
        end
    end
        
if isempty(PrefixRow)
    error('Entry not found in MovieDatabase.xlsx')
end

ExperimentType=XLSRaw{PrefixRow,ExperimentTypeColumn};
ExperimentAxis=XLSRaw{PrefixRow,ExperimentAxisColumn};
APResolution = XLSRaw{PrefixRow,APResolutionColumn};
Channel1=XLSTxt(PrefixRow,Channel1Column);
Channel2=XLSTxt(PrefixRow,Channel2Column)

nc9=cell2mat(XLSRaw(PrefixRow,nc9Column));
nc10=cell2mat(XLSRaw(PrefixRow,nc10Column));
nc11=cell2mat(XLSRaw(PrefixRow,nc11Column));
nc12=cell2mat(XLSRaw(PrefixRow,nc12Column));
nc13=cell2mat(XLSRaw(PrefixRow,nc13Column));
nc14=cell2mat(XLSRaw(PrefixRow,nc14Column));
%This is in case the last column for CF is all nan and is not part of
%the Num matrix
if ~isempty(CFColumn)    
    CF=cell2mat(XLSRaw(PrefixRow,CFColumn));
else
    CF=nan;
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

NewCyclePos=[nc9,nc10,nc11,nc12,nc13,nc14];
NewCyclePos=NewCyclePos(~(NewCyclePos==0));
NewCyclePos=NewCyclePos(~isnan(NewCyclePos));




%Add the APPosition to Particles if they don't exist yet
if (~isfield(schnitzcells,'APpos'))&(strcmp(lower(ExperimentAxis),'ap'))
    %error('This part of the code still needs to be implemented')
    AddNuclearPosition(Prefix)
    load([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'])
end



% Put together CompiledNuclei

%Get the actual time corresponding to each frame
if isfield(FrameInfo,'FileMode')
    if strcmp(FrameInfo(1).FileMode,'TIF')%Is this a TIF file if not is it a LSM or LIFE Export
        for j=1:length(FrameInfo)
            ElapsedTime(j)=etime(datevec(FrameInfo(j).TimeString),datevec(FrameInfo(1).TimeString));
        end
    elseif strcmp(FrameInfo(1).FileMode,'LSM')|strcmp(FrameInfo(1).FileMode,'LIFExport')%If it is a LIFEexport Sum over Fram e info
        for j=1:length(FrameInfo)
            ElapsedTime(j)=FrameInfo(j).Time-FrameInfo(1).Time;%Finds the elapsed time by subtracting each time point by the initial time point
        end
    else
        error('File mode not supported. Cannot extract time information. Include format in ExportDataForFISH.m')
    end
else
    warning('No FileMode information found. Assuming that this is TIF from the 2-photon.')
    for j=1:length(FrameInfo)
        ElapsedTime(j)=etime(datevec(FrameInfo(j).TimeString),datevec(FrameInfo(1).TimeString));
    end
end
    

ElapsedTime=ElapsedTime/60;     %Time is in minutes

%If there is no Approved field then create it
if ~isfield(schnitzcells,'Approved')
    for i=1:length(schnitzcells)
        schnitzcells(i).Approved=logical(1);
    end
end

%If there is no FrameApproved field then create it
if ~isfield(schnitzcells,'FrameApproved')
    for i=1:length(schnitzcells)
        schnitzcells(i).FrameApproved=logical(ones(size(schnitzcells(i).frames)));
    end
end




%Now get the nuclear information for those that were approved
NZSclices=size(schnitzcells(1).Fluo,2);

%CompiledNuclei(length(schnitzcells))=struct;

h=waitbar(0,'Compiling nuclear traces');
k=1;
for i=1:length(schnitzcells)
    

    waitbar(i/length(schnitzcells),h)
    
    if (schnitzcells(i).Approved==1)
        %Which frames were approved manually?
        FrameFilter=schnitzcells(i).FrameApproved;
        
        %Check that for the remaining frames we got a good z-profile
        for j=1:length(schnitzcells(i).frames)
            [MaxValue,MaxPos]=max(schnitzcells(i).Fluo(j,:));
            if (MaxPos==2)|(MaxPos==NZSclices-1)
                FrameFilter(j)=0;
            end
        end
        
        if sum(FrameFilter)
        
            %Copy the filtered information
            CompiledNuclei(k).P=schnitzcells(i).P;
            CompiledNuclei(k).E=schnitzcells(i).E;
            CompiledNuclei(k).D=schnitzcells(i).D;
            CompiledNuclei(k).Frames=schnitzcells(i).frames(FrameFilter);
            CompiledNuclei(k).xPos=schnitzcells(i).cenx(FrameFilter);
            CompiledNuclei(k).yPos=schnitzcells(i).ceny(FrameFilter);
            CompiledNuclei(k).Radius=schnitzcells(i).len(FrameFilter);
            CompiledNuclei(k).cellno=schnitzcells(i).cellno(FrameFilter);
            CompiledNuclei(k).nc=[];
            
            %Save the information about the original schnitz
            CompiledNuclei(k).schnitz=i;

            
            if strcmp(lower(ExperimentAxis),'ap')
                %Determine the particles average and median AP position
                CompiledNuclei(k).MeanAP=mean(schnitzcells(i).APpos(FrameFilter));
                CompiledNuclei(k).MedianAP=median(schnitzcells(i).APpos(FrameFilter));
            end


            
            %Copy and extract the fluorescence information
            CompiledNuclei(k).FluoMax=max(schnitzcells(i).Fluo(FrameFilter,:),[],2);

            k=k+1;
        end
    end
end
close(h)     


%% Create filters

%nc filters:

%ncFilterID just tells you the identity of the different
%filters stored in the cell ncFilter
ncFilterID=[];
if nc9~=0
    ncFilterID=9;
end
if nc10~=0
    ncFilterID=[ncFilterID,10];
end
if nc11~=0
    ncFilterID=[ncFilterID,11];
end
if nc12~=0
    ncFilterID=[ncFilterID,12];
end
if nc13~=0
    ncFilterID=[ncFilterID,13];
end
if nc14~=0
    ncFilterID=[ncFilterID,14];
end
%Add the first nc
ncFilterID=[min(ncFilterID)-1,ncFilterID];


%Create the filter
ncFilter=logical(zeros(length(CompiledNuclei),length(ncFilterID)));
for i=1:length(CompiledNuclei)
    if ~isempty(CompiledNuclei(i).Frames)
    
        if ~isempty(CompiledNuclei(i).nc)
            ncFilter(i,find(CompiledNuclei(i).nc==ncFilterID))=logical(1);
        else
            ncsFound=find(CompiledNuclei(i).Frames(1)>=[nc9,nc10,nc11,nc12,nc13,nc14]);
            if ncsFound(end)==1
                CompiledNuclei(i).nc=9;
                ncFilter(i,ncFilterID==9)=logical(1);
            elseif ncsFound(end)==2
                CompiledNuclei(i).nc=10;
                ncFilter(i,ncFilterID==10)=logical(1);
            elseif ncsFound(end)==3
                CompiledNuclei(i).nc=11;
                ncFilter(i,ncFilterID==11)=logical(1);
            elseif ncsFound(end)==4
                CompiledNuclei(i).nc=12;
                ncFilter(i,ncFilterID==12)=logical(1);
            elseif ncsFound(end)==5
                CompiledNuclei(i).nc=13;
                ncFilter(i,ncFilterID==13)=logical(1);
            elseif ncsFound(end)==6
                CompiledNuclei(i).nc=14;
                ncFilter(i,ncFilterID==14)=logical(1);
            end

        end
    end
end





if strcmp(lower(ExperimentAxis),'ap')
    %AP filters:

    %Divide the AP axis into boxes of a certain AP size. We'll see which
    %particle falls where.

    APResolution=APResolution;
    APbinID=0:APResolution:1;

    APFilter=logical(zeros(length(CompiledNuclei),length(APbinID)));
    for i=1:length(CompiledNuclei)
        APFilter(i,max(find(APbinID<=CompiledNuclei(i).MeanAP)))=1;
    end
end





%% Binning and averaging data




if strcmp(lower(ExperimentAxis),'ap')

   
    %Get the data for the individual particles in a matrix that has the frame
    %number and the particle number as dimensions. Also, get a vector that
    %reports the mean AP position.
    [AllTracesVector,AllTracesAP]=AllTracesNuclei(FrameInfo,CompiledNuclei);

    %Mean plot for different AP positions

    %Figure out the AP range to use
    MinAPIndex=1;%min(find(sum(APFilter)));
    MaxAPIndex=size(APFilter,2);%max(find(sum(APFilter)));

    %Get the corresponding mean information
    k=1;
    for i=MinAPIndex:MaxAPIndex
        [MeanVectorAPTemp,SDVectorAPTemp,NParticlesAPTemp]=AverageTracesNuclei(FrameInfo,...
            CompiledNuclei(APFilter(:,i)));
        MeanVectorAPCell{k}=MeanVectorAPTemp';
        SDVectorAPCell{k}=SDVectorAPTemp';
        NParticlesAPCell{k}=NParticlesAPTemp';
        k=k+1;
    end
    MeanVectorAP=cell2mat(MeanVectorAPCell);
    SDVectorAP=cell2mat(SDVectorAPCell);
    NParticlesAP=cell2mat(NParticlesAPCell);
end
    

%Calculate the mean for all of them
[MeanVectorAll,SDVectorAll,NParticlesAll]=AverageTracesNuclei(FrameInfo,CompiledNuclei);
%Now find the different maxima in each nc

MaxFrame=[];
for i=1:length(NewCyclePos)
    if i==1
        [Dummy,MaxIndex]=max(MeanVectorAll(1:NewCyclePos(1)));
        MaxFrame=[MaxFrame,MaxIndex];
    elseif i<=length(NewCyclePos)
        [Dummy,MaxIndex]=max(MeanVectorAll(NewCyclePos(i-1):NewCyclePos(i)));
        MaxFrame=[MaxFrame,NewCyclePos(i-1)+MaxIndex-1];
    end
end
[Dummy,MaxIndex]=max(MeanVectorAll(NewCyclePos(i):end));
MaxFrame=[MaxFrame,NewCyclePos(i)+MaxIndex-1];





%% Save everything


if strcmp(lower(ExperimentAxis),'ap')
    save([DropboxFolder,filesep,Prefix,filesep,'CompiledNuclei.mat'],...
        'CompiledNuclei','ElapsedTime','NewCyclePos','nc9','nc10','nc11',...
        'nc12','nc13','nc14','ncFilterID','ncFilter','APbinID','APFilter',...
        'MeanVectorAP','SDVectorAP','NParticlesAP','MeanVectorAll',...
        'SDVectorAll','NParticlesAll','MaxFrame','MinAPIndex','MaxAPIndex',...
        'AllTracesVector','AllTracesAP')
else
    save([DropboxFolder,filesep,Prefix,filesep,'CompiledNuclei.mat'],...
        'CompiledNuclei','ElapsedTime','NewCyclePos','nc9','nc10','nc11',...
        'nc12','nc13','nc14','ncFilterID','ncFilter','MeanVectorAll',...
        'SDVectorAll','NParticlesAll','MaxFrame')
end


save([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'],'schnitzcells')


