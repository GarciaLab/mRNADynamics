function CompileParticles(varargin)




%Parameters:
%First, the prefix.
%There after:
%ForceAP  -  Force AP detection even if it's there already.
%SkipTraces - Don't output the individual traces.
%SkipFluctuations - Don't generate the plots of the correlation of signal
%                   and offset
%SkipFits - Don't do the fits
%SkipMovie - Don't do the movie
%SkipAll - Skip all htat can be skipped
%ApproveAll - Approves all particles. This is useful if we want to do a
%             quick check of, for example, the AP profiles

%This function puts togetether all the information we have about particles.
%Things we want in here are:

close all

%Information about about folders
[SourcePath,FISHPath,DefaultDropboxFolder,MS2CodePath,SchnitzcellsFolder]=...
    DetermineLocalFolders;

%Look at the input parameters and use defaults if missing
%Prefix=[];
ForceAP=0;      %Force AP detection even if it's already there
SkipTraces=0;   %Do not output the individual traces.
SkipFluctuations=0;  %Do not generate the plots of correlations of fluctuations and offset
SkipFits=0;         %Do not generate the fit output (but still does the fit)
SkipMovie=0;        %Do not generate the movie
ApproveAll=0;

if isempty(varargin)
    FolderTemp=uigetdir(DefaultDropboxFolder,'Select folder with data to analyze');
    Dashes=strfind(FolderTemp,'\');
    Prefix=FolderTemp((Dashes(end)+1):end);
else
    Prefix=varargin{1};

    for i=2:length(varargin)
        if strcmp(varargin{i},'ForceAP')
            ForceAP=1;
        elseif strcmp(varargin{i},'SkipTraces')
            SkipTraces=1;
        elseif strcmp(varargin{i},'SkipFluctuations')
            SkipFluctuations=1;
        elseif strcmp(varargin{i},'SkipFits')    
            SkipFits=1;
        elseif strcmp(varargin{i},'SkipMovie')    
            SkipMovie=1;
        elseif strcmp(varargin{i},'SkipAll')        
            SkipTraces=1;
            SkipFluctuations=1;
            SkipFits=1;
            SkipMovie=1;
        elseif strcmp(varargin{i},'ApproveAll')    
            ApproveAll=1;            
        end
    end

end
FilePrefix=[Prefix,'_'];

%Now get the actual Dropbox folder
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,SchnitzcellsFolder]=...
    DetermineLocalFolders(Prefix);



%What type of experiment are we dealing with? Get this out of
%MovieDatabase.xlsx

[XLSNum,XLSTxt,XLSRaw]=xlsread([DefaultDropboxFolder,filesep,'MovieDatabase.xlsx']);
ExperimentTypeColumn=find(strcmp(XLSRaw(1,:),'ExperimentType'));
ExperimentAxisColumn=find(strcmp(XLSRaw(1,:),'ExperimentAxis'));

DataFolderColumn=find(strcmp(XLSRaw(1,:),'DataFolder'));
Dashes=findstr(Prefix,'-');
PrefixRow=find(strcmp(XLSRaw(:,DataFolderColumn),[Prefix(1:Dashes(3)-1),'\',Prefix(Dashes(3)+1:end)]));

if isempty(PrefixRow)
    error('Entry not found in MovieDatabase.xlsx')
end

ExperimentType=XLSRaw{PrefixRow,ExperimentTypeColumn};
ExperimentAxis=XLSRaw{PrefixRow,ExperimentAxisColumn};



%Load all the information
load([DropboxFolder,filesep,Prefix,filesep,'Particles.mat'])
%Check that FrameInfo exists
if exist([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'])
    load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'])
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


%Delete the files in folder where we'll write again.
if ~SkipTraces
    delete([DropboxFolder,filesep,Prefix,filesep,'ParticleTraces',filesep,'*.*'])
end
if ~SkipFits
    delete([DropboxFolder,filesep,Prefix,filesep,'Fits',filesep,'*.*'])
end



%See if we had any lineage/nuclear information
if exist([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'])
    load([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'])
    HistoneChannel=1;
else
    display('No lineage / nuclear information found. Proceeding without it.');
    HistoneChannel=0;
end




%Load the information about the nc from the XLS file
[Num,Txt]=xlsread([DefaultDropboxFolder,filesep,'MovieDatabase.xlsx']);
XLSHeaders=Txt(1,:);
Txt=Txt(2:end,:);

%Find the different columns.
DataFolderColumn=7;

%Convert the prefix into the string used in the XLS file
Dashes=findstr(FilePrefix(1:end-1),'-');



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

%Convert the prefix into the string used in the XLS file
Dashes=findstr(Prefix,'-');

%Find the corresponding entry in the XLS file
if (~isempty(findstr(Prefix,'Bcd')))&(isempty(findstr(Prefix,'BcdE1')))&...
        (isempty(findstr(Prefix,'NoBcd')))
    XLSEntry=find(strcmp(Txt(:,DataFolderColumn),...
        [Date,'\BcdGFP-HisRFP']));
else
    XLSEntry=find(strcmp(XLSRaw(:,DataFolderColumn),...
        [Prefix(1:Dashes(3)-1),'\',Prefix(Dashes(3)+1:end)]));
end


if strcmp(XLSRaw(XLSEntry,Channel2Column),'His-RFP')
    nc9=cell2mat(XLSRaw(XLSEntry,nc9Column));
    nc10=cell2mat(XLSRaw(XLSEntry,nc10Column));
    nc11=cell2mat(XLSRaw(XLSEntry,nc11Column));
    nc12=cell2mat(XLSRaw(XLSEntry,nc12Column));
    nc13=cell2mat(XLSRaw(XLSEntry,nc13Column));
    nc14=cell2mat(XLSRaw(XLSEntry,nc14Column));
    %This is in case the last column for CF is all nan and is not part of
    %the Num matrix
    if ~isempty(CFColumn)    
        CF=cell2mat(XLSRaw(XLSEntry,CFColumn));
    else
        CF=nan;
    end
else
    error('nc information not define in MovieDatabase.xlsx')
end

% Read in which end the stem loops are at, if this information is available
% (ES 2014-03-20)
StemLoopEndColumn = find(strcmp(XLSRaw(1, :), 'StemLoopEnd'));
if ~isempty(StemLoopEndColumn)
    StemLoopEnd = XLSRaw{XLSEntry, StemLoopEndColumn};
else
    StemLoopEnd = '';
end


NewCyclePos=[nc9,nc10,nc11,nc12,nc13,nc14];
NewCyclePos=NewCyclePos(~(NewCyclePos==0));
NewCyclePos=NewCyclePos(~isnan(NewCyclePos));



%Add the APPosition to Particles if they don't exist yet. Do this only if
%we took AP data. Otherwise just add XY.

if strcmp(ExperimentAxis,'AP')
    if (~isfield(Particles,'APpos'))|ForceAP
        if HistoneChannel
            AddParticlePosition(Prefix);
        else
            AddParticlePosition(Prefix,'SkipAlignment')
        end

    else
        display('Using saved AP information')
    end   
else
    AddParticlePosition(Prefix,'NoAP');
end

    
load([DropboxFolder,filesep,Prefix,filesep,'Particles.mat'])

if HistoneChannel
    load([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'])
end


%Folders for reports
if strcmp(ExperimentAxis,'AP')
    mkdir([DropboxFolder,filesep,Prefix,filesep,'APMovie'])
end
mkdir([DropboxFolder,filesep,Prefix,filesep,'ParticleTraces'])
mkdir([DropboxFolder,filesep,Prefix,filesep,'TracesFluctuations'])
mkdir([DropboxFolder,filesep,Prefix,filesep,'Offset'])
mkdir([DropboxFolder,filesep,Prefix,filesep,'Fits'])
mkdir([DropboxFolder,filesep,Prefix,filesep,'Probabilities'])
mkdir([DropboxFolder,filesep,Prefix,filesep,'Various'])


%% Put together CompiledParticles


%Approve all particles if the mode has been selected
if ApproveAll
    for i=1:length(Particles)
        %Make sure the particle has an associated nucleus if we are in
        %HistoneChannel mode
        if HistoneChannel
            if ~isempty(Particles(i).Nucleus)
                %If a particle has been explicitly rejected then don't
                %approve it!
                if Particles(i).Approved~=-1
                    Particles(i).Approved=1;
                end
            end
        else
            Particles(i).Approved=1;
        end
    end
end


if HistoneChannel
    %Get information about which particle is associated to which nucleus in
    %schnitzcells
    for i=1:length(Particles)
        if ~isempty(Particles(i).Nucleus)
            AssignedNuclei(i)=Particles(i).Nucleus;
        else
            AssignedNuclei(i)=nan;
        end
    end
end


if HistoneChannel&strcmp(ExperimentAxis,'AP')
    %First, figure out the AP position of each of the nuclei.

    %Load the AP detection information
    load([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'])
    %Angle between the x-axis and the AP-axis
    APAngle=atan((coordPZoom(2)-coordAZoom(2))/(coordPZoom(1)-coordAZoom(1)));
    APLength=sqrt((coordPZoom(2)-coordAZoom(2))^2+(coordPZoom(1)-coordAZoom(1))^2);


    %The information in Ellipses is
    %(x, y, a, b, theta, maxcontourvalue, time, particle_id)
    for i=1:length(Ellipses)
        for j=1:size(Ellipses{i},1)

            %Angle between the x-axis and the particle using the A position as a
            %zero
            Angles=atan((Ellipses{i}(j,2)-coordAZoom(2))./(Ellipses{i}(j,1)-coordAZoom(1)));

            %Distance between the points and the A point
            Distances=sqrt((coordAZoom(2)-Ellipses{i}(j,2)).^2+(coordAZoom(1)-Ellipses{i}(j,1)).^2);
            APPositions=Distances.*cos(Angles-APAngle);
            EllipsePos{i}(j)=APPositions/APLength;
        end
    end
end


%Get the actual time corresponding to each frame
if isfield(FrameInfo,'FileMode')
    if strcmp(FrameInfo(end).FileMode,'TIF')
        for j=1:length(FrameInfo)
            ElapsedTime(j)=etime(datevec(FrameInfo(j).TimeString),datevec(FrameInfo(1).TimeString));
        end
    elseif strcmp(FrameInfo(end).FileMode,'LSM')
        for j=1:length(FrameInfo)
            ElapsedTime(j)=FrameInfo(j).Time-FrameInfo(1).Time;
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
    

%Some parameters
IntArea=109;        %Are of integration
MinAPArea=12500;%700;    %Minimum area in pixels in order to consider an AP bin as valid.


if strcmp(ExperimentAxis,'AP')
    %Divide the image into AP bins. The size of the bin will depend on the
    %experiment
    if strfind(lower(Prefix),'eve')     %Eve2 experiments
        APResolution=0.01;
    %hb or kni BAC experiments
    elseif ~isempty(strfind(lower(Prefix),'hbbac'))|...
            ~isempty(strfind(lower(Prefix),'knibac'))     
        APResolution=0.015;
    else                                %All other experiments
        APResolution=0.025;
    end
       
    APbinID=0:APResolution:1;

    %Create an image for the different AP bins
    APPosImage=zeros(FrameInfo(1).LinesPerFrame,FrameInfo(1).PixelsPerLine);
    [Rows,Columns]=size(APPosImage);

    for i=1:Rows
        for j=1:Columns
            Angle=atan((i-coordAZoom(2))./(j-coordAZoom(1)));
            Distance=sqrt((coordAZoom(2)-i).^2+(coordAZoom(1)-j).^2);
            APPosition=Distance.*cos(Angle-APAngle);
            APPosImage(i,j)=APPosition/APLength;
        end
    end


    APPosBinImage=zeros(size(APPosImage));
    for i=1:(length(APbinID)-1)
        FilteredMask=(APbinID(i)<=APPosImage)&(APbinID(i+1)>APPosImage);
        APPosBinImage=APPosBinImage+FilteredMask*i;
    end


    %Calculate the area in pixels corresponding to each AP bin. We will use
    %this to get rid of small AP bins in the image and also to calculate
    %probabilities of nuclei being active.
    for i=1:length(APbinID)
        APbinArea(i)=sum(sum(APPosBinImage==i));

        %Discard anything that is below MinAPArea
        if APbinArea(i)<(MinAPArea/0.025*APResolution)
            APbinArea(i)=nan;
        end
    end
end


%Now get the particle information for those that were approved
k=1;
h=waitbar(0,'Compiling traces');
for i=1:length(Particles)

      
    waitbar(i/length(Particles),h)
    if (Particles(i).Approved==1)
 
        
        %Which frames were approved manually?
        FrameFilter=Particles(i).FrameApproved;
        %What is the first frame that was found, regardless of the column
        %condition?
        FirstFrame=Particles(i).Frame(min(find(Particles(i).FrameApproved)));
        
        %Check that for the remaining frames we got a good z-profile
        for j=1:length(Particles(i).Frame)
            ZProfile=fad.channels(Particles(i).Frame(j)).fits.shadowsDog{Particles(i).Index(j)};
            [Dummy,ZMax]=max(ZProfile);
            if (ZMax==1)|(ZMax==length(ZProfile))
                FrameFilter(j)=0;
            end
        end
        
        %Should I only keep traces of a certain length? We also just keep
        %the ones that have a real schnitz associated with them
        AnalyzeThisParticle=1;      %Flag to see if this particle should be analyzed.
        
        if HistoneChannel
            if ~((sum(FrameFilter)>0)&...
                (~isempty(schnitzcells(Particles(i).Nucleus).frames)))
                AnalyzeThisParticle=0;
            end
        elseif ~(sum(FrameFilter)>0)
            AnalyzeThisParticle=0;
        end
        
        
        %See if this particle is in one of the approved AP bins
        if strcmp(ExperimentAxis,'AP')
            CurrentAPbin=max(find(APbinID<mean(Particles(i).APpos(FrameFilter))));
            if isnan(APbinArea(CurrentAPbin))
                AnalyzeThisParticle=0;
            end
        end
        
        
        
        if AnalyzeThisParticle
        
            %Reference to the original Particles index
            CompiledParticles(k).OriginalParticle=i;            
            
            %Copy the filtered information
            CompiledParticles(k).Frame=Particles(i).Frame(FrameFilter);
            CompiledParticles(k).Index=Particles(i).Index(FrameFilter);
            CompiledParticles(k).xPos=Particles(i).xPos(FrameFilter);
            CompiledParticles(k).yPos=Particles(i).yPos(FrameFilter);
            
            if strcmp(ExperimentAxis,'AP')
                CompiledParticles(k).APpos=Particles(i).APpos(FrameFilter);
                
                %Determine the particles average and median AP position
                CompiledParticles(k).MeanAP=mean(Particles(i).APpos(FrameFilter));
                CompiledParticles(k).MedianAP=median(Particles(i).APpos(FrameFilter));
            elseif strcmp(ExperimentAxis,'DV')&isfield(Particles,'APpos')
                %AP information:
                CompiledParticles(k).APpos=Particles(i).APpos(FrameFilter);
                CompiledParticles(k).MeanAP=mean(Particles(i).APpos(FrameFilter));
                CompiledParticles(k).MedianAP=median(Particles(i).APpos(FrameFilter));
                %DV information:
                CompiledParticles(k).DVpos=Particles(i).DVpos(FrameFilter);
                CompiledParticles(k).MeanDV=mean(Particles(i).DVpos(FrameFilter));
                CompiledParticles(k).MedianDV=median(Particles(i).DVpos(FrameFilter));

            end
            
            %If we have the histone channel we will actually replace the AP
            %position by the position of the nucleus where the particle was
            %found. If there is on nucleus (like when a particle survives
            %past the nuclear division) we will still use the actual particle
            %position.
            if HistoneChannel&strcmp(ExperimentAxis,'AP')
                %Save the original particle position
                CompiledParticles(k).APposParticle=CompiledParticles(k).APpos;                
                
                FramesToCheck=schnitzcells(Particles(i).Nucleus).frames(...
                    ismember(schnitzcells(Particles(i).Nucleus).frames,Particles(i).Frame(FrameFilter)));
                EllipsesToCheck=schnitzcells(Particles(i).Nucleus).cellno(...
                    ismember(schnitzcells(Particles(i).Nucleus).frames,Particles(i).Frame(FrameFilter)));
               
                for j=1:length(FramesToCheck)
                    IndexToChange=find(CompiledParticles(k).Frame==FramesToCheck(j));
                    CompiledParticles(k).APPos(IndexToChange)=EllipsePos{FramesToCheck(j)}(EllipsesToCheck(j));
                end
            end

            %First frame it was detected at
            CompiledParticles(k).FirstFrame=FirstFrame;

            CompiledParticles(k).Approved=Particles(i).Approved;
            
            %Copy the fit results if they are there
            if isfield(Particles,'Fit')
                CompiledParticles(k).Fit=Particles(i).Fit;
            end
 
            
            %Extract information form fad about fluorescence and background
            [Frame,Amp,Off,Off2,Amp2,AmpOld,AmpRaw,Error,optFit1,FitType]=GetParticleTrace(k,CompiledParticles,fad);
            CompiledParticles(k).Fluo=Amp;
            CompiledParticles(k).Off=Off;
            CompiledParticles(k).Off2=Off2;
            CompiledParticles(k).Fluo2=Amp2;
            CompiledParticles(k).FluoOld=AmpOld;
            CompiledParticles(k).FluoRaw=AmpRaw;
            CompiledParticles(k).FluoError=Error;
            CompiledParticles(k).optFit1=optFit1;
            CompiledParticles(k).FitType=FitType;


            
            
            %Determine the nc where this particle was born
            CompiledParticles(k).nc=FrameInfo(CompiledParticles(k).Frame(1)).nc;

            
            if HistoneChannel
                CompiledParticles(k).Nucleus=Particles(i).Nucleus;

                %Save lineage information in terms of particles
                if ~isempty(schnitzcells(Particles(i).Nucleus).P)
                    if isempty(find(AssignedNuclei==schnitzcells(Particles(i).Nucleus).P))
                        CompiledParticles(k).PParticle=0;
                    else
                        CompiledParticles(k).PParticle=find(AssignedNuclei==schnitzcells(Particles(i).Nucleus).P);
                    end
                else
                    CompiledParticles(k).PParticle=[];
                end

                if ~isempty(schnitzcells(Particles(i).Nucleus).D)
                    if isempty(find(AssignedNuclei==schnitzcells(Particles(i).Nucleus).D))
                        CompiledParticles(k).DParticle=0;
                    else
                        CompiledParticles(k).DParticle=find(AssignedNuclei==schnitzcells(Particles(i).Nucleus).D);
                    end
                else
                    CompiledParticles(k).DParticle=[];
                end

                if ~isempty(schnitzcells(Particles(i).Nucleus).E)
                    if isempty(find(AssignedNuclei==schnitzcells(Particles(i).Nucleus).E))
                        CompiledParticles(k).EParticle=0;
                    else
                        CompiledParticles(k).EParticle=find(AssignedNuclei==schnitzcells(Particles(i).Nucleus).E);
                    end
                else
                    CompiledParticles(k).EParticle=[];
                end

                %Save information about the nucleus birth and death
                CompiledParticles(k).NucStart=schnitzcells(Particles(i).Nucleus).frames(1);
                CompiledParticles(k).NucEnd=schnitzcells(Particles(i).Nucleus).frames(end);
                
   
            end

            
            
            %Plot and save this trace together with its offset value
            
            if ~SkipTraces     
                
                figure(2)
                %Size of the snippet for each frame
                SnippetSize=31;
                %Width of the screen
                ScreenWidth=get( 0, 'ScreenSize' );
                ScreenWidth=ScreenWidth(3);

                %Figure out the arrangement of snippets
                NFrames=length(CompiledParticles(k).Frame);

                NRows=ceil(NFrames/ScreenWidth*SnippetSize)*2;
                NCols=max([2,ceil(NFrames/NRows)]);

                %Actual total number of rows
                TotalRows=14;

                subplot(TotalRows,NCols,[1:((TotalRows-NRows)*NCols)])



                %Top left plot
                FilterMatrix=zeros((TotalRows-NRows),NCols);
                FilterMatrix(:,1:ceil(NCols/2))=1;
                subplot(TotalRows,NCols,find(FilterMatrix'))

                errorbar(ElapsedTime(CompiledParticles(k).Frame),...
                    CompiledParticles(k).Fluo,ones(size(CompiledParticles(k).Fluo))*...
                    CompiledParticles(k).FluoError,...
                    '.-r');
                hold on


                plot(ElapsedTime(CompiledParticles(k).Frame),...
                    CompiledParticles(k).Off*IntArea,'.-g');
                if ~isempty(CompiledParticles(k).optFit1)

                    if strcmp(CompiledParticles(k).FitType,'spline')
                        SplineValues=ppval(CompiledParticles(k).optFit1,double(CompiledParticles(k).Frame));
                    elseif strcmp(CompiledParticles(k).FitType,'mean')
                        SplineValues=ones(size(CompiledParticles(k).Frame))*CompiledParticles(k).optFit1;
                    elseif strcmp(CompiledParticles(k).FitType,'line')
                        SplineValues=polyval(CompiledParticles(k).optFit1,CompiledParticles(k).Frame);     
                    end

                    plot([ElapsedTime(CompiledParticles(k).Frame)],SplineValues*IntArea,'-k')
                    title(['Particle ',num2str(k),'(',num2str(i),'), nc',num2str(CompiledParticles(k).nc)])
                else
                    title(['Particle ',num2str(k),'(',num2str(i),'), nc',num2str(CompiledParticles(k).nc),...
                        ' - WARNING: No offset fit'])
                end
                hold off
                %legend({'Particle','Offset','Offset fit'},'Location','Best')
                xlabel('Time (min)')
                ylabel('Fluorescence (au)')
                axis square
                set(gca, 'Position', get(gca, 'OuterPosition') - ...
                    get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);

                drawnow




                %Top right plot
                if HistoneChannel
                    subplot(TotalRows,NCols,find(~FilterMatrix'))

                    if length(CompiledParticles(k).Frame)>1
                        colormap(jet(128));
                        cmap=colormap;

                        ColorTime=[];
                        for j=1:length(CompiledParticles(k).Frame)
                            ColorTime(j,:)= cmap(round((j-1)/...
                                (length(CompiledParticles(k).Frame)-1)*127+1),:);
                        end
                    else
                        ColorTime=[];
                        ColorTime(1,:)=[1,0,0];
                    end



                    hold on
                    for j=1:length(CompiledParticles(k).Frame)
                        PosSchnitz=find((schnitzcells(CompiledParticles(k).Nucleus).frames)==...
                            CompiledParticles(k).Frame(j));
                        PosEllipse=schnitzcells(CompiledParticles(k).Nucleus).cellno(PosSchnitz);
                        CurrEllipse=Ellipses{CompiledParticles(k).Frame(j)}(PosEllipse,:);

                        if ~isempty(CurrEllipse)
                            EllipseHandle=ellipse(CurrEllipse(3),...
                                CurrEllipse(4),...
                                CurrEllipse(5),...
                                0,0);
                            set(EllipseHandle,'color',ColorTime(j,:))
                            plot(CompiledParticles(k).xPos(j)-CurrEllipse(1),...
                                CompiledParticles(k).yPos(j)-CurrEllipse(2),'o','color',...
                                ColorTime(j,:))
                        else
                            PosSchnitz=length(schnitzcells(CompiledParticles(k).Nucleus).frames);
                            PosEllipse=schnitzcells(CompiledParticles(k).Nucleus).cellno(PosSchnitz);
                            CurrEllipse=...
                                Ellipses{schnitzcells(CompiledParticles(k).Nucleus).frames(PosSchnitz)}(PosEllipse,:);



                            plot(CompiledParticles(k).xPos(j)-CurrEllipse(1),...
                                CompiledParticles(k).yPos(j)-CurrEllipse(2),'o','color',...
                                ColorTime(j,:))
                        end



                    end
                    hold off
                    box on
                    xlabel('x position (pixels)')
                    ylabel('y position (pixels)')
                    axis square
                    set(gca, 'Position', get(gca, 'OuterPosition') - ...
                        get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
                    set(gca, 'YDir', 'reverse')
                    BarHandle = colorbar;
                    set(BarHandle,'YTick',[])
                    BarHandle=cbfreeze(BarHandle);
                    ylabel(BarHandle,'Time')
                    if strcmp(ExperimentAxis,'AP')
                        title(['Mean AP: ',num2str(CompiledParticles(k).MeanAP)])
                    end
                    drawnow
                end


                %Snippets
                for j=1:NFrames
                    subplot(TotalRows,NCols,(TotalRows-NRows)*NCols+j)

                    [x,y,z]=fad2xyzFit(CompiledParticles(k).Frame(j),fad, 'addMargin'); 
                    xTrace=x(CompiledParticles(k).Index(j));
                    yTrace=y(CompiledParticles(k).Index(j));
                    zTrace=z(CompiledParticles(k).Index(j));

                    Image=imread([FISHPath,filesep,'Data',filesep,FilePrefix(1:end-1),filesep,...
                        FilePrefix,iIndex(CompiledParticles(k).Frame(j),3),'_z',iIndex(zTrace,2),'.tif']);
                    [ImRows,ImCols]=size(Image);

                    ImageSnippet=zeros(SnippetSize,SnippetSize);



                    yRange=round(yTrace)+[-(SnippetSize-1)/2:(SnippetSize-1)/2];
                    yFilter=(yRange>0)&(yRange<=ImRows);


                    xRange=round(xTrace)+[-(SnippetSize-1)/2:(SnippetSize-1)/2];
                    xFilter=(xRange>0)&(xRange<=ImCols);



                    ImageSnippet(yFilter,xFilter)=Image(yRange(yFilter),...
                        xRange(xFilter));

                    imshow(ImageSnippet,[],'Border','Tight','InitialMagnification',200)
                    set(gca, 'Position', get(gca, 'OuterPosition') - ...
                        get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);

                    if HistoneChannel
                        %Plot the corresponding nucleus
                        CurrentSchnitz=schnitzcells(CompiledParticles(k).Nucleus);
                        if sum((CurrentSchnitz.frames)==CompiledParticles(k).Frame(j))==1
                            hold on
                            EllipseNumber=CurrentSchnitz.cellno(...
                                find((CurrentSchnitz.frames)==CompiledParticles(k).Frame(j)));

                            CurrEllipse=Ellipses{CompiledParticles(k).Frame(j)}(EllipseNumber,:);

                            EllipseHandle=ellipse(CurrEllipse(3),...
                                CurrEllipse(4),...
                                CurrEllipse(5),...
                                CurrEllipse(1)-xTrace+(SnippetSize-1)/2,...
                                CurrEllipse(2)-yTrace+(SnippetSize-1)/2);
                            %set(EllipseHandle,'color',ColorTime(j,:))
                            set(EllipseHandle,'color','g')
                            hold off
                        end
                    end


                end
                set(gcf,'Position',[1,41,1280,684])  

                drawnow
                saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'ParticleTraces',filesep,iIndex(k,3),...
                    '(',num2str(i),')-nc',...
                    num2str(CompiledParticles(k).nc),'.tif'])
                close(2)
            end
            
           
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
ncFilter=logical(zeros(length(CompiledParticles),length(ncFilterID)));
for i=1:length(CompiledParticles)
    %Sometimes CompiledParticles(i).nc is empty. This is because of some
    %problem with FrameInfo! In that case we'll pull the information out of
    %the XLS file.
    if ~isempty(CompiledParticles(i).nc)
        ncFilter(i,find(CompiledParticles(i).nc==ncFilterID))=logical(1);
    else
        ncsFound=find(CompiledParticles(i).Frame(1)>=[nc9,nc10,nc11,nc12,nc13,nc14]);
        if ncsFound(end)==1
            CompiledParticles(i).nc=9;
            ncFilter(i,ncFilterID==9)=logical(1);
        elseif ncsFound(end)==2
            CompiledParticles(i).nc=10;
            ncFilter(i,ncFilterID==10)=logical(1);
        elseif ncsFound(end)==3
            CompiledParticles(i).nc=11;
            ncFilter(i,ncFilterID==11)=logical(1);
        elseif ncsFound(end)==4
            CompiledParticles(i).nc=12;
            ncFilter(i,ncFilterID==12)=logical(1);
        elseif ncsFound(end)==5
            CompiledParticles(i).nc=13;
            ncFilter(i,ncFilterID==13)=logical(1);
        elseif ncsFound(end)==6
            CompiledParticles(i).nc=14;
            ncFilter(i,ncFilterID==14)=logical(1);
        end
    
    end
end



%AP filters:
if strcmp(ExperimentAxis,'AP')
    %Divide the AP axis into boxes of a certain AP size. We'll see which
    %particle falls where.

    APFilter=logical(zeros(length(CompiledParticles),length(APbinID)));
    for i=1:length(CompiledParticles)
        APFilter(i,max(find(APbinID<=CompiledParticles(i).MeanAP)))=1;
    end
end


%% Binning and averaging data

%Get the data for the individual particles in a matrix that has the frame
%number and the particle number as dimensions. Also, get a vector that
%reports the mean AP position.
[AllTracesVector,AllTracesAP]=AllTraces(FrameInfo,CompiledParticles,'NoAP');

if strcmp(ExperimentAxis,'AP')
    %Mean plot for different AP positions

    %Figure out the AP range to use
    MinAPIndex=1;%min(find(sum(APFilter)));
    MaxAPIndex=size(APFilter,2);%max(find(sum(APFilter)));

    %Get the corresponding mean information
    k=1;
    for i=MinAPIndex:MaxAPIndex
        [MeanVectorAPTemp,SDVectorAPTemp,NParticlesAPTemp]=AverageTraces(FrameInfo,...
            CompiledParticles(APFilter(:,i)));
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
[MeanVectorAll,SDVectorAll,NParticlesAll]=AverageTraces(FrameInfo,CompiledParticles);
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
if ~isempty(NewCyclePos)        %Why is this empty sometimes?
                                %I think this only occurs with suboptimal
                                %data
    MaxFrame=[MaxFrame,NewCyclePos(i)+MaxIndex-1];
end




%% Instantaneous rate of change
FrameWindow=5;

%Calculate the derivative as a function of time for each particle
for j=1:length(CompiledParticles)
    [CompiledParticles(j).SlopeTrace,...
        CompiledParticles(j).SDSlopeTrace]=...
        TraceDerivative(CompiledParticles(j),...
        ElapsedTime,FrameWindow);
end


if strcmp(ExperimentAxis,'AP')
    %Calculate the average slope over an AP window
    MeanSlopeVectorAP=nan(size(MeanVectorAP));
    SDSlopeVectorAP=nan(size(MeanVectorAP));
    NSlopeAP=nan(size(MeanVectorAP));


    APBins=find(sum(APFilter));



    for j=1:length(APBins)
        TraceCell=cell(length(ElapsedTime),1);

        ParticlesToAverage=find(APFilter(:,APBins(j)));

        for k=1:length(ParticlesToAverage)

            for m=1:length(CompiledParticles(ParticlesToAverage(k)).SlopeTrace)

                TraceCell{CompiledParticles(ParticlesToAverage(k)).Frame(m)}=...
                    [TraceCell{CompiledParticles(ParticlesToAverage(k)).Frame(m)},...
                    CompiledParticles(ParticlesToAverage(k)).SlopeTrace(m)];
            end
        end

        %Get rid of the nan in certain time points
        TraceCell=cellfun(@(x) x(~isnan(x)),TraceCell,'UniformOutput',false);

        MeanTrace=cellfun(@mean,TraceCell,'UniformOutput',false);
        SDTrace=cellfun(@std,TraceCell,'UniformOutput',false);
        NParticlesTrace=cellfun(@length,TraceCell,'UniformOutput',false);

        MeanSlopeVectorAP(:,APBins(j))=[MeanTrace{:}];
        SDSlopeVectorAP(:,APBins(j))=[SDTrace{:}];
        NSlopeAP(:,APBins(j))=[NParticlesTrace{:}];
    end
end



%% Integrating each particle

%In order to take this seriously I need to come up with a way to deal with
%the particles that survive during mitosis.

for i=1:length(CompiledParticles)
    if length(ElapsedTime(CompiledParticles(i).Frame))>1
        CompiledParticles(i).TotalmRNA=trapz(ElapsedTime(CompiledParticles(i).Frame),CompiledParticles(i).Fluo);
        
        %Estimate the error
        if length(CompiledParticles(i).Frame)==2
            CompiledParticles(i).TotalmRNAError=(ElapsedTime(CompiledParticles(i).Frame(2))-...
                ElapsedTime(CompiledParticles(i).Frame(1)))/2*...
                CompiledParticles(i).FluoError*sqrt(2);
        else
              
            ErrorTemp=[];
            %Calculate the error of the inner points
            for j=2:(length(CompiledParticles(i).Frame)-1)
                 ErrorTemp(j)=(ElapsedTime(CompiledParticles(i).Frame(j+1))-...
                     ElapsedTime(CompiledParticles(i).Frame(j-1)))*...
                     CompiledParticles(i).FluoError;
            end
            
            %Calculate the error of the outer points
            ErrorTemp(1)=(ElapsedTime(CompiledParticles(i).Frame(2))-...
                ElapsedTime(CompiledParticles(i).Frame(1)))/2*...
                CompiledParticles(i).FluoError;
            
            ErrorTemp(length(CompiledParticles(i).Frame))=...
                (ElapsedTime(CompiledParticles(i).Frame(end))-...
                ElapsedTime(CompiledParticles(i).Frame(end-1)))/2*...
                CompiledParticles(i).FluoError;
            
            %Now, add it all up
            CompiledParticles(i).TotalmRNAError=sqrt(sum(ErrorTemp.^2));
   
            
        end
    else
        CompiledParticles(i).TotalmRNA=[];
    end
end


%% Information about the cytoplasm
%If the nuclear masks are present then use them. Otherwise just calculate
%the median of the images as a function of time
if HistoneChannel&strcmp(ExperimentAxis,'AP')
    [MeanCyto,SDCyto,MedianCyto,MaxCyto]=GetCytoMCP(Prefix);
else
    MeanCyto=[];
    SDCyto=[];
    MaxCyto=[];

    h=waitbar(0,'Calculating the median cyto intentisy');
    for i=1:length(FrameInfo)
        waitbar(i/length(FrameInfo),h)
        for j=1:FrameInfo(1).NumberSlices
            Image(:,:,j)=imread([FISHPath,filesep,'Data',filesep,Prefix,filesep,Prefix,'_',iIndex(i,3),'_z',iIndex(j,2),'.tif']);
        end
        ImageMax=max(Image,[],3);
        MedianCyto(i)=median(double(ImageMax(:)));
    end
    close(h)
end    
    
    
    
    

%% Offset and fluctuations


%Is there a correlation between the fluctuations coming from the offset and
%those observed in the traces? In order to figure this out I'll fit the
%nc13 and nc14 intensity data with splines and compute the deviations with
%respect to them. I'll look into different versions of the data such as
%with and without the offset subtracted


if ~SkipFluctuations

    IntArea=109;


    FilteredParticles=find(ncFilter(:,end)|ncFilter(:,end-1));

    OffsetFluct=[];
    DataRawFluct=[];
    DataOldFluct=[];
    DataSplineFluct=[];

    for j=1:length(FilteredParticles)

        try
            %Deviation from offset with respect to spline
            optFit = adaptiveSplineFit(double([ElapsedTime(CompiledParticles(FilteredParticles(j)).Frame)]),...
                    double([CompiledParticles(FilteredParticles(j)).Off*IntArea]),5);
            SplineValues=ppval(optFit,double([ElapsedTime(CompiledParticles(FilteredParticles(j)).Frame)]));    



            %Deviation of the raw data, without background subtraction, with
            %respect to a spline.
            DataFitRaw = adaptiveSplineFit(double([ElapsedTime(CompiledParticles(FilteredParticles(j)).Frame)]),...
                double(CompiledParticles(FilteredParticles(j)).FluoRaw),10);
            DataFitRawValues=ppval(DataFitRaw,double([ElapsedTime(CompiledParticles(FilteredParticles(j)).Frame)]));


            %Deviation of the raw data minus the actual offset with respect to a
            %spline
            DataFitOld = adaptiveSplineFit(double([ElapsedTime(CompiledParticles(FilteredParticles(j)).Frame)]),...
                double(CompiledParticles(FilteredParticles(j)).FluoOld),10);
            DataSplineValuesOld=ppval(DataFitOld,double([ElapsedTime(CompiledParticles(FilteredParticles(j)).Frame)]));



            %Deviation of the raw data minues the spline offset
            DataFit = adaptiveSplineFit(double([ElapsedTime(CompiledParticles(FilteredParticles(j)).Frame)]),...
                double(CompiledParticles(FilteredParticles(j)).Fluo),10);
            DataSplineValues=ppval(DataFit,double([ElapsedTime(CompiledParticles(FilteredParticles(j)).Frame)]));



            %Put all the data together for the plot
            OffsetFluct=[OffsetFluct,CompiledParticles(FilteredParticles(j)).Off*IntArea-SplineValues];
            DataRawFluct=[DataRawFluct,double(CompiledParticles(FilteredParticles(j)).FluoRaw)-DataFitRawValues];
            DataOldFluct=[DataOldFluct,double(CompiledParticles(FilteredParticles(j)).FluoOld)-DataSplineValuesOld];
            DataSplineFluct=[DataSplineFluct,double(CompiledParticles(FilteredParticles(j)).Fluo)-DataSplineValues];
        end
    end


    xRange=linspace(-4500,4500);

    figure(4)
    plot(OffsetFluct,DataRawFluct,'.k')
    xlabel('Offset fluctuation')
    ylabel('Fluctuations in raw data')
    axis square
    xlim([-4500,4500])
    ylim([-4500,4500])
    R = corrcoef(OffsetFluct,DataRawFluct);
    title(['Correlation: ',num2str(R(2,1))])
    saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'TracesFluctuations',filesep,'Fluct-OffsetVsRawData.tif'])

    figure(5)
    plot(OffsetFluct,DataOldFluct,'.k')
    xlabel('Offset fluctuation')
    ylabel('Fluctuations with instantaneous offset subtraction')
    axis square
    xlim([-4500,4500])
    ylim([-4500,4500])
    R = corrcoef(OffsetFluct,DataOldFluct)
    title(['Correlation: ',num2str(R(2,1))])
    saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'TracesFluctuations',filesep,'Fluct-OffsetVsInstData.tif'])



    figure(6)
    plot(OffsetFluct,DataSplineFluct,'.k')
    xlabel('Offset fluctuation')
    ylabel('Fluctuations with spline offset subtraction')
    axis square
    xlim([-4500,4500])
    ylim([-4500,4500])
    R = corrcoef(OffsetFluct,DataSplineFluct)
    title(['Correlation: ',num2str(R(2,1))])
    saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'TracesFluctuations',filesep,'Fluct-OffsetVsSplineData.tif'])
end



%Look at the offset of each particle. Do they all look the same? This is
%only for AP for now
if strcmp(ExperimentAxis,'AP')
    figure(7)
    subplot(1,2,1)
    MaxAP=0;
    MinAP=inf;
    hold all
    for i=1:length(CompiledParticles)
        if sum(CompiledParticles(i).Frame==MaxFrame(end-1))
            MaxAP=max([CompiledParticles(i).MeanAP,MaxAP]);
            MinAP=min([CompiledParticles(i).MeanAP,MinAP]);
            FramePos=find(CompiledParticles(i).Frame==MaxFrame(end-1));
            plot(CompiledParticles(i).MeanAP,CompiledParticles(i).Off(FramePos),'.k')
        end
    end
    hold off
    if MinAP<MaxAP
        xlim([MinAP*0.8,MaxAP*1.2])
    end
    xlabel('AP position')
    ylabel('Offset fluorescence')
    title('Offset at maximum in nc13')
    axis square

    subplot(1,2,2)
    MaxAP=0;
    MinAP=inf;
    hold all
    for i=1:length(CompiledParticles)
        if sum(CompiledParticles(i).Frame==MaxFrame(end))
            MaxAP=max([CompiledParticles(i).MeanAP,MaxAP]);
            MinAP=min([CompiledParticles(i).MeanAP,MinAP]);
            FramePos=find(CompiledParticles(i).Frame==MaxFrame(end));
            plot(CompiledParticles(i).MeanAP,CompiledParticles(i).Off(FramePos),'.k')
        end
    end
    hold off
    xlim([MinAP*0.8,MaxAP*1.2])
    xlabel('AP position')
    ylabel('Offset fluorescence')
    title('Offset at maximum in nc14')
    axis square
    saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'Offset',filesep,'OffsetVsAP.tif'])
end


%Average over all time points

%Average and SD over each time point. In order to do this we'll generate a
%cell array with all the values for a given time point

OffsetCell=cell(length(FrameInfo),1);


for i=1:length(CompiledParticles)
    for j=1:length(CompiledParticles(i).Frame)
        OffsetCell{CompiledParticles(i).Frame(j)}=[OffsetCell{CompiledParticles(i).Frame(j)},...
            CompiledParticles(i).Off(j)];
    end
end


MeanOffsetTrace=cellfun(@mean,OffsetCell,'UniformOutput',false);
SDOffsetTrace=cellfun(@std,OffsetCell,'UniformOutput',false);
NParticlesOffsetTrace=cellfun(@length,OffsetCell,'UniformOutput',false);


MeanOffsetVector=[MeanOffsetTrace{:}];
SDOffsetVector=[SDOffsetTrace{:}];
NOffsetParticles=[NParticlesOffsetTrace{:}];


if strcmp(ExperimentAxis,'AP')
    figure(8)
    IntArea=109;
    errorbar(1:length(MeanOffsetVector),MeanOffsetVector*IntArea,...
        SDOffsetVector*IntArea,'.-r')
    hold on
    errorbar(1:length(MeanVectorAll),MeanVectorAll,...
        SDVectorAll,'.-k')
    hold off
    xlabel('Frame')
    ylabel('Fluorescence (au)')
    xlim([0,length(MeanOffsetVector)*1.1])
    legend('Offset','Spots','Location','Best')
    saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'Offset',filesep,'OffsetAndFluoTime.tif'])


    figure(9)
    errorbar(1:length(MeanOffsetVector),MeanOffsetVector*IntArea-min(MeanOffsetVector*IntArea)+...
        min(MeanVectorAll),...
        SDOffsetVector*IntArea,'.-r')
    hold on
    errorbar(1:length(MeanVectorAll),MeanVectorAll,...
        SDVectorAll,'.-k')
    hold off
    xlabel('Frame')
    ylabel('Fluorescence (au)')
    xlim([0,length(MeanOffsetVector)*1.1])
    legend('Offset (displaced)','Spots','Location','Best')
    axis square
    saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'Offset',filesep,'OffsetAndFluoTime-Displaced.tif'])
end



%% Rate of mRNA production

%Plot the results from fitting the individual traces. Notice that I'm
%redoing the fits here using the range given by FrameRange. This is because
%the fluorescence values were calculated slightly differently sucht hat the
%offset could vary.

if isfield(CompiledParticles,'Fit')

    %First, find the maximum in intensity over traces for each cycle
    nc13Max=max([CompiledParticles(ncFilter(:,end-1)).Fluo]);
    nc14Max=max([CompiledParticles(ncFilter(:,end)).Fluo]);

    figure(10)
    clf

    for i=1:length(CompiledParticles)
        if ~isempty(CompiledParticles(i).Fit)

            %Redo the fit and obtain the parameters

            FrameRange=CompiledParticles(i).Fit.FrameRange;

            FramesRangeFilter=ismember(CompiledParticles(i).Frame,[FrameRange(1):FrameRange(2)]);

            [a, b, sigma_a, sigma_b] = york_fit(ElapsedTime(CompiledParticles(i).Frame(FramesRangeFilter)),...
                        CompiledParticles(i).Fluo(FramesRangeFilter),...
                        ones(1,sum(FramesRangeFilter))*mean(diff(ElapsedTime))/2,...
                        ones(1,sum(FramesRangeFilter))*CompiledParticles(i).FluoError);


            CompiledParticles(i).Fit.Intercept=a;
            CompiledParticles(i).Fit.SDIntercept=sigma_a;

            CompiledParticles(i).Fit.Slope=b;
            CompiledParticles(i).Fit.SDSlope=sigma_b;



            if ncFilter(i,end-1)
                StartFrame=nc13;
                EndFrame=nc14;
                MaxFluo=nc13Max;
            elseif ncFilter(i,end)
                StartFrame=nc14;
                EndFrame=length(ElapsedTime);
                MaxFluo=nc14Max;
            end




            if ~SkipFits



                xRange=linspace(ElapsedTime(FrameRange(1)),...
                    ElapsedTime(FrameRange(end)));


                plot(xRange,b*xRange+a,'-k','LineWidth',3)

                hold on
                errorbar(ElapsedTime(CompiledParticles(i).Frame),CompiledParticles(i).Fluo,...
                    ones(1,length(CompiledParticles(i).Frame))*CompiledParticles(i).FluoError,'.-r')
                hold off
                xlabel('Time (min)')
                ylabel('Fluorescence (au)')
                title(['Compiled particle ',num2str(i),', nc',num2str(CompiledParticles(i).nc)])

                xlim([ElapsedTime(StartFrame),ElapsedTime(EndFrame)])
                ylim([0,MaxFluo])

                saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'Fits',filesep,'Fit',iIndex(i,3),'-nc',...
                    num2str(CompiledParticles(i).nc),'.tif'])
            end
        end
    end

    close(10)
end
     

%% First frames

%How does the first frame from schnitzcells compare to the general one set
%by just looking at the movie? Do this only for nuclei where we have
%identified the parent nucleus.

if HistoneChannel

    figure(11)
    clf
    xRange=linspace(13.5,14.5);
    plot(xRange,ones(size(xRange))*nc14,'-k')
    hold on
    for i=1:length(CompiledParticles)
        if (CompiledParticles(i).nc==14)&(CompiledParticles(i).PParticle>0)
            plot(14*(1+(rand(1)-0.5)/100),CompiledParticles(i).NucStart,'.k')
        end
    end

    xRange=linspace(12.5,13.5);
    plot(xRange,ones(size(xRange))*nc13,'-k')
    hold on
    for i=1:length(CompiledParticles)
        if (CompiledParticles(i).nc==13)&(CompiledParticles(i).PParticle>0)
            plot(13*(1+(rand(1)-0.5)/100),CompiledParticles(i).NucStart,'.k')
        end
    end
    hold off
    ylim([nc13-5,nc14+5])
    set(gca,'XTick',[13,14])
    xlabel('nc')
    ylabel('Frame')
    title('Division time set by eye vs. actual division times of nuclei')
    saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'Various',filesep,'DivisionTimes.tif'])



    %Is there any correlation between the first frame and the time of division?


    figure(12)
    subplot(1,2,1)
    hold on
    for i=1:length(CompiledParticles)
        if (CompiledParticles(i).nc==13)&(CompiledParticles(i).PParticle>0)
            plot(CompiledParticles(i).NucStart*(1+(rand(1)-0.5)/100),CompiledParticles(i).FirstFrame,'.k')
        end
    end
    hold off
    axis square
    xlabel('Nuclear birth (frame)')
    ylabel('First particle frame')
    title('nc13')


    subplot(1,2,2)
    hold on
    for i=1:length(CompiledParticles)
        if (CompiledParticles(i).nc==14)&(CompiledParticles(i).PParticle>0)
            plot(CompiledParticles(i).NucStart*(1+(rand(1)-0.5)/100),CompiledParticles(i).FirstFrame,'.k')
        end
    end
    hold off
    axis square
    title('nc14')
    xlabel('Nuclear birth (frame)')
    ylabel('First particle frame')
    saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'Various',filesep,'DivisionVsFirstFrame.tif'])





    %First frame and AP position
    if strcmp(ExperimentAxis,'AP')
        figure(13)
        clf
        hold on
        for i=1:length(CompiledParticles)
            plot(CompiledParticles(i).MeanAP,....
                ElapsedTime(CompiledParticles(i).FirstFrame)-...
                ElapsedTime(nc14),'.k')
        end
        hold off
        box on
        xlabel('AP position (x/L)')
        ylabel('Particle first frame (min)')
        if length(ElapsedTime) > nc14+20
            ylim([0,ElapsedTime(nc14+20)-ElapsedTime(nc14)])
        else
            ylim([0, ElapsedTime(end) - ElapsedTime(nc14)])
        end
        % ES 2014-01-05 Testing early-nc14-only movies
        saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'Various',filesep,'FirstFrameVsAP.tif'])
    end
end


%% AP position of particle vs nucleus

if HistoneChannel&strcmp(ExperimentAxis,'AP')

    %How different are the AP positions of the nuclei to the particles as a
    %function of time? Let's save the information about nuclear position in the
    %CompiledParticles structure.

    for i=1:length(CompiledParticles)
        for j=1:length(CompiledParticles(i).Frame)
            CurrentNucleus=CompiledParticles(i).Nucleus;
            CurrentFrame=CompiledParticles(i).Frame(j);

            CurrentEllipse=schnitzcells(CurrentNucleus).cellno(...
                find((schnitzcells(CurrentNucleus).frames)==CurrentFrame));

            if ~isempty(CurrentEllipse)
                CompiledParticles(i).NuclearAP(j)=EllipsePos{CurrentFrame}(CurrentEllipse);
            else
                CompiledParticles(i).NuclearAP(j)=nan;
            end
        end
    end

    figure(14)
    clf
    hold all
    for i=1:length(CompiledParticles)
        plot(CompiledParticles(i).APpos-CompiledParticles(i).NuclearAP)
    end
    hold off
    box on
    xlabel('Frame')
    ylabel('AP difference between spot and nucleus (x/L)')
    saveas(gca,[DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'APNucVsParticle.tif'])

end



%% Probability of being on

%I'm going to measure the probability of a nucleus having detectable
%expressiona as a function of time and AP. In order to do this I'll use
%Particles that have both the Approved flag set to 1 and 2. However, I'll
%also check that the nuclei are not too close to the edges.

%NOTE: I need a way to go back and check the nuclei that weren't on. Maybe
%I should move this to Check particles


%Create an image that is partitioned according to the AP bins. We will use
%this to calculate the area per AP bin.

if HistoneChannel&strcmp(ExperimentAxis,'AP')

    %I'll use the Ellipses structure to count nuclei. This is because
    %schnitzcells sometimes misses things at the edges.


    %First, add the corresponding Ellipse number to each frame of the
    %particles. Also save the information in a cell array. This should make
    %searching easier.
    for i=1:length(Particles)
        if ~isempty(Particles(i).Nucleus)
            ParticleNuclei{i}=...
                schnitzcells(Particles(i).Nucleus).cellno(ismember(schnitzcells(Particles(i).Nucleus).frames,...
                Particles(i).Frame));
            ParticleFrames{i}=...
                schnitzcells(Particles(i).Nucleus).frames(ismember(schnitzcells(Particles(i).Nucleus).frames,...
                Particles(i).Frame));
        else
            ParticleNuclei{i}=[];
            ParticleFrames{i}=[];
        end
    end
    
    
    EdgeWidth=10;
    %For each frame find the number of ellipses that are outside of an area
    %delimited from the edge of the image.
    %The information in Ellipses is
    %(x, y, a, b, theta, maxcontourvalue, time, particle_id)
    
    %Initialize matrices where we will store the information of number of
    %particles vs. AP vs. time
    NEllipsesAP=zeros(length(Ellipses),length(APbinID));
    NParticlesEllipsesAP=zeros(length(Ellipses),length(APbinID));
    for i=1:length(Ellipses)
        CurrentEllipses=Ellipses{i};

        Radius=max(CurrentEllipses(:,3:4)')';

        EllipseFilter=(CurrentEllipses(:,1)-Radius-EdgeWidth>0)&...
            (CurrentEllipses(:,2)-Radius-EdgeWidth>0)&...
            (CurrentEllipses(:,1)+Radius+EdgeWidth<FrameInfo(1).PixelsPerLine)&...
            (CurrentEllipses(:,2)+Radius+EdgeWidth<FrameInfo(1).LinesPerFrame);
        

        %Figure out which particles are in this frame and have an approved
        %flag of 1 or 2. Note that we haven't yet checked if their
        %corresponding Ellipses have been approved.
        CurrentParticlesIndex=cellfun(@(x) find(x==i),ParticleFrames,...
            'UniformOutput',false);
        CurrentParticlesFilter=~cellfun(@isempty,CurrentParticlesIndex);
        ParticlesToCheck=find(CurrentParticlesFilter);
        
        
        %Find which of the particles in this frame are related to filtered
        %ellipses and save their corresonding AP information
        FilteredParticlesPos=[];
        for j=1:length(ParticlesToCheck)
            if EllipseFilter(ParticleNuclei{ParticlesToCheck(j)}(CurrentParticlesIndex{ParticlesToCheck(j)}))&...
                    (Particles(ParticlesToCheck(j)).Approved==1 | Particles(ParticlesToCheck(j)).Approved==2)
                FilteredParticlesPos=[FilteredParticlesPos,...
                    EllipsePos{i}(ParticleNuclei{ParticlesToCheck(j)}(CurrentParticlesIndex{ParticlesToCheck(j)}))];
            end
        end
        
        %Count the number of filtered ellipses per AP bin
        EllipsesFilteredPos{i}=EllipsePos{i}(EllipseFilter);
        for j=1:length(EllipsesFilteredPos{i})
            NEllipsesAP(i,max(find(APbinID<=EllipsesFilteredPos{i}(j))))=...
                NEllipsesAP(i,max(find(APbinID<=EllipsesFilteredPos{i}(j))))+1;
        end
 
        
        
        %Count the number of filtered particles per AP bin.
        for j=1:length(FilteredParticlesPos)
            NParticlesEllipsesAP(i,max(find(APbinID<=FilteredParticlesPos(j))))=...
                NParticlesEllipsesAP(i,max(find(APbinID<=FilteredParticlesPos(j))))+1;
        end
        
        
        EllipsesFiltered{i}=Ellipses{i}(EllipseFilter,:);

        NEllipsesFiltered(i)=sum(EllipseFilter);

    end


    %Minimum number of nuclei to actually do the calculation
    MinNuclei=3;
    MinAPIndexProb=min(find(sum(NEllipsesAP>=MinNuclei)));
    MaxAPIndexProb=max(find(sum(NEllipsesAP>=MinNuclei)));

    MinNucleiFilter=NEllipsesAP>=MinNuclei;

    %Calculate the ratio
    OnRatioAP=NParticlesEllipsesAP./NEllipsesAP;
    %Filter out the elements that correspond to a number of nuclei below our
    %limit of MinNuclei
    OnRatioAP(~MinNucleiFilter)=nan;
    OnRatioAP=reshape(OnRatioAP,size(MinNucleiFilter));


    if MaxAPIndexProb>MinAPIndexProb
        colormap(jet(128));
        cmap=colormap;

        Color=cmap(round((APbinID(MinAPIndexProb:MaxAPIndexProb)-...
            APbinID(MinAPIndexProb))/...
            (APbinID(MaxAPIndexProb)-APbinID(MinAPIndexProb))*127)+1,:);
        figure(15)
        clf
        PlotHandle=[];
        hold on
        for j=MinAPIndexProb:MaxAPIndexProb
            PlotHandle=[PlotHandle,...
                plot(ElapsedTime,OnRatioAP(:,j),'color',Color(j-MinAPIndexProb+1,:))];
        end
        hold off
        xlabel('Time (min)')
        ylabel('Fraction of on nuclei')
        h = colorbar;
        caxis([APbinID(MinAPIndexProb),APbinID(MaxAPIndexProb)])
        ylabel(h,'AP Position (x/L)')
        StandardFigure(PlotHandle,gca)
        xlim([0,ElapsedTime(end)])
        ylim([0,1.01])
        saveas(gca,[DropboxFolder,filesep,Prefix,filesep,'Probabilities',filesep,'ProbVsTimeVsAP.tif'])
    end




    %Now I want to compute the probability of nuclei being on in at least one
    %frame over the whole nc. This is a ltitle bit tricky because I haven't
    %checked the tracking of the nuclei without particles. As a result, those
    %lineages are not complete and we could possibly overcount off nuclei if we
    %just looked at the schnitzcells we have right now. I'll fix this later,
    %but for now I'm going calculate the number of on nuclei per AP bin, where
    %I'll actually check how much area corresponds to that AP bin.





    %What's the probability of a nucleus being on for at least one frame in an nc?
    %This will be an array with rows for nc13 and nc14 and columns
    %corresponding to each AP bin.
    %This only works if we trust the tracking within one nc
    ParticleCountAP=zeros(3,length(APbinID));

    try
        for i=1:length(Particles)
            %See if the particle has either the flag 1 or 2
            if (Particles(i).Approved==1)|(Particles(i).Approved==2)

                %Determine the nc so we can add to the right position of ParticleCountAP
                if (FrameInfo(min(Particles(i).Frame(Particles(i).FrameApproved))).nc)>=12
                    CurrentNC=FrameInfo(min(Particles(i).Frame(Particles(i).FrameApproved))).nc;

                    %Now determine the AP bin this particle is on. We'll use
                    %the nucleus positioning. Also, in order to count it we
                    %need to make sure that it fulfills the same criteria we
                    %used to count NEllipsesAP
                    EllipseAPPosTemp=[];
                    EllipsePosXTemp=[];
                    EllipsePosYTemp=[];
                    RadiusTemp=[];
                    xPosTemp=[];
                    yPosTemp=[];
                    for j=1:length(schnitzcells(Particles(i).Nucleus).frames)
                        CurrentEllipse=Ellipses{schnitzcells(Particles(i).Nucleus).frames(j)}(schnitzcells(Particles(i).Nucleus).cellno(j),:);
                        RadiusTemp=[RadiusTemp,max(CurrentEllipse(3:4))];

                        %Determine the AP position
                        EllipseAPPosTemp=[EllipseAPPosTemp,...
                            EllipsePos{schnitzcells(Particles(i).Nucleus).frames(j)}(schnitzcells(Particles(i).Nucleus).cellno(j))];

                        %Determine the x and y positions
                        xPosTemp=[xPosTemp,CurrentEllipse(1)];
                        yPosTemp=[yPosTemp,CurrentEllipse(2)];
                    end

                    %Make sure that this nucleus was inside the limits at all
                    %points

                    if sum((xPosTemp-RadiusTemp-EdgeWidth>0)&...
                        (yPosTemp-RadiusTemp-EdgeWidth>0)&...
                        (xPosTemp+RadiusTemp+EdgeWidth<FrameInfo(1).PixelsPerLine)&...
                        (yPosTemp+RadiusTemp+EdgeWidth<FrameInfo(1).LinesPerFrame))==...
                        length(xPosTemp)

                        MeanAP=mean(EllipseAPPosTemp);

                        ParticleCountAP(CurrentNC-11,max(find(APbinID<=MeanAP)))=...
                            ParticleCountAP(CurrentNC-11,max(find(APbinID<=MeanAP)))+1;
                    end
                end
            end
        end
    end

    %Calculate the probability using the mean number of nuclei per AP bin
    %in each nc. Note that we look at a reduced range within the nc to
    %reduce variability in counting at mitosis.
    ParticleCountProbAP(:,1)=ParticleCountAP(1,:)./mean(NEllipsesAP(nc12+5:nc13-5,:));
    ParticleCountProbAP(:,2)=ParticleCountAP(2,:)./mean(NEllipsesAP(nc13+5:nc14-5,:));
    ParticleCountProbAP(:,3)=ParticleCountAP(3,:)./...
        mean(NEllipsesAP(max(1,nc14-5):length(FrameInfo)-5,:));
    % ES 2014-01-08: accounting for movies started fewer than 5 frames before
    % mitosis 13
    
    figure(16)   
    plot(APbinID,ParticleCountProbAP(:,1),'.-b')
    hold on
    plot(APbinID,ParticleCountProbAP(:,2),'.-k')
    plot(APbinID,ParticleCountProbAP(:,3),'.-r')
    hold off
    %ylim([0,max(ParticleCountAP(1,:))*2*1.1])
    xlabel('AP position (x/L)')
    ylabel('Number of on nuclei per unit AP position')
    saveas(gca,[DropboxFolder,filesep,Prefix,filesep,'Probabilities',filesep,'ProbVsAP.tif'])

    
    %Use the alternative approach I used for the movies. We are oging to
    %look at each nucleus towards the end of each nc and ask if they
    %correspond to an on or off particle in any frame previous to that one.
    
    
    %We'll go for 2.5 minutes before the next mitosis. I might relate this
    %later to the elongation time as a way to say that these particles
    %won't contribute to the total amount of mRNA produced anyway.
    FramesBack=ceil(2.5/mean(diff(ElapsedTime)));
    
    TotalEllipsesAP=zeros(length(APbinID),3);
    EllipsesOnAP=zeros(length(APbinID),3);
    for nc=12:14
        
        %Figure out which frame we'll look at
        if nc==14
            FrameToUse=length(FrameInfo)-FramesBack;
        else
            FrameToUse=eval(['nc',num2str(nc+1)])-FramesBack;
        end
        
        if FrameToUse>0
            %Filter ellipses that are within the image
            CurrentEllipses=Ellipses{FrameToUse};

            Radius=max(CurrentEllipses(:,3:4)')';

            EllipseFilter=(CurrentEllipses(:,1)-Radius-EdgeWidth>0)&...
                (CurrentEllipses(:,2)-Radius-EdgeWidth>0)&...
                (CurrentEllipses(:,1)+Radius+EdgeWidth<FrameInfo(1).PixelsPerLine)&...
                (CurrentEllipses(:,2)+Radius+EdgeWidth<FrameInfo(1).LinesPerFrame);

            %Check if the filtered ellipses had an associated particle
            EllipsesToCheck=find(EllipseFilter);




            for j=1:length(EllipsesToCheck)
                %Find which AP bind we're in
                CurrentAPbin=max(find(APbinID<EllipsePos{FrameToUse}(EllipsesToCheck(j))));
                %Count the total amount of ellipses in the right AP bin
                TotalEllipsesAP(CurrentAPbin,nc-11)=TotalEllipsesAP(CurrentAPbin,nc-11)+1;



                %Find the schnitz this corresponds to
                for k=1:length(schnitzcells)

                    IndexToUse=find((schnitzcells(k).frames)==FrameToUse);
                    if ~isempty(IndexToUse)

                        %Check this schnitz for consistency with cellno.
                        %Otherwise fix it. I obtained the fixing code from
                        %TrackmRNADynamicsV2.m
                        if length(schnitzcells(k).frames)~=length(schnitzcells(k).cellno)
                           %If there number of frames is different from the number of
                           %cellno then use the cenx and ceny to find the cellno in
                           %Ellipses an repopulate this schnitz
                           if (length(schnitzcells(k).frames)==length(schnitzcells(k).cenx))&...
                                   (length(schnitzcells(k).frames)==length(schnitzcells(k).ceny))
                               for m=1:length(schnitzcells(k).frames)
                                    %The information in Ellipses is
                                    %(x, y, a, b, theta, maxcontourvalue, time, particle_id)
                                    MaxDistance=2;  %Maximum pixel distance to identify an
                                                    %ellipse with a schnitz
                                    Distances=sqrt((Ellipses{schnitzcells(k).frames(m)}(:,1)-...
                                       schnitzcells(k).cenx(m)).^2+...
                                       (Ellipses{schnitzcells(k).frames(m)}(:,2)-...
                                       schnitzcells(k).ceny(m)).^2);
                                    [MinValue,MinIndex]=min(Distances);

                                   %Make sure no other schnitz is associated to this
                                   %ellipse
                                   EllipseFoundElsewhere=0;
                                   for n=[1:k-1,k+1:length(schnitzcells)]
                                       %Only consider it if the schnitzcell is also valid!
                                       if (length(schnitzcells(n).frames)==length(schnitzcells(n).cellno))
                                           if sum(schnitzcells(n).frames==schnitzcells(k).frames(m))
                                              IndexToCheck=find(schnitzcells(n).frames==schnitzcells(k).frames(m));

                                              %The schnitz I'm comparing to
                                              %might also be screwed up.
                                              %I'd have to compare its cenx
                                              %and ceny to be sure
                                              try
                                                  if schnitzcells(k).cellno(IndexToCheck)==MinIndex
                                                        error('duplicated schnitz?')
                                                        DistancesK=sqrt((Ellipses{schnitzcells(k).frames(m)}(:,1)-...
                                                           schnitzcells(n).cenx(IndexToCheck)).^2+...
                                                           (Ellipses{schnitzcells(k).frames(m)}(:,2)-...
                                                           schnitzcells(n).ceny(IndexToCheck)).^2);          

                                                        [MinValueK,MinIndexK]=min(DistancesK);
                                                        if MinValue<MinValueK
                                                            schnitzcells(n).cellno(IndexToCheck)=[];                                        
                                                        end
                                                    end
                                              end
                                           end
                                       end
                                    end

                                    if ~EllipseFoundElsewhere
                                       schnitzcells(k).cellno(m)=MinIndex;
                                    else
                                       MinValue
                                       1+1; error('What to do here?')
                                    end
                               end

                           else
                               error('Cannnot rescue schnitz')
                           end
                        end

                    end

                    if schnitzcells(k).cellno(IndexToUse)==EllipsesToCheck(j)
                        %Now see if there is an associated particle with it
                        for m=1:length(CompiledParticles)
                            if CompiledParticles(m).Nucleus==k
                                [j,k,m];
                                EllipsesOnAP(CurrentAPbin,nc-11)=EllipsesOnAP(CurrentAPbin,nc-11)+1;
                            end                        
                        end
                    end
               end
            end
        end
    end
        
    figure(17)
	plot(APbinID,EllipsesOnAP(:,1)./TotalEllipsesAP(:,1),'.-b')
    hold on
    plot(APbinID,EllipsesOnAP(:,2)./TotalEllipsesAP(:,2),'.-k')
    plot(APbinID,EllipsesOnAP(:,3)./TotalEllipsesAP(:,3),'.-r')
    hold off
     
end




%% Movie of AP profile

%I want to make a movie of the average fluorescence as a function of AP as
%a function of time. In order to make life easier I'll just export to a
%folder. I can then load everything in ImageJ.

if ~SkipMovie&strcmp(ExperimentAxis,'AP')
    figure(17)

    MaxValue=max(max(MeanVectorAP));

    MinParticles=4;     %Minimum number of particles in a bin to take it seriously
    NParticlesAPFilter=NParticlesAP>=MinParticles;

    for i=1:length(FrameInfo)
        PlotHandle=errorbar(APbinID(NParticlesAPFilter(i,:)),...
            MeanVectorAP(i,NParticlesAPFilter(i,:)),SDVectorAP(i,NParticlesAPFilter(i,:)),'.-k');
        hold on
        PlotHandle=[PlotHandle,errorbar(APbinID(NParticlesAPFilter(i,:)),...
            MeanVectorAP(i,NParticlesAPFilter(i,:)),...
            SDVectorAP(i,NParticlesAPFilter(i,:))./sqrt(NParticlesAP(i,NParticlesAPFilter(i,:))),'-k')];
        hold off
        xlim([0.1,0.8])
        ylim([0,MaxValue])
        xlabel('AP position (x/L)')
        ylabel('Mean fluorescence')

        if exist(['nc',num2str(FrameInfo(i).nc)])
            if eval(['nc',num2str(FrameInfo(i).nc)])>0
                title(['nc',num2str(FrameInfo(i).nc),'. Time into nc: ',num2str(round((ElapsedTime(i)-...
                    ElapsedTime(eval(['nc',num2str(FrameInfo(i).nc)])))*10)/10),' min. Total time: ',...
                    num2str(round(ElapsedTime(i)*10)/10),' min (Frame ',num2str(i),').'])
            else
                title(['nc',num2str(FrameInfo(i).nc),'. Total time: ',...
                    num2str(round(ElapsedTime(i)*10)/10),' min (Frame ',num2str(i),').'])
            end
        else
            title(['nc',num2str(FrameInfo(i).nc),'. Total time: ',...
                num2str(round(ElapsedTime(i)*10)/10),' min (Frame ',num2str(i),').'])
        end

        StandardFigure(PlotHandle,gca)

        saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'APMovie',filesep,iIndex(i,3),'.tif']);   
    end
    
end








%% Checking correlations with position and expression
% 
% 1+1;
% 
% i=100
% 
% for j=1:length(CompiledParticles(i).Frame)
% 
%     CurrentFrame=CompiledParticles(i).Frame(j);
%     EllipsePos=find((schnitzcells(CompiledParticles(i).Nucleus).frames)==CurrentFrame);
%     CurrentEllipse=Ellipses{CompiledParticles(i).Frame};
%     CurrentEllipse=CurrentEllipse(EllipsePos,:);
% 
%     Position(j)=(CompiledParticles(i).xPos(j)-CurrentEllipse(1))^2+...
%         (CompiledParticles(i).yPos(j)-CurrentEllipse(2))^2;
% end
% 
% figure(9)
% plot(  Position ,CompiledParticles(i).Fluo,'.k')
%     
%     
    





%% Save everything
if HistoneChannel&strcmp(ExperimentAxis,'AP')
    save([DropboxFolder,filesep,Prefix,filesep,'CompiledParticles.mat'],...
        'CompiledParticles','ElapsedTime','NewCyclePos','nc9','nc10','nc11',...
        'nc12','nc13','nc14','ncFilterID','StemLoopEnd','ncFilter','APbinID','APFilter',...
        'MeanVectorAP','SDVectorAP','NParticlesAP','MeanVectorAll',...
        'SDVectorAll','NParticlesAll','MaxFrame','MinAPIndex','MaxAPIndex',...
        'AllTracesVector','AllTracesAP','MeanCyto','SDCyto','MedianCyto','MaxCyto',...
        'MeanOffsetVector','SDOffsetVector','NOffsetParticles',...
        'MeanSlopeVectorAP','SDSlopeVectorAP','NSlopeAP',...
        'ParticleCountAP','APbinArea','OnRatioAP','NEllipsesAP',...
        'ParticleCountProbAP',...
        'EllipsesOnAP','TotalEllipsesAP',...
        'EllipsePos','EllipsesFilteredPos','FilteredParticlesPos')
elseif HistoneChannel&strcmp(ExperimentAxis,'DV')
    save([DropboxFolder,filesep,Prefix,filesep,'CompiledParticles.mat'],...
        'CompiledParticles','ElapsedTime','NewCyclePos','nc9','nc10','nc11',...
        'nc12','nc13','nc14','StemLoopEnd','ncFilterID','ncFilter',...
        'MeanVectorAll',...
        'SDVectorAll','NParticlesAll','MaxFrame',...
        'AllTracesVector','MeanCyto','SDCyto','MedianCyto','MaxCyto',...
        'MeanOffsetVector','SDOffsetVector','NOffsetParticles')
else
    save([DropboxFolder,filesep,Prefix,filesep,'CompiledParticles.mat'],...
        'CompiledParticles','ElapsedTime','NewCyclePos','nc9','nc10','nc11',...
        'nc12','nc13','nc14','StemLoopEnd','ncFilterID','ncFilter','APbinID','APFilter',...
        'MeanVectorAP','SDVectorAP','NParticlesAP','MeanVectorAll',...
        'SDVectorAll','NParticlesAll','MaxFrame','MinAPIndex','MaxAPIndex',...
        'AllTracesVector','AllTracesAP','MeanCyto','SDCyto','MedianCyto','MaxCyto',...
        'MeanOffsetVector','SDOffsetVector','NOffsetParticles',...
        'MeanSlopeVectorAP','SDSlopeVectorAP','NSlopeAP')
end

