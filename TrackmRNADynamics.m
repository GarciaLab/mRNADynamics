function [Particles,schnitzcells]=TrackmRNADynamics(varargin)

%This function sets up particle tracking using the FISH analysis code. If
%nuclei have been tracked it uses that information for the particle
%tracking.

%To do: The no-histone part of the code doesn't take into account the
%Approved field of the Particles structure.


[SourcePath,FISHPath,DefaultDropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders;


%Look at the input parameter and use defaults if missing
if isempty(varargin)
    %Folders
    FolderTemp=uigetdir(DefaultDropboxFolder,'Select folder with data to analyze');
    Dashes=strfind(FolderTemp,filesep);
    Prefix=FolderTemp((Dashes(end)+1):end);
    
    %Thresholds
    Threshold1=100;     %This is the first threshold we apply. We'll use it for
                        %the initial search and then switch to Threshold2.
    Threshold2=30;
elseif ~ischar(varargin{1})
    %Folders
    FolderTemp=uigetdir(DefaultDropboxFolder,'Select folder with data to analyze');
    Dashes=strfind(FolderTemp,filesep);
    Prefix=FolderTemp((Dashes(end)+1):end);
    
    %Thresholds
    Threshold1=varargin{1};
    Threshold2=varargin{2};    
else
    %Folders
    Prefix=varargin{1};
    
    %Thresholds
    if length(varargin)==3
        Threshold1=varargin{2};
        Threshold2=varargin{3}; 
    elseif length(varargin)==5
        Threshold1(1)=varargin{2};
        Threshold2(1)=varargin{3};
        Threshold1(2)=varargin{4};
        Threshold2(2)=varargin{5};
    else
        error('Number of parameters not recognized.')
    end

end

%Save a backup of the threshold. We'll check that it didn't change if we're
%doing retracking. This is definitely not an elegant solution
Threshold1Backup=Threshold1;
Threshold2Backup=Threshold2;


%Get the actual folder now that we have the Prefix
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders(Prefix);

[~,~,~,~,~,~,~,ExperimentType, Channel1, Channel2,~] =...
    readMovieDatabase(Prefix);

%What type of experiment are we dealing with? Get this out of
%MovieDatabase.xlsx
[XLSNum,XLSTxt]=xlsread([DefaultDropboxFolder,filesep,'MovieDatabase.xlsx']);
ExperimentTypeColumn=find(strcmp(XLSTxt(1,:),'ExperimentType'));
DataFolderColumn=find(strcmp(XLSTxt(1,:),'DataFolder'));
Channel1Column=find(strcmp(XLSTxt(1,:),'Channel1'));
Channel2Column=find(strcmp(XLSTxt(1,:),'Channel2'));


Dashes=findstr(Prefix,'-');
PrefixRow=find(strcmp(XLSTxt(:,DataFolderColumn),[Prefix(1:Dashes(3)-1),'\',Prefix(Dashes(3)+1:end)]));
if isempty(PrefixRow)
    PrefixRow=find(strcmp(XLSTxt(:,DataFolderColumn),[Prefix(1:Dashes(3)-1),'/',Prefix(Dashes(3)+1:end)]));
    if isempty(PrefixRow)
        error('Could not find data set in MovieDatabase.XLSX. Check if it is defined there.')
    end
end

ExperimentType=XLSTxt(PrefixRow,ExperimentTypeColumn);
Channel1=XLSTxt(PrefixRow,Channel1Column);
Channel2=XLSTxt(PrefixRow,Channel2Column);

%Set the source folders
Folder=[FISHPath,filesep,Prefix,'_',filesep,'preanalysis',filesep];
FileName=['CompactResults_',Prefix,'_.mat'];

%Set the destination folders
OutputFolder=[DropboxFolder,filesep,Prefix];

FilePrefix=FileName(16:end-4);
    
DataFolder=[PreProcPath,filesep,FilePrefix(1:end-1)];


%Determine the search radius based on the imaging conditions
SearchRadiusMicrons=3;       %Search radius in um



%Check whether we're using the same threshold in the case of retracking
if exist([OutputFolder,filesep,'Particles.mat'])
    load([OutputFolder,filesep,'Particles.mat'],'Threshold1','Threshold2')
        
    if (~sum(Threshold1==Threshold1Backup)==length(Threshold1))&...
        (~sum(Threshold2==Threshold2Backup)==length(Threshold2))
            Answer=input('Thresholds changed, will delete previous tracking. Proceed? (y/n):','s');
            if strcmp(lower(Answer),'y')
                Threshold1=Threshold1Backup;
                Threshold2=Threshold2Backup;
                delete([OutputFolder,filesep,'Particles.mat'])
            else
                error('Cannot retrack if the threshold changed')
            end
  
    end

end

%Load the information about this image

%Check if we have FrameInfo otherwise try to get the information straight
%from the file.
if exist([OutputFolder,filesep,'FrameInfo.mat'])
    load([OutputFolder,filesep,'FrameInfo.mat'])
    
    %See if this came from the 2-photon, which is the default
    if ~isfield(FrameInfo,'FileMode')
    
        if (FrameInfo(1).ZoomFactor==8)&(FrameInfo(1).PixelsPerLine==256)
            PixelSize=0.22;     %This is the pixel size in um for a zoom of 8.
        elseif (FrameInfo(1).ZoomFactor==4)&...
                (FrameInfo(1).PixelsPerLine==512)&...
                (FrameInfo(1).ScanAmplitudeY==1)
            PixelSize=0.22;
        elseif (FrameInfo(1).ZoomFactor==16)&(FrameInfo(1).PixelsPerLine==128)
            PixelSize=0.22;
        else
            display('Warning: Imaging setting not defined. Using a pixel size of 0.22um')
            PixelSize=0.22;
        end
        
    %For data sets from the two photon with the actual TIF file mode
    elseif strcmp(FrameInfo(end).FileMode,'TIF')
         if (FrameInfo(1).ZoomFactor==8)&(FrameInfo(1).PixelsPerLine==256)
            PixelSize=0.22;     %This is the pixel size in um for a zoom of 8.
        elseif (FrameInfo(1).ZoomFactor==4)&...
                (FrameInfo(1).PixelsPerLine==512)&...
                (FrameInfo(1).ScanAmplitudeY==1)
            PixelSize=0.22;
        elseif (FrameInfo(1).ZoomFactor==16)&(FrameInfo(1).PixelsPerLine==128)
            PixelSize=0.22;
        else
            display('Warning: Imaging setting not defined. Using a pixel size of 0.22um')
            PixelSize=0.22;
        end
        
    elseif strcmp(FrameInfo(1).FileMode,'LSM')|strcmp(FrameInfo(1).FileMode,'LSMExport')
        PixelSize=FrameInfo(1).PixelSize;
    elseif strcmp(FrameInfo(1).FileMode,'LIFExport') || strcmp(FrameInfo(1).FileMode,'LAT') || strcmp(FrameInfo(1).FileMode,'DSPIN')  %CS20170907
        PixelSize=FrameInfo(1).PixelSize;
    end
else
    warning('No FrameInfo.mat detected. Trying to pull out magnification information from the TIF file')

    DZoom=dir([PreProcPath,filesep,Prefix,filesep,'*z*.tif']);
    ImageInfo=imfinfo([PreProcPath,filesep,Prefix,filesep,DZoom(1).name]);
    PixelSize=1/ImageInfo.XResolution;
end
    
SearchRadius=ceil(SearchRadiusMicrons/PixelSize);   


%Check if we have tracked the lineages of the nuclei
if exist([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'])
    UseHistone=1;
    
	%Load the nuclei segmentation information
    load([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'])
    
    %Load the nuclei tracking information
    load([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'])
    
    %Do a bunch of check on schnitzcells only if tracking for the first
    %time.
       
else
    UseHistone=0;
    warning('Warning: No nuclei lineage tracking found. Proceeding with tracking particles only.')
end


%Single color mode and 2spot2color mode
if strcmpi(ExperimentType,'1spot')||strcmpi(ExperimentType,'2spot')||...
        strcmpi(ExperimentType,'inputoutput')||strcmpi(ExperimentType,'2spot2color')

    %Figure out which channel has the spots
    if strcmp(ExperimentType,'1spot')||strcmp(ExperimentType,'2spot')
        SpotsChannel=1;
    elseif strcmp(ExperimentType,'inputoutput')
        SpotsChannel=find(~cellfun(@isempty,strfind(lower([Channel1,Channel2]),'mcp'))|...
            ~cellfun(@isempty,strfind(lower([Channel1,Channel2]),'pcp')));
        if length(SpotsChannel)>1
            error('only one output channel currently supported in inputoutput mode')
        end
    elseif strcmpi(ExperimentType,'2spot2color')
        SpotsChannel=[1,2];
    end
    
    %Get the number of channels
    NCh=length(SpotsChannel);    
    
    %Check if particle tracking has already been done on this dataset
    if exist([OutputFolder,filesep,'Particles.mat'])

        load([OutputFolder,filesep,'Particles.mat'])
        
        %If there's only one channel, Particles, Spots and other structures are
        %not saved as cells. We turn them into a cell to ensure
        %compatibility.
        if NCh==1;
            Particles={Particles};
        end
        
        for Channel=1:NCh
            if isfield(Particles{1},'Approved')
                Retracking{Channel}=1;           %Flag for whether we are performing retracking
                display(['Performing retracking on channel ',num2str(Channel)])

                %Only keep the approved particles and start the tracking from there
                    k=1;
                    for i=1:length(Particles{Channel})
                        if Particles{Channel}(i).Approved~=0
                            NewParticles{Channel}(k)=Particles{Channel}(i);
                            k=k+1;
                        end
                    end
                    if exist('NewParticles')
                        Particles{Channel}=NewParticles{Channel};
                    else
                        Particles{Channel}=[];
                        Retracking{Channel}=0;
                    end

            else
                Retracking{Channel}=0;
                Particles{Channel}=[];
            end
        end
    else
        for Channel=1:NCh
            Particles{Channel}=[];   %This is the structure where we'll be tracking all particles.
            Retracking{Channel}=0;
        end
    end

    %First, generate a structure array with a flag that determines whether
    %each spots is above Threshold1
    if ~exist('Spots')
        load([DropboxFolder,filesep,Prefix,filesep,'Spots.mat'])
        
        %If there's only one channel, Particles, Spots and other structures are
        %not saved as cells. We turn them into a cell to ensure
        %compatibility.
        if NCh==1
            Spots={Spots};
        end
        
        for Channel=1:NCh
            %Determine the maximum number of spots in a given frame for the
            %whole movie
            MaxSpots{Channel}=0;
            for i=1:length(Spots{Channel})
                MaxSpots{Channel}=max([MaxSpots{Channel},...
                    length(Spots{Channel}(i).Fits)]);
            end
            %This filter tells us whether a spot is above the threshold.
            SpotFilter{Channel}=nan(length(Spots{Channel}),MaxSpots{Channel});
            %Populate the filter
            for i=1:length(Spots{Channel})
                for j=1:length(Spots{Channel}(i).Fits)
                    if sum(Spots{Channel}(i).Fits(j).DOGIntensity>Threshold1(Channel))
                        SpotFilter{Channel}(i,j)=1;
                    end
                end
            end
        end
    end

    %Start by numbering the particles found
    ParticlesFig=figure;
    NucleiFig=figure;

    
    %See how  many frames we have and adjust the index size of the files to
    %load accordingly
    if length(FrameInfo)<1E3
        NDigits=3;
    elseif length(FrameInfo)<1E4
        NDigits=4;
    else
        error('No more than 10,000 frames supported. Change this in the code')
    end
    
    for Channel=1:NCh %Iterate over the channels
        %Initially, only track particles that are above Threshold1
        for i=1:length(Spots{Channel}) %iterate over all frames
            %I should probably not have both i and CurrentFrame here
            CurrentFrame=i;


            figure(ParticlesFig)
            set(gcf,'units', 'normalized', 'position',[0.01, .55, .33, .33]);

            %Get the filter for this frame
            CurrentFrameFilter=...
                logical(SpotFilter{Channel}(CurrentFrame,~isnan(SpotFilter{Channel}(CurrentFrame,:))));

            %Get the positions of the spots in this frame
            [x,y,~]=SpotsXYZ(Spots{Channel}(CurrentFrame));

            %Z plane to be displayed. We use the median of all particles found
            %in this frame
            if ~isempty(Spots{Channel}(CurrentFrame).Fits)
                CurrentZ=round(median([Spots{Channel}(CurrentFrame).Fits.brightestZ]));
            else
                CurrentZ=round(FrameInfo(1).NumberSlices/2);     
            end

            %Load the corresponding mRNA image. Check whether we have multiple
            %channels saved or not.
            D=dir([PreProcPath,filesep,Prefix,filesep,FilePrefix,iIndex(CurrentFrame,3),'_z',iIndex(CurrentZ,2),'*.tif']);
            Image=imread([PreProcPath,filesep,Prefix,filesep,FilePrefix,iIndex(CurrentFrame,3),'_z',iIndex(CurrentZ,2),'_ch',iIndex(SpotsChannel(Channel),2),'.tif']);

            %TO-DO: Show spots above and below threshold differently
            imshow(Image,[])
            hold on
            plot(x(CurrentFrameFilter),y(CurrentFrameFilter),'or', 'MarkerSize', 20)
            plot(x(~CurrentFrameFilter),y(~CurrentFrameFilter),'ow', 'MarkerSize', 20)
            hold off
            title(i)
            set(gcf,'Name',['Frame: ',num2str(CurrentFrame),'/',num2str(length(Spots))])

            figure(NucleiFig)
            set(gcf,'units', 'normalized', 'position',[0.65, .5, .2, .2])

            if UseHistone
                try
                    Image=imread([PreProcPath,filesep,Prefix,filesep,FilePrefix(1:end-1),'-His_',iIndex(CurrentFrame,NDigits),'.tif']);
                catch
                    try
                        Image=imread([PreProcPath,filesep,Prefix,filesep,FilePrefix(1:end-1),'_His_',iIndex(CurrentFrame,NDigits),'.tif']);
                    catch
                        Image=0;
                    end
                end
                imshow(Image,[],'Border','Tight')
                hold on
                PlotHandle=[];
                [NEllipses,Dummy]=size(Ellipses{CurrentFrame});
                for j=1:NEllipses
                    PlotHandle=[PlotHandle,ellipse(Ellipses{CurrentFrame}(j,3),...
                        Ellipses{CurrentFrame}(j,4),...
                        Ellipses{CurrentFrame}(j,5),Ellipses{CurrentFrame}(j,1)+1,...
                        Ellipses{CurrentFrame}(j,2)+1)];
                        text(Ellipses{CurrentFrame}(j,1)+1,...
                            Ellipses{CurrentFrame}(j,2)+1,num2str(j),'BackgroundColor',[.7 .9 .7]);

                end
                set(PlotHandle,'Color','r')
                hold off
                title(i)
            end
            drawnow


            %If we don't have nuclear tracking then track the particles based on
            %proximity

            if ~UseHistone
                %Get the particles detected in the frame
                if ~isempty(x)

                    %Find the approved spots in this frame
                    ApprovedSpots=find(SpotFilter{Channel}(CurrentFrame,~isnan(SpotFilter{Channel}(CurrentFrame,:))));

                    %Get the positions of ALL spots (approved and
                    %disapproved)
                    [NewSpotsX,NewSpotsY]=SpotsXYZ(Spots{Channel}(CurrentFrame));
                    %Filter for the approved spots
                    NewSpotsX=NewSpotsX(ApprovedSpots);
                    NewSpotsY=NewSpotsY(ApprovedSpots);

                    if isempty(Particles{Channel})
                        %Initialize the Particles structure if it doesn't exist yet
                        for j=1:length(ApprovedSpots)
                            Particles{Channel}(j).Frame=CurrentFrame;
                            Particles{Channel}(j).Index=ApprovedSpots(j);
                            Particles{Channel}(j).Approved=0;
                        end


                    else
                        %If we already have recorded particles, we need to compare
                        %them to the new ones found and try to assign them.



                        %Get a list of the particles that were present in
                        %the previous frame and of their positions.
                        PreviousFrameParticles=[];
                        xPreviousFrameParticles=[];
                        yPreviousFrameParticles=[];
                        [PreviousSpotsX,PreviousSpotsY]=SpotsXYZ(Spots{Channel}(CurrentFrame-1));
                        for j=1:length(Particles{Channel})
                            if Particles{Channel}(j).Frame(end)==(CurrentFrame-1)
                                PreviousFrameParticles=[PreviousFrameParticles,j];
                                PreviousFrameIndex=Particles{Channel}(j).Index(end);
                                xPreviousFrameParticles=[xPreviousFrameParticles,...
                                    PreviousSpotsX(PreviousFrameIndex)];
                                yPreviousFrameParticles=[yPreviousFrameParticles,...
                                    PreviousSpotsY(PreviousFrameIndex)];
                            end
                        end

                        %If there were particles present in the previous frame,
                        %then we find their distances to the spots present in
                        %the current frame. Otherwise, we create new particles
                        %for each new spot.
                        %We keep track of which spots goes to a new or old
                        %particle using the NewParticle array.
                        NewParticleFlag=true(size(ApprovedSpots));


                        if ~isempty(PreviousFrameParticles)
                            %Get the distances between the spots in this frame
                            %and those within the particles in the previous
                            %frame
                            clear Distance
                            %The rows of Distance correspond to the new spots.
                            %The columns correspond to the particles present in
                            %the previous frame. Each element is the distance.
                            for j=1:length(NewSpotsX)
                                Distance(j,:)=sqrt((NewSpotsX(j)*PixelSize-...
                                    xPreviousFrameParticles*PixelSize).^2+...
                                    (NewSpotsY(j)*PixelSize-...
                                    yPreviousFrameParticles*PixelSize).^2);
                            end

                            %We want to make sure there is only one match of
                            %new spot to previous particle. To make this
                            %possible, we'll go through each column in the
                            %matrix Distance and set all pairwise distances
                            %that are not the minimum one to infinity.
                            for j=1:length(PreviousFrameParticles)
                                DistanceFilter=false(size(Distance(:,j)));
                                [DistMinValue,DistMindIndex]=min(Distance(:,j));
                                DistanceFilter(DistMindIndex)=true;
                                Distance(~DistanceFilter,j)=inf;
                            end

                            %The rows of distance correspond to the new particles.
                            %The columns correspond to their distance to the old
                            %particles.


                            %The followign tracking works well if we have more
                            %than one previous particle. If not, we need to be
                            %more careful.
                            if (size(Distance,2)>1)
                                %MinIndex is a row vector. The element position
                                %correspond to the new spot and the value within it
                                %correspond to the previous particle that it's
                                %closest to.
                                [MinValues,MinIndex]=min(Distance');
                                %Note that inf can be a distance as well. In those
                                %cases, turn MinIndex to 0.
                                MinIndex(MinValues==inf)=0;
                                %Now, check that the distances are smaller than
                                %SearchRadius
                                MinIndex(~(MinValues<SearchRadius))=0;

                                %Assign the new spots to their
                                %corresponding particles.
                                if sum(MinIndex)
                                    for j=1:length(MinIndex)
                                        if MinIndex(j)>0
                                            Particles(PreviousFrameParticles(MinIndex(j))).Frame(end+1)=CurrentFrame;
                                            Particles(PreviousFrameParticles(MinIndex(j))).Index(end+1)=ApprovedSpots(j);

                                            %We don't want this new spot to generate a
                                            %new particle further below
                                            NewParticleFlag(j)=false;
                                        end
                                    end
                                end
                            else
                                %Find the new spot that is closest to the one
                                %previous particle
                                [MinValues,MinIndex]=min(Distance);
                                %Note that inf can be a distance as well. In those
                                %cases, turn MinIndex to 0.
                                MinIndex(MinValues==inf)=0;
                                %Now, check that the distances are smaller than
                                %SearchRadius
                                MinIndex(~(MinValues<SearchRadius))=0;

                                if sum(MinIndex)
                                    Particles(PreviousFrameParticles).Frame(end+1)=CurrentFrame;
                                    Particles(PreviousFrameParticles).Index(end+1)=MinIndex;
                                    %We don't want this new spot to generate a
                                    %new particle further below
                                    NewParticleFlag(MinIndex)=false;
                                end
                            end
                        end


                        %See which spots weren't assigned and add them
                        %to the structure as new particles
                        NewParticles=find(NewParticleFlag);
                        for j=1:length(NewParticles)
                            TotalParticles=length(Particles{Channel});
                            Particles{Channel}(TotalParticles+1).Frame=CurrentFrame;
                            Particles{Channel}(TotalParticles+1).Index=...
                                         ApprovedSpots(NewParticles(j));
                            Particles{Channel}(TotalParticles+1).Approved=0;
                        end
                    end
                end

            %If we do have the histone channel    
            else
                if strcmpi(ExperimentType,'1spot')|strcmpi(ExperimentType,'inputoutput')|...
                        strcmpi(ExperimentType,'2spot2color')
                    SpotsPerNucleus=1;
                elseif strcmp(ExperimentType,'2spot')
                    SpotsPerNucleus=2;
                end

                [Particles{Channel},SpotFilter{Channel}]=...
                    AssignParticle2Nucleus(schnitzcells,Ellipses,Particles{Channel},Spots{Channel},SpotFilter{Channel},...
                    CurrentFrame,PixelSize,SearchRadius,SpotsPerNucleus);
            end
        end
    end
    close(ParticlesFig)
    close(NucleiFig)
%   
% %Double particle color mode
% elseif strcmp(ExperimentType,'2spot2color')
% 
%     %Load the FISH files. Split fad into particles that are above Threshold1 and
%     %those between that threshold and Threshold2.
% 
%     if ~exist('fad')
%         load([Folder,FileName])
%         fad=fishAnalysisData;
%         clear fishAnalysisData
% 
%         %We are assuming that all fad "frames" first correspond to all
%         %frames of channel one and the rest to all frames of channel 2.
%         
%         %Split the fad into the two channels
%         fadNew(1)=fad;
%         fadNew(2)=fad;
%         fadNew(1).channels=fadNew(1).channels(1:length(FrameInfo));
%         fadNew(2).channels=fadNew(2).channels(length(FrameInfo)+1:end);
%         fad=fadNew;
%         clear fadNew;
%         
%         fadTemp=fad;
%         clear fad
%         fad=fadTemp;
%         fad2=fadTemp;
%         
%         
%         for ChN=1:length(fadTemp)
%             for i=1:length(fadTemp(ChN).channels)
% 
% 
%                 Filter1=fadTemp(ChN).channels(i).fits.third>Threshold1(ChN);
%                 Filter2=(fadTemp(ChN).channels(i).fits.third<Threshold1(ChN))&...
%                     (fadTemp(ChN).channels(i).fits.third>Threshold2(ChN));
%                 Fields=fieldnames(fadTemp(ChN).channels(i).fits);
%                 NParticles=length(fadTemp(ChN).channels(i).fits.third);  %I'll use this to
%                                                                          %detemine the
%                                                                          %dimension of the
%                                                                          %different elements
%                                                                          %in the structure
% 
%                 for j=1:length(Fields)
%                     %In the end this seems to work, the fields that have higher
%                     %dimensions are still devided well.
%                     if strcmp(Fields{j},'snippets')
%                         Temp=getfield(fadTemp(ChN).channels(i).fits,Fields{j});
%                         SnippetDimentions=size(Temp);
%                         if length(SnippetDimentions)==3
%                             fad(ChN).channels(i).fits=setfield(fad(ChN).channels(i).fits,Fields{j},Temp(:,:,Filter1));
%                             fad2(ChN).channels(i).fits=setfield(fad2(ChN).channels(i).fits,Fields{j},Temp(:,:,Filter2));
%                         elseif length(SnippetDimentions)==2
%                             if Filter1
%                                 fad(ChN).channels(i).fits=setfield(fad(ChN).channels(i).fits,Fields{j},Temp(:,:,Filter1));
%                             elseif Filter2
%                                 fad2(ChN).channels(i).fits=setfield(fad2(ChN).channels(i).fits,Fields{j},Temp(:,:,Filter2));
%                             end
%                         end
%                     elseif strcmp(Fields{j},'maskUsedForTotalInt')
%                         Temp=getfield(fadTemp(ChN).channels(i).fits,Fields{j});
%                         fad(ChN).channels(i).fits=setfield(fad(ChN).channels(i).fits,Fields{j},Temp);
%                         fad2(ChN).channels(i).fits=setfield(fad2(ChN).channels(i).fits,Fields{j},Temp); 
%                     elseif length(getfield(fadTemp(ChN).channels(i).fits,Fields{j}))==NParticles
%                         Temp=getfield(fadTemp(ChN).channels(i).fits,Fields{j});
%                         fad(ChN).channels(i).fits=setfield(fad(ChN).channels(i).fits,Fields{j},Temp(Filter1));
%                         fad2(ChN).channels(i).fits=setfield(fad2(ChN).channels(i).fits,Fields{j},Temp(Filter2));
%                     else
%                         error('Missing field when separating fad and fad2')
%                     end
%                 end
%             end
%         end
%     end
% 
%     
%     %Check if particle tracking has already been done on this dataset
%     if exist([OutputFolder,filesep,'Particles.mat'])
% 
%         %error('Trying to re-track particles. This is not yet supported for 2spot2color mode.')
%         
%         load([OutputFolder,filesep,'Particles.mat'])
%         
%        
%         if isfield(Particles{1},'Approved')
%             Retracking=1;           %Flag for whether we are performing retracking
% 
%             %Only keep the approved particles and start the tracking from there
%             for ChN=1:length(fadTemp)
%                 k=1;
%                 clear NewParticles
%                 for i=1:length(Particles{ChN})
%                     if Particles{ChN}(i).Approved~=0
%                         NewParticles(k)=Particles{ChN}(i);
%                         k=k+1;
%                     end
%                 end
%                 if exist('NewParticles')
%                     Particles{ChN}=NewParticles;
%                 else
%                     Particles{ChN}=[];
%                     Retracking=0;
%                 end
%             end
%         else
%             Retracking=0;
%         end
%     else
%         Particles{length(fad)}=[];   %This is the structure where we'll be tracking all particles.
%         Retracking=0;
%     end
% 
%     
%     %Start by numbering the particles found
%     ParticlesFig=figure;
%     NucleiFig=figure;
% 
%     %Initially, only track particles that are above Threshold1
%     
%     for ChN=1:length(fad)
%         for i=1:length(fad(ChN).channels)
% 
% 
%             figure(ParticlesFig)
%             set(gcf,'Position',[ 16   369   676   342])
%             CurrentFrame=i;
%             [x,y]=fad2xyz(CurrentFrame,fad(ChN), 'addMargin'); 
%             CurrentZ=3;     
%             Image=imread([DataFolder,filesep,FilePrefix,iIndex(CurrentFrame,NDigits),'_z',iIndex(CurrentZ,2),...
%                 '_ch',iIndex(ChN,2),'.tif']);
%             imshow(Image,[])
%             hold on
%             plot(x,y,'or')
%             hold off
%             title(i)
%             set(gcf,'Name',['Frame: ',num2str(i),'/',num2str(length(fad(ChN).channels))])
% 
%             figure(NucleiFig)
%             set(gcf,'Position',[ 728   373   512   256])
%             CurrentFrame=i;
%             [x,y]=fad2xyz(CurrentFrame,fad(ChN), 'addMargin');
%             if UseHistone
%                 if exist([DataFolder,filesep,FilePrefix(1:end-1),'-His_',iIndex(CurrentFrame,NDigits),'.tif'])
%                     Image=imread([DataFolder,filesep,FilePrefix(1:end-1),'-His_',iIndex(CurrentFrame,NDigits),'.tif']);
%                 elseif exist([DataFolder,filesep,FilePrefix(1:end-1),'_',iIndex(CurrentFrame,NDigits),'_FakeHis.tif'])
%                     Image=imread([DataFolder,filesep,FilePrefix(1:end-1),'_',iIndex(CurrentFrame,NDigits),'_FakeHis.tif']);
%                 else
%                     Image=0;
%                 end
%                 imshow(Image,[],'Border','Tight')
%                 hold on
%                 PlotHandle=[];
%                 [NEllipses,Dummy]=size(Ellipses{CurrentFrame});
%                 for j=1:NEllipses
%                     PlotHandle=[PlotHandle,ellipse(Ellipses{CurrentFrame}(j,3),...
%                         Ellipses{CurrentFrame}(j,4),...
%                         Ellipses{CurrentFrame}(j,5),Ellipses{CurrentFrame}(j,1)+1,...
%                         Ellipses{CurrentFrame}(j,2)+1)];
%                         text(Ellipses{CurrentFrame}(j,1)+1,...
%                             Ellipses{CurrentFrame}(j,2)+1,num2str(j),'BackgroundColor',[.7 .9 .7]);
% 
%                 end
%                 set(PlotHandle,'Color','r')
%                 hold off
%                 title(i)
%                 drawnow
%             end
% 
% 
% 
% 
% 
% 
%             %If we don't have nuclear tracking then track the particles based on
%             %proximity
% 
%             if ~UseHistone
%                               
%                 error('This case needs to be adapted to multiple channels')
%                 
%                 %Get the particles detected in the frame
%                 if ~isempty(fad.channels(i).fits.x)
%                     if isempty(Particles)
%                         %Initialize the Particles structure if it doesn't exist yet
%                         for j=1:length(fad.channels(i).fits.x)
%                             Particles(j).Frame=i;
%                             Particles(j).Index=j;
%                             LastFrame(j)=i;
%                             Particles(j).Approved=0;
%                         end
% 
% 
%                     else
%                     %If we already have recorded particles, we need to compare
%                     %them to the new ones found and try to assign them. We also
%                     %need to make sure that they stay within the nuclei.
% 
%                         %Get the positions of the potentially new particles
%                         [NewParticlesX,NewParticlesY]=fad2xyz(i,fad, 'addMargin');
%                         NewParticlesZ=single(fad.channels(i).fits.z);
% 
%                         NewParticlesFlag=ones(size(NewParticlesX));
% 
% 
% 
%                         %Get a list of the particles that existed in the previous
%                         %frame
% 
%                         FilterPreviousFrame=(LastFrame==(i-1));
% 
%                         PreviousParticlesIndex=[];
%                         for j=1:length(FilterPreviousFrame)
%                             if FilterPreviousFrame(j)
%                                 PreviousParticlesIndex=...
%                                     [PreviousParticlesIndex,Particles(j).Index(end)];
%                             end
%                         end
% 
%                         PreviousParticles=find(FilterPreviousFrame);
% 
%                         if (~isempty(PreviousParticlesIndex))&(sum((PreviousParticlesIndex)))
%                             %Now, compare the positions of all of the old and new
%                             %particles
%                             clear Distance
%                             [PreviousParticlesX,PreviousParticlesY]=fad2xyz(i-1,fad, 'addMargin');
% 
%                             for j=1:length(NewParticlesX)
%                                 Distance(j,:)=sqrt((NewParticlesX(j)*PixelSize-...
%                                     PreviousParticlesX*PixelSize).^2+...
%                                     (NewParticlesY(j)*PixelSize-...
%                                     PreviousParticlesY*PixelSize).^2);
%                             end
%                             %The rows of distance correspond to the new particles.
%                             %The columns correspond to their distance to the old
%                             %particles.
%                             %MinIndex is a row vector. The element position
%                             %correspond to the new particle and the value within it
%                             %correspond to the old particle that it's closest to.
%                             [MinValues,MinIndex]=min(Distance');
% 
%                             %Now, use the information to figure out if we're dealing
%                             %with an old or a new particle. This will be based on who's
%                             %closest and if that distance is smaller than the threshold
%                             %given above in the parameters.
% 
%                             %These indices refer to the old particles. We're
%                             %finding which ones are good candidates                    
%                             UniqueMinima=unique(MinIndex);
% 
%                             for j=1:length(UniqueMinima)
%                                 %If we have only one previous and one new particle
%                                 %the assignment is trivial
%                                 if (length(PreviousParticles)==1)&(length(NewParticlesX)==1)&...
%                                         (MinValues<SearchRadius)
%                                     Particles(PreviousParticles).Frame(end+1)=i;
%                                     Particles(PreviousParticles).Index(end+1)=...
%                                             1;    
%                                     NewParticlesFlag(1)=0;
%                                     LastFrame(PreviousParticles)=i;
% 
%                                 %If we have only one previous particle MinIndex points at
%                                 %the new particle that is closest to it    
%                                 elseif (length(PreviousParticles)==1)&(MinValues<SearchRadius)
%                                     Particles(PreviousParticles).Frame(end+1)=i;
%                                     Particles(PreviousParticles).Index(end+1)=MinIndex;
%                                     NewParticlesFlag(MinIndex)=0;
%                                     LastFrame(PreviousParticles)=i;
% 
%                                 %If there is only one new particle that a previous
%                                 %particle is closest to
%                                 elseif (sum(MinIndex==UniqueMinima(j))==1)
%                                     if MinValues(find(MinIndex==UniqueMinima(j)))<(SearchRadius)
%                                         Particles(PreviousParticles(PreviousParticlesIndex==UniqueMinima(j))).Frame(end+1)=i;
%                                         Particles(PreviousParticles(PreviousParticlesIndex==UniqueMinima(j))).Index(end+1)=...
%                                             find(MinIndex==UniqueMinima(j));
%                                         NewParticlesFlag(find(MinIndex==UniqueMinima(j)))=0;
%                                         LastFrame(PreviousParticles(PreviousParticlesIndex==UniqueMinima(j)))=i;
%                                     end
% 
%                                 %If there are many new particles that are close to a previous one.    
%                                 else
%                                     MinFilter=find(MinIndex==UniqueMinima(j));
%                                     MinValues2=MinValues(MinFilter);
%                                     [MinMinValue2,MinIndex2]=min(MinValues2);
% 
%                                     if MinMinValue2<(SearchRadius)
%                                         ClosestParticle=MinFilter(MinIndex2);
% 
%                                         %This is just in case there are two
%                                         %previous particles at the exact same
%                                         %distance
%                                         if length(ClosestParticle)>1
%                                             ClosestParticle=ClosestParticle(1);
%                                         end
% 
%                                         Particles(PreviousParticles(PreviousParticlesIndex==UniqueMinima(j))).Frame(end+1)=i;
%                                         Particles(PreviousParticles(PreviousParticlesIndex==UniqueMinima(j))).Index(end+1)=...
%                                             ClosestParticle;
%                                         NewParticlesFlag(ClosestParticle)=0;
% 
%                                         LastFrame(PreviousParticles(PreviousParticlesIndex==UniqueMinima(j)))=i;
%                                     end
%                                 end
%                             end
%                         end
% 
%                         %Finally, see which particles weren't assigned and add them
%                         %to the structure
%                         NewParticles=find(NewParticlesFlag);
%                         for j=1:length(NewParticles)
%                             TotalParticles=length(Particles);
%                             Particles(TotalParticles+1).Frame=i;
%                             Particles(TotalParticles+1).Index=...
%                                          NewParticles(j);
%                             LastFrame(TotalParticles+1)=i;
%                             Particles(TotalParticles+1).Approved=0;
%                         end
%                     end
%                 end
%             else
%                 %If we do have the histone channel    
%                 [Particles{ChN},fad(ChN),fad2(ChN),schnitzcells]=AssignParticle2Nucleus(schnitzcells,Ellipses,Particles{ChN},fad(ChN),fad2(ChN),...
%                             i,PixelSize,SearchRadius);
% 
%             end
%         end
%     end
% 
% 
%     close(ParticlesFig)
%     close(NucleiFig)
% 
% 
%     if ~UseHistone
%         
%         error('This mode of the code has not been adapted to two-particles, two-colors')
%         
%         %Go through the particles we found and try to fill in the gaps with the
%         %Threshold2 ones in fad2.
% 
%         for i=1:length(Particles)
% 
%             %Move forward in time
%             CurrentParticleLength=length(Particles(i).Frame)-1;
% 
%             while (CurrentParticleLength<length(Particles(i).Frame))&...
%                     (Particles(i).Frame(end)<length(fad.channels))
%                 CurrentParticleLength=length(Particles(i).Frame);
%                 CurrentFrame=Particles(i).Frame(end);
% 
%                 [fad,fad2,Particles] = ConnectToThreshold(fad,fad2,Particles,...
%                     i,CurrentFrame,CurrentFrame+1,SearchRadius*0.5);
%             end
% 
%             %Move backward in time
%             CurrentParticleLength=length(Particles(i).Frame)-1;
% 
%             while (CurrentParticleLength<length(Particles(i).Frame))&(Particles(i).Frame(1)>1)
%                 CurrentParticleLength=length(Particles(i).Frame);
%                 CurrentFrame=Particles(i).Frame(1);
% 
%                 %[2,CurrentFrame]
% 
%                 [fad,fad2,Particles] = ConnectToThreshold(fad,fad2,Particles,...
%                     i,CurrentFrame,CurrentFrame-1,SearchRadius*0.5);
%             end
% 
%         end
% 
%         mkdir([OutputFolder,filesep]);
% 
%         save([OutputFolder,filesep,'Particles.mat'],'Particles','SpotFilter',...
%             'Threshold1','Threshold2');
% 
% 
%     else
%         
%         %If we have the histone channel we now have the extra information of
%         %the nuclei. In this case we'll admit one particle per nucleus
%         %and we'll connect it to the closest particles within that nucleus
% 
%         %As a result we will scan through the nuclei trying to fill in any gaps
%         for ChN=1:length(fad)
%             for i=1:length(Particles{ChN})
%                 clear ParticleNuclei
%                 if ~isempty(Particles{ChN}(i).Nucleus)
%                     ParticleNuclei(i)=Particles{ChN}(i).Nucleus;
%                 else
%                     ParticleNuclei(i)=0;
%                 end
%             end
% 
% 
% 
%             h=waitbar(0,'Checking secondary threshold');
%             for i=1:length(schnitzcells)
% 
% 
%                 waitbar(i/length(schnitzcells),h);
% 
%                 [i,length(schnitzcells)];
% 
%                 %See how many particles were detected within this nucleus
%                 ParticlesInNucleus=find(ParticleNuclei==i);
% 
% 
%                 for k=1:length(ParticlesInNucleus)
%                     %If there is only one particle in that nucleus then this is
%                     %relatively easy
%                     CurrentParticle=ParticlesInNucleus(k);
% 
%                     %Flag to see whether we'll check this particle
%                     CheckParticle=0;
%                     if ~Retracking
%                         CheckParticle=1;
%                     elseif ~Particles{ChN}(CurrentParticle).Approved
%                         CheckParticle=1;
%                     end
% 
% 
% 
%                     if CheckParticle
%                         %Move forward in time
%                         CurrentParticleLength=length(Particles{ChN}(CurrentParticle).Frame)-1;
%                         while (CurrentParticleLength<length(Particles{ChN}(CurrentParticle).Frame))&...
%                                 (Particles{ChN}(CurrentParticle).Frame(end)<length(Ellipses))
%                             CurrentParticleLength=length(Particles{ChN}(CurrentParticle).Frame);
%                             CurrentFrame=Particles{ChN}(CurrentParticle).Frame(end);
% 
%                             [fad(ChN),fad2(ChN),Particles{ChN}] = ConnectToThresholdHistone(fad(ChN),fad2(ChN),Particles{ChN},...
%                                 CurrentParticle,CurrentFrame,CurrentFrame+1,schnitzcells,Ellipses,SearchRadius,...
%                                 PixelSize);
%                         end
% 
% 
%                         %Move backwards in time
%                         CurrentParticleLength=length(Particles{ChN}(CurrentParticle).Frame)-1;
%                         while (CurrentParticleLength<length(Particles{ChN}(CurrentParticle).Frame))&...
%                                 (Particles{ChN}(CurrentParticle).Frame(1)>1)
%                             CurrentParticleLength=length(Particles{ChN}(CurrentParticle).Frame);
%                             CurrentFrame=Particles{ChN}(CurrentParticle).Frame(1);
% 
%                             [fad(ChN),fad2(ChN),Particles{ChN}] = ConnectToThresholdHistone(fad(ChN),fad2(ChN),Particles{ChN},...
%                                 CurrentParticle,CurrentFrame,CurrentFrame-1,schnitzcells,Ellipses,SearchRadius,...
%                                 PixelSize);
%                         end
% 
% 
%                         %Fill in any gaps if there are any
% 
%                         %Find which ones are the missing frames
%                         MissingFrames=Particles{ChN}(CurrentParticle).Frame(1):Particles{ChN}(CurrentParticle).Frame(end);
%                         MissingFrames=MissingFrames(~ismember(MissingFrames,Particles{ChN}(CurrentParticle).Frame));
% 
%                         if ~isempty(MissingFrames)            
%                             for j=1:length(MissingFrames)
%                                 [fad(ChN),fad2(ChN),Particles{ChN}] = ConnectToThresholdHistone(fad(ChN),fad2(ChN),Particles{ChN},...
%                                     CurrentParticle,[],MissingFrames(j),schnitzcells,Ellipses,SearchRadius,...
%                                     PixelSize);
%                             end
%                         end
%                     end
%                 end
% 
% 
%                 %After the for loop, try to join particles in the same nucleus, but
%                 %different frames. Do this only if not retracking.
%                 if ~Retracking
%                     if length(ParticlesInNucleus)>1
% 
%                         PreviousParticleNucleiLength=length(ParticleNuclei)+1;
% 
%                         while PreviousParticleNucleiLength>length(ParticleNuclei)
%                             PreviousParticleNucleiLength=length(ParticleNuclei);
%                             Particles{ChN}=FillParticleGaps(ParticlesInNucleus,Particles{ChN});
%                             ParticleNuclei=[Particles{ChN}.Nucleus];
%                             ParticlesInNucleus=find(ParticleNuclei==i);
%                         end
%                     end
%                 end
% 
% 
% 
% 
%             end
% 
% 
%             close(h)
%         end
%             
%         mkdir([OutputFolder,filesep])
% 
%         save([OutputFolder,filesep,'Particles.mat'],'Particles','SpotFilter',...
%             'Threshold1','Threshold2')
% 
%     end
%     
%     
% % elseif strcmp(lower(ExperimentType),'inputoutput')
% % 
% %     %Check if particle tracking has already been done on this dataset
% %     if exist([OutputFolder,filesep,'Particles.mat'])
% % 
% %         load([OutputFolder,filesep,'Particles.mat'])
% %         Particles=Particles{1};
% %         
% %         if isfield(Particles,'Approved')
% %             Retracking=1;           %Flag for whether we are performing retracking
% % 
% %             %Only keep the approved particles and start the tracking from there
% %             k=1;
% %             for i=1:length(Particles)
% %                 if Particles(i).Approved~=0
% %                     NewParticles(k)=Particles(i);
% %                     k=k+1;
% %                 end
% %             end
% %             if exist('NewParticles')
% %                 Particles=NewParticles;
% %             else
% %                 Particles=[];
% %                 Retracking=0;
% %             end
% %         else
% %             Retracking=0;
% %             Particles=[];
% %         end
% %     else
% %         Particles=[];   %This is the structure where we'll be tracking all particles.
% %         Retracking=0;
% %     end
% % 
% % 
% %     %Load the FISH files. Split fad into particles that are above Threshold1 and
% %     %those between that threshold and Threshold2.
% % 
% %     if ~exist('fad')
% %         load([Folder,FileName])
% %         fad=fishAnalysisData;
% %         clear fishAnalysisData
% % 
% %         fadTemp=fad;
% %         clear fad
% %         fad=fadTemp;
% %         fad2=fadTemp;
% %         for i=1:length(fadTemp.channels)
% % 
% % 
% %             Filter1=fadTemp.channels(i).fits.third>Threshold1;
% %             Filter2=(fadTemp.channels(i).fits.third<Threshold1)&...
% %                 (fadTemp.channels(i).fits.third>Threshold2);
% %             Fields=fieldnames(fadTemp.channels(i).fits);
% %             NParticles=length(fadTemp.channels(i).fits.third);  %I'll use this to
% %                                                                 %detemine the
% %                                                                 %dimension of the
% %                                                                 %different elements
% %                                                                 
% %             for j=1:length(Fields)
% %                 %In the end this seems to work, the fields that have higher
% %                 %dimensions are still devided well.
% %                 if strcmp(Fields{j},'snippets')
% %                     Temp=getfield(fadTemp.channels(i).fits,Fields{j});
% %                     SnippetDimentions=size(Temp);
% %                     if length(SnippetDimentions)==3
% %                         fad.channels(i).fits=setfield(fad.channels(i).fits,Fields{j},Temp(:,:,Filter1));
% %                         fad2.channels(i).fits=setfield(fad2.channels(i).fits,Fields{j},Temp(:,:,Filter2));
% %                     elseif length(SnippetDimentions)==2
% %                         if Filter1
% %                             fad.channels(i).fits=setfield(fad.channels(i).fits,Fields{j},Temp(:,:,Filter1));
% %                         elseif Filter2
% %                             fad2.channels(i).fits=setfield(fad2.channels(i).fits,Fields{j},Temp(:,:,Filter2));
% %                         end
% %                     end
% %                 elseif strcmp(Fields{j},'maskUsedForTotalInt')
% %                     Temp=getfield(fadTemp.channels(i).fits,Fields{j});
% %                     fad.channels(i).fits=setfield(fad.channels(i).fits,Fields{j},Temp);
% %                     fad2.channels(i).fits=setfield(fad2.channels(i).fits,Fields{j},Temp); 
% %                 elseif length(getfield(fadTemp.channels(i).fits,Fields{j}))==NParticles
% %                     Temp=getfield(fadTemp.channels(i).fits,Fields{j});
% %                     fad.channels(i).fits=setfield(fad.channels(i).fits,Fields{j},Temp(Filter1));
% %                     fad2.channels(i).fits=setfield(fad2.channels(i).fits,Fields{j},Temp(Filter2));
% %                 else
% %                     error('Missing field when separating fad and fad2')
% %                 end
% %             end
% %         end
% %     end
% % 
% %     %Start by numbering the particles found
% %     ParticlesFig=figure;
% %     NucleiFig=figure;
% % 
% %     
% %     %Figure out which channel has the spots
% %     if ~isempty(strfind(lower(Channel1{1}),'mcp'))|~isempty(strfind(lower(Channel1{1}),'pcp'))
% %         SpotChannel=1;
% %     elseif ~isempty(strfind(lower(Channel2{1}),'mcp'))|~isempty(strfind(lower(Channel2{1}),'pcp'))
% %         SpotChannel=2;
% %     else
% %         error('Could not identify spot channel in MovieDatabase.xlsx. Are they named using "MCP" or "PCP"?')
% %     end
% %     
% %     
% %     %Initially, only track particles that are above Threshold1
% %     for i=1:length(fad.channels)
% % 
% % 
% %         figure(ParticlesFig)
% %         set(gcf,'Position',[ 16   369   676   342])
% %         CurrentFrame=i;
% %         [x,y]=fad2xyz(CurrentFrame,fad, 'addMargin'); 
% %         CurrentZ=3;     
% %         Image=imread([PreProcPath,filesep,Prefix,filesep,FilePrefix,iIndex(CurrentFrame,NDigits),'_z',iIndex(CurrentZ,2),'_ch',iIndex(SpotChannel,2),'.tif']);
% %         
% %         
% %         imshow(Image,[])
% %         hold on
% %         plot(x,y,'or')
% %         hold off
% %         title(i)
% %         set(gcf,'Name',['Frame: ',num2str(i),'/',num2str(length(fad.channels))])
% % 
% %         figure(NucleiFig)
% %         set(gcf,'Position',[ 728   373   512   256])
% %         CurrentFrame=i;
% %         [x,y]=fad2xyz(CurrentFrame,fad, 'addMargin');
% %         if UseHistone
% %             try
% %                 Image=imread([PreProcPath,filesep,Prefix,filesep,FilePrefix(1:end-1),'-His_',iIndex(CurrentFrame,NDigits),'.tif']);
% %             catch
% %                 try
% %                     Image=imread([PreProcPath,filesep,Prefix,filesep,FilePrefix(1:end-1),'_His_',iIndex(CurrentFrame,NDigits),'.tif']);
% %                 catch
% %                     Image=0;
% %                 end
% %             end
% %             imshow(Image,[],'Border','Tight')
% %             hold on
% %             PlotHandle=[];
% %             [NEllipses,Dummy]=size(Ellipses{CurrentFrame});
% %             for j=1:NEllipses
% %                 PlotHandle=[PlotHandle,ellipse(Ellipses{CurrentFrame}(j,3),...
% %                     Ellipses{CurrentFrame}(j,4),...
% %                     Ellipses{CurrentFrame}(j,5),Ellipses{CurrentFrame}(j,1)+1,...
% %                     Ellipses{CurrentFrame}(j,2)+1)];
% %                     text(Ellipses{CurrentFrame}(j,1)+1,...
% %                         Ellipses{CurrentFrame}(j,2)+1,num2str(j),'BackgroundColor',[.7 .9 .7]);
% % 
% %             end
% %             set(PlotHandle,'Color','r')
% %             hold off
% %             title(i)
% %         end
% %         drawnow
% %         
% %         %If we don't have nuclear tracking then track the particles based on
% %         %proximity
% % 
% %         if ~UseHistone
% %             %Get the particles detected in the frame
% %             if ~isempty(fad.channels(i).fits.x)
% %                 if isempty(Particles)
% %                     %Initialize the Particles structure if it doesn't exist yet
% %                     for j=1:length(fad.channels(i).fits.x)
% %                         Particles(j).Frame=i;
% %                         Particles(j).Index=j;
% %                         LastFrame(j)=i;
% %                         Particles(j).Approved=0;
% %                     end
% % 
% % 
% %                 else
% %                 %If we already have recorded particles, we need to compare
% %                 %them to the new ones found and try to assign them. We also
% %                 %need to make sure that they stay within the nuclei.
% % 
% %                     %Get the positions of the potentially new particles
% %                     [NewParticlesX,NewParticlesY]=fad2xyz(i,fad, 'addMargin');
% %                     NewParticlesZ=single(fad.channels(i).fits.z);
% % 
% %                     NewParticlesFlag=ones(size(NewParticlesX));
% % 
% % 
% % 
% %                     %Get a list of the particles that existed in the previous
% %                     %frame
% % 
% %                     FilterPreviousFrame=(LastFrame==(i-1));
% % 
% %                     PreviousParticlesIndex=[];
% %                     for j=1:length(FilterPreviousFrame)
% %                         if FilterPreviousFrame(j)
% %                             PreviousParticlesIndex=...
% %                                 [PreviousParticlesIndex,Particles(j).Index(end)];
% %                         end
% %                     end
% % 
% %                     PreviousParticles=find(FilterPreviousFrame);
% % 
% %                     if (~isempty(PreviousParticlesIndex))&(sum((PreviousParticlesIndex)))
% %                         %Now, compare the positions of all of the old and new
% %                         %particles
% %                         clear Distance
% %                         [PreviousParticlesX,PreviousParticlesY]=fad2xyz(i-1,fad, 'addMargin');
% % 
% %                         for j=1:length(NewParticlesX)
% %                             Distance(j,:)=sqrt((NewParticlesX(j)*PixelSize-...
% %                                 PreviousParticlesX*PixelSize).^2+...
% %                                 (NewParticlesY(j)*PixelSize-...
% %                                 PreviousParticlesY*PixelSize).^2);
% %                         end
% %                         %The rows of distance correspond to the new particles.
% %                         %The columns correspond to their distance to the old
% %                         %particles.
% %                         %MinIndex is a row vector. The element position
% %                         %correspond to the new particle and the value within it
% %                         %correspond to the old particle that it's closest to.
% %                         [MinValues,MinIndex]=min(Distance');
% % 
% %                         %Now, use the information to figure out if we're dealing
% %                         %with an old or a new particle. This will be based on who's
% %                         %closest and if that distance is smaller than the threshold
% %                         %given above in the parameters.
% % 
% %                         %These indices refer to the old particles. We're
% %                         %finding which ones are good candidates                    
% %                         UniqueMinima=unique(MinIndex);
% % 
% %                         if i==31
% %                             1+1;
% %                         end
% % 
% % 
% % 
% %                         for j=1:length(UniqueMinima)
% %                             %If we have only one previous and one new particle
% %                             %the assignment is trivial
% %                             if (length(PreviousParticles)==1)&(length(NewParticlesX)==1)&...
% %                                     (MinValues<SearchRadius)
% %                                 Particles(PreviousParticles).Frame(end+1)=i;
% %                                 Particles(PreviousParticles).Index(end+1)=...
% %                                         1;    
% %                                 NewParticlesFlag(1)=0;
% %                                 LastFrame(PreviousParticles)=i;
% % 
% %                             %If we have only one previous particle MinIndex points at
% %                             %the new particle that is closest to it    
% %                             elseif (length(PreviousParticles)==1)&(MinValues<SearchRadius)
% %                                 Particles(PreviousParticles).Frame(end+1)=i;
% %                                 Particles(PreviousParticles).Index(end+1)=MinIndex;
% %                                 NewParticlesFlag(MinIndex)=0;
% %                                 LastFrame(PreviousParticles)=i;
% % 
% %                             %If there is only one new particle that a previous
% %                             %particle is closest to
% %                             elseif (sum(MinIndex==UniqueMinima(j))==1)
% %                                 if MinValues(find(MinIndex==UniqueMinima(j)))<(SearchRadius)
% %                                     Particles(PreviousParticles(PreviousParticlesIndex==UniqueMinima(j))).Frame(end+1)=i;
% %                                     Particles(PreviousParticles(PreviousParticlesIndex==UniqueMinima(j))).Index(end+1)=...
% %                                         find(MinIndex==UniqueMinima(j));
% %                                     NewParticlesFlag(find(MinIndex==UniqueMinima(j)))=0;
% %                                     LastFrame(PreviousParticles(PreviousParticlesIndex==UniqueMinima(j)))=i;
% %                                 end
% % 
% %                             %If there are many new particles that are close to a previous one.    
% %                             else
% %                                 MinFilter=find(MinIndex==UniqueMinima(j));
% %                                 MinValues2=MinValues(MinFilter);
% %                                 [MinMinValue2,MinIndex2]=min(MinValues2);
% % 
% %                                 if MinMinValue2<(SearchRadius)
% %                                     ClosestParticle=MinFilter(MinIndex2);
% % 
% %                                     %This is just in case there are two
% %                                     %previous particles at the exact same
% %                                     %distance
% %                                     if length(ClosestParticle)>1
% %                                         ClosestParticle=ClosestParticle(1);
% %                                     end
% % 
% %                                     Particles(PreviousParticles(PreviousParticlesIndex==UniqueMinima(j))).Frame(end+1)=i;
% %                                     Particles(PreviousParticles(PreviousParticlesIndex==UniqueMinima(j))).Index(end+1)=...
% %                                         ClosestParticle;
% %                                     NewParticlesFlag(ClosestParticle)=0;
% % 
% %                                     LastFrame(PreviousParticles(PreviousParticlesIndex==UniqueMinima(j)))=i;
% %                                 end
% %                             end
% %                         end
% %                     end
% % 
% %                     %Finally, see which particles weren't assigned and add them
% %                     %to the structure
% %                     NewParticles=find(NewParticlesFlag);
% %                     for j=1:length(NewParticles)
% %                         TotalParticles=length(Particles);
% %                         Particles(TotalParticles+1).Frame=i;
% %                         Particles(TotalParticles+1).Index=...
% %                                      NewParticles(j);
% %                         LastFrame(TotalParticles+1)=i;
% %                         Particles(TotalParticles+1).Approved=0;
% %                     end
% %                 end
% %             end
% % 
% %         %If we do have the histone channel    
% %         else
% %             [Particles,fad,fad2,schnitzcells]=AssignParticle2Nucleus(schnitzcells,Ellipses,Particles,fad,fad2,...
% %                 i,PixelSize,SearchRadius);
% %         end
% %     end
% % 
% %     close(ParticlesFig)
% %     close(NucleiFig)
% % 
% % 
% %     if ~UseHistone
% %         %Go through the particles we found and try to fill in the gaps with the
% %         %Threshold2 ones in fad2.
% % 
% %         for i=1:length(Particles)
% % 
% %             if i==199
% %                 1+1
% %             end
% % 
% %             %Move forward in time
% %             CurrentParticleLength=length(Particles(i).Frame)-1;
% % 
% %             while (CurrentParticleLength<length(Particles(i).Frame))&...
% %                     (Particles(i).Frame(end)<length(fad.channels))
% %                 CurrentParticleLength=length(Particles(i).Frame);
% %                 CurrentFrame=Particles(i).Frame(end);
% % 
% %                 [fad,fad2,Particles] = ConnectToThreshold(fad,fad2,Particles,...
% %                     i,CurrentFrame,CurrentFrame+1,SearchRadius*0.5);
% %             end
% % 
% %             %Move backward in time
% %             CurrentParticleLength=length(Particles(i).Frame)-1;
% % 
% %             while (CurrentParticleLength<length(Particles(i).Frame))&(Particles(i).Frame(1)>1)
% %                 CurrentParticleLength=length(Particles(i).Frame);
% %                 CurrentFrame=Particles(i).Frame(1);
% % 
% %                 %[2,CurrentFrame]
% % 
% %                 [fad,fad2,Particles] = ConnectToThreshold(fad,fad2,Particles,...
% %                     i,CurrentFrame,CurrentFrame-1,SearchRadius*0.5);
% %             end
% % 
% %         end
% % 
% %         mkdir([OutputFolder,filesep]);
% % 
% %         save([OutputFolder,filesep,'Particles.mat'],'Particles','SpotFilter',...
% %             'Threshold1','Threshold2');
% % 
% % 
% %     else
% % 
% %         if strcmp(ExperimentType,'1spot')
% % 
% % 
% %             %if ~Retracking
% %                 %If we have the histone channel we now have the extra information of
% %                 %the nuclei. In this case we'll admit one particle per nucleus
% %                 %and we'll connect it to the closest particles within that nucleus
% % 
% %                 %As a result we will scan through the nuclei trying to fill in any gaps
% %                 for i=1:length(Particles)
% %                     if ~isempty(Particles(i).Nucleus)
% %                         ParticleNuclei(i)=Particles(i).Nucleus;
% %                     else
% %                         ParticleNuclei(i)=0;
% %                     end
% %                 end
% % 
% % 
% % 
% %                 h=waitbar(0,'Checking secondary threshold');
% %                 for i=1:length(schnitzcells)
% % 
% % 
% %                     waitbar(i/length(schnitzcells),h);
% % 
% %                     [i,length(schnitzcells)];
% % 
% %                     %See how many particles were detected within this nucleus
% %                     ParticlesInNucleus=find(ParticleNuclei==i);
% % 
% % 
% %                     for k=1:length(ParticlesInNucleus)
% %                         %If there is only one particle in that nucleus then this is
% %                         %relatively easy
% %                         CurrentParticle=ParticlesInNucleus(k);
% % 
% %                         %Flag to see whether we'll check this particle
% %                         CheckParticle=0;
% %                         if ~Retracking
% %                             CheckParticle=1;
% %                         elseif ~Particles(CurrentParticle).Approved
% %                             CheckParticle=1;
% %                         end
% % 
% %                         if CheckParticle
% %                             %Move forward in time
% %                             CurrentParticleLength=length(Particles(CurrentParticle).Frame)-1;
% %                             while (CurrentParticleLength<length(Particles(CurrentParticle).Frame))&...
% %                                     (Particles(CurrentParticle).Frame(end)<length(Ellipses))
% %                                 CurrentParticleLength=length(Particles(CurrentParticle).Frame);
% %                                 CurrentFrame=Particles(CurrentParticle).Frame(end);
% % 
% %                                 [fad,fad2,Particles] = ConnectToThresholdHistone(fad,fad2,Particles,...
% %                                     CurrentParticle,CurrentFrame,CurrentFrame+1,schnitzcells,Ellipses,SearchRadius,...
% %                                     PixelSize);
% %                             end
% % 
% % 
% %                             %Move backwards in time
% %                             CurrentParticleLength=length(Particles(CurrentParticle).Frame)-1;
% %                             while (CurrentParticleLength<length(Particles(CurrentParticle).Frame))&...
% %                                     (Particles(CurrentParticle).Frame(1)>1)
% %                                 CurrentParticleLength=length(Particles(CurrentParticle).Frame);
% %                                 CurrentFrame=Particles(CurrentParticle).Frame(1);
% % 
% %                                 [fad,fad2,Particles] = ConnectToThresholdHistone(fad,fad2,Particles,...
% %                                     CurrentParticle,CurrentFrame,CurrentFrame-1,schnitzcells,Ellipses,SearchRadius,...
% %                                     PixelSize);
% %                             end
% % 
% % 
% %                             %Fill in any gaps if there are any
% % 
% %                             %Find which ones are the missing frames
% %                             MissingFrames=Particles(CurrentParticle).Frame(1):Particles(CurrentParticle).Frame(end);
% %                             MissingFrames=MissingFrames(~ismember(MissingFrames,Particles(CurrentParticle).Frame));
% % 
% %                             if ~isempty(MissingFrames)            
% %                                 for j=1:length(MissingFrames)
% %                                     [fad,fad2,Particles] = ConnectToThresholdHistone(fad,fad2,Particles,...
% %                                         CurrentParticle,[],MissingFrames(j),schnitzcells,Ellipses,SearchRadius,...
% %                                         PixelSize);
% %                                 end
% %                             end
% %                         end
% %                     end
% % 
% % 
% %                     %After the for loop, try to join particles in the same nucleus, but
% %                     %different frames. Do this only if not retracking.
% %                     if ~Retracking
% %                         if length(ParticlesInNucleus)>1
% % 
% %                             PreviousParticleNucleiLength=length(ParticleNuclei)+1;
% % 
% %                             while PreviousParticleNucleiLength>length(ParticleNuclei)
% %                                 PreviousParticleNucleiLength=length(ParticleNuclei);
% %                                 Particles=FillParticleGaps(ParticlesInNucleus,Particles);
% %                                 ParticleNuclei=[Particles.Nucleus];
% %                                 ParticlesInNucleus=find(ParticleNuclei==i);
% %                             end
% %                         end
% %                     end
% % 
% % 
% % 
% % 
% %                 end
% % 
% % 
% %                 close(h)
% %             %end    
% % 
% % 
% %         end
% %  
% %         mkdir([OutputFolder,filesep])
% % 
% %         save([OutputFolder,filesep,'Particles.mat'],'Particles','SpotFilter',...
% %             'Threshold1','Threshold2')
% % 
% %     end    
else
    error('Experiment type in MovieDatabase.xlsx not recognized')    
end

%If we only have one channel, then convert SpotFilter and Particles to a
%standard structure.
if NCh==1
   SpotFilter=SpotFilter{1};
   Particles=Particles{1};
end

mkdir([OutputFolder,filesep]);

save([OutputFolder,filesep,'Particles.mat'],'Particles','SpotFilter',...
    'Threshold1','Threshold2');