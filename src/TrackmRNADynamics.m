function [Particles,schnitzcells]=TrackmRNADynamics(varargin)
%[Particles,schnitzcells]=TrackmRNADynamics(varargin)
%
% DESCRIPTION
% %This function tracks transcription loci over time after 
% segmentation and z-tracking have been performed. If nuclei have been 
% tracked it uses that information for the particle tracking.
%
% ARGUMENTS
% Prefix: Prefix of the data set to analyze
% Threshold1 : Primary. This must be an array of two thresholds for 2 spot 2
% color experiments. 
% Threshold2 : Secondary. This must be an array of two thresholds for 2 spot 2
% color experiments. 
% [Options]: See below.
%
% OPTIONS
% None.
%
% OUTPUT
% Particles.mat : List of time traces found in the movie. 
% *_lin.mat : Schnitzcell nuclear tracking is modified here. 
%
% Author (contact): Hernan Garcia (hggarcia@berkeley.edu)
% Created: 01/01/2013 ish. 
% Last Updated: 9/11/2018
%
% Documented by: Armando Reimer (areimer@berkeley.edu)
%
%
%
%To do: The no-histone part of the code doesn't take into account the
%Approved field of the Particles structure. 
%^ AR 9/3/18: has this been done? 
%

[SourcePath,FISHPath,DefaultDropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders;

app = {};
bypassUserPrompt = false;

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
    
    %2 spot channel input should be like TrackmRNADynamics(Prefix,
    %[thresh1, thresh2], [thresh3, thresh4])
    Threshold1=varargin{2};
    Threshold2=varargin{3}; 
    
    for i = 4:length(varargin)
        if strcmpi(varargin{i}, 'app')
          app{1} = varargin{i+1};
          app{2} = varargin{i+2};
        elseif strcmpi(varargin{i}, 'bypassUserPrompt')
          bypassUserPrompt = true;
        end
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

%What type of experiment are we dealing with? Get this out of MovieDatabase
[Date, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution,...
Channel1, Channel2, Objective, Power, DataFolder, DropboxFolderName, Comments,...
nc9, nc10, nc11, nc12, nc13, nc14, CF] = getExperimentDataFromMovieDatabase(Prefix, DefaultDropboxFolder)

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
            if ~bypassUserPrompt
              Answer=input('Thresholds changed, will delete previous tracking. Proceed? (y/n):','s');
            else 
              Answer = 'y';
            end
            
            if strcmpi(Answer,'y')
                Threshold1=Threshold1Backup;
                Threshold2=Threshold2Backup;
                %save a backup of the older version just in case. 
                mkdir([OutputFolder,filesep,'OldParticlesVersions']);
                save([OutputFolder,filesep,'OldParticlesVersions',filesep,'Particles.mat']);
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
    
        if (FrameInfo(1).ZoomFactor==8)&(FrameInfo(1).PixelsPerLine==256) %#ok<*AND2>
            PixelSize=0.22;     %This is the pixel size in um for a zoom of 8.
        elseif (FrameInfo(1).ZoomFactor==4)&...
                (FrameInfo(1).PixelsPerLine==512)&...
                (FrameInfo(1).ScanAmplitudeY==1)
            PixelSize=0.22;
        elseif (FrameInfo(1).ZoomFactor==16)&(FrameInfo(1).PixelsPerLine==128)
            PixelSize=0.22;
        else
            disp('Warning: Imaging setting not defined. Using a pixel size of 0.22um')
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
            disp('Warning: Imaging setting not defined. Using a pixel size of 0.22um')
            PixelSize=0.22;
        end
        
    elseif strcmp(FrameInfo(1).FileMode,'LSM')|strcmp(FrameInfo(1).FileMode,'LSMExport') %#ok<*OR2>
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
if exist([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'], 'file')
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
%(MT, 2018-02-11) Added support for lattice imaging, maybe temporary -
%FIX LATER
if strcmpi(ExperimentType,'1spot')||strcmpi(ExperimentType,'2spot')||...
        strcmpi(ExperimentType,'inputoutput')||...
        strcmpi(ExperimentType,'2spot2color')||...
        strcmpi(ExperimentType,'lattice')

    %Figure out which channel has the spots
    if strcmp(ExperimentType,'1spot')||strcmp(ExperimentType,'2spot')
        SpotsChannel=1;
    %(MT, 2018-02-11) Added support for lattice imaging, maybe temporary -
    %FIX LATER
    elseif strcmp(ExperimentType,'inputoutput') || strcmp(ExperimentType,'lattice')
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
    if exist([OutputFolder,filesep,'Particles.mat'], 'file')

        load([OutputFolder,filesep,'Particles.mat'])
        
        %If there's only one channel, Particles, Spots and other structures are
        %not saved as cells. We turn them into a cell to ensure
        %compatibility.
        if NCh==1
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
    if ~exist('Spots', 'var')
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
            if exist('SpotFilter','var')
                if ~iscell(SpotFilter)
                    SpotFilter = {SpotFilter};
                end
            end
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
    if isempty(app)
        ParticlesFig=figure;
        particlesAxes = axes(ParticlesFig);
        if UseHistone
            NucleiFig=figure;
            set(NucleiFig,'units', 'normalized', 'position',[0.65, .5, .2, .2])
            nucAxes = axes(NucleiFig);
        end
    end

    
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

            if isempty(app)
                figure(ParticlesFig)
                set(gcf,'units', 'normalized', 'position',[0.01, .55, .33, .33]);
            end

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
            if ~isempty(app)
                ax1 = app{1};
                title(ax1, ['Ch', num2str(Channel),'  Frame: ',num2str(CurrentFrame),'/',num2str(length(Spots{Channel}))])
            else
                ax1 = particlesAxes;
                set(ParticlesFig,'Name',['Ch', num2str(Channel),'  Frame: ',num2str(CurrentFrame),'/',num2str(length(Spots{Channel}))]); 
                title(ax1, i)
            end
            imshow(Image,[], 'Parent', ax1, 'InitialMagnification', 'fit')
            hold(ax1, 'on')
            plot(ax1,x(CurrentFrameFilter),y(CurrentFrameFilter),'or', 'MarkerSize', 10)
            plot(ax1,x(~CurrentFrameFilter),y(~CurrentFrameFilter),'ow', 'MarkerSize', 10)
            hold(ax1, 'off')
            
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
                if ~isempty(app)
                    ax2 = app{2};
                else
                    ax2 = nucAxes;
                end
                imshow(Image,[],'Border','Tight', 'Parent', ax2,'InitialMagnification', 'fit')
                hold(ax2, 'on')
                PlotHandle=[];
                [NEllipses,~]=size(Ellipses{CurrentFrame});
                for j=1:NEllipses
                    PlotHandle=[PlotHandle,ellipse(Ellipses{CurrentFrame}(j,3),...
                        Ellipses{CurrentFrame}(j,4),...
                        Ellipses{CurrentFrame}(j,5),Ellipses{CurrentFrame}(j,1)+1,...
                        Ellipses{CurrentFrame}(j,2)+1, [],[],ax2)];
                        
                    text(ax2, Ellipses{CurrentFrame}(j,1)+1,...
                            Ellipses{CurrentFrame}(j,2)+1,num2str(j),'BackgroundColor',[.7 .9 .7]);

                end
                set(PlotHandle,'Color','r')
                hold(ax2, 'off')
                title(ax2, i)
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
                                            Particles{Channel}(PreviousFrameParticles(MinIndex(j))).Frame(end+1)=CurrentFrame;
                                            Particles{Channel}(PreviousFrameParticles(MinIndex(j))).Index(end+1)=ApprovedSpots(j);

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
                                    Particles{Channel}(PreviousFrameParticles).Frame(end+1)=CurrentFrame;
                                    Particles{Channel}(PreviousFrameParticles).Index(end+1)=MinIndex;
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
                %(MT, 2018-02-11) Added support for lattice imaging, 
                %maybe temporary - FIX LATER
                if strcmpi(ExperimentType,'1spot')||...
                        strcmpi(ExperimentType,'inputoutput')||...
                        strcmpi(ExperimentType,'2spot2color')||...
                        strcmpi(ExperimentType,'lattice')
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
    if isempty(app)
        close(ParticlesFig)
        if UseHistone
            close(NucleiFig)
        end
    end
else
    error('Experiment type in MovieDatabase not recognized')    
end

%If we only have one channel, then convert SpotFilter and Particles to a
%standard structure.
 
for currentChannel=1:NCh
    if ~isfield(Particles{currentChannel},'FrameApproved')
        for i=1:length(Particles{currentChannel})
            Particles{currentChannel}(i).FrameApproved=true(size(Particles{currentChannel}(i).Frame));
        end
    else
        for i=1:length(Particles{currentChannel})
            if isempty(Particles{currentChannel}(i).FrameApproved)
                Particles{currentChannel}(i).FrameApproved=true(size(Particles{currentChannel}(i).Frame));
            end
        end
    end
end

if NCh==1
   SpotFilter=SpotFilter{1};
   Particles=Particles{1};
end

mkdir([OutputFolder,filesep]);


save([OutputFolder,filesep,'Particles.mat'],'Particles','SpotFilter',...
    'Threshold1','Threshold2', '-v7.3');

% creating the field nc for FrameInfo
if exist([OutputFolder,filesep,'FrameInfo.mat'], 'file')
    numberOfFrames = length(FrameInfo);
    for currentFrame=1:numberOfFrames
        if currentFrame<nc9
            FrameInfo(currentFrame).nc=8;
        elseif (currentFrame>=nc9) &(currentFrame<nc10)
            FrameInfo(currentFrame).nc=9;
        elseif (currentFrame>=nc10)&(currentFrame<nc11)
            FrameInfo(currentFrame).nc=10;
        elseif (currentFrame>=nc11)&(currentFrame<=nc12)
            FrameInfo(currentFrame).nc=11;
        elseif (currentFrame>=nc12)&(currentFrame<=nc13)
            FrameInfo(currentFrame).nc=12;
        elseif (currentFrame>=nc13)&(currentFrame<=nc14)
            FrameInfo(currentFrame).nc=13;
        elseif currentFrame>=nc14
            FrameInfo(currentFrame).nc=14;
        end
    end
    save([OutputFolder,filesep,'FrameInfo.mat'],'FrameInfo')
else
    warning('Tried to save nc frame information, but could not since there is no FrameInfo.mat')
end



