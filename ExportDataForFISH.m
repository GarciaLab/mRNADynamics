function Prefix=ExportDataForFISH(varargin)

%This function grabs individual z-stacks and splits them in
%multiple channels so that it can be analyzed by the FISH code.
%It adds a blank image at the beginning and the end so that the code
%doesn't discard columns that peak at the edges of the Z-stack. 

%Options:
%medianprojection: Uses a median projection in the nuclear channel rather
%                  than the default maximum projection


%Note:
%Flatfield image: We assume there is no background of fluorescence / dark
%current.


%Where do things go?
%
%1) The original raw data: RawDynamics Folder
%2) The data processed such that it can be analyzed by the FISH code:
%   the PreProcessedData Folder
%3) The data analyzed by the FISH code: ProcessedData folder
%4) The resulting structures from the particle tracking:
%   DynamicsResults folder
%
%The idea of (4) being in Dropbox is that I don't need to be synchronizing
%the part related to the manual analysis.


%Parameters:
NIndices=3;     %Number of indices ScanImage used to save the files
MaxShift=9;     %Maximum shift in pixels corresponding to image shift and
                %alignment
MaxHistone=1000;    %Maximum intensity for the histone channel. Anything above
                    %this will be capped.
ProjectionType = 'maxprojection'; %Default setting for z-projection is maximum-based...
                    %This may fail when high intensity reflections are present
                

%Look at parameters
PrefixOverrideFlag = 0;
SkipFrames=[];
k=1;
while k<=length(varargin)
    if strcmpi(varargin{k},'skipframes')
        SkipFrames=varargin{k+1};
        k=k+1;
        warning('SkipFrame mode.')
    elseif strcmpi(varargin{k},'medianprojection')
        ProjectionType = 'medianprojection';
    else
        Prefix = varargin{k};
        PrefixOverrideFlag = 1;
    end
    k=k+1;
end

[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
    Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder, Channel3...
    ] = readMovieDatabase(PrefixOverrideFlag);

[D, FileMode] = DetermineFileMode(Folder);

%Create the output folder
OutputFolder=[PreProcPath,filesep,Prefix];
mkdir(OutputFolder)

%Generate FrameInfo
FrameInfo=struct('LinesPerFrame',{},'PixelsPerLine',{},...
    'NumberSlices',{},'ZStep',{},'FileMode',{},...
    'PixelSize',{});


%Extract channel information
%This information will be stored in FrameInfo for use by subsequent parts
%of the code. Note, however, that the channels are also extracted in this
%code for each data type. I should integrate this.




if strcmp(FileMode,'TIF') && ~strcmp(FileMode,'DSPIN')

  FrameInfo = process2PhotonPrincetonData(Folder, D, FrameInfo, Channel2, OutputFolder);

elseif strcmp(FileMode, 'LAT')
      
  FrameInfo = processLatticeLightSheetData(Folder, D, Channel1, Channel2, ProjectionType, Prefix, OutputFolder);
                  
%LSM mode
elseif strcmp(FileMode,'LSM')
    
  FrameInfo = processZeissConfocalLSMData(Folder, D, ExperimentType, Channel1, Channel2, Prefix, OutputFolder);
    
%LIFExport mode
elseif strcmp(FileMode,'LIFExport')
    
  FrameInfo = processLIFExportMode(Folder, ExperimentType, ProjectionType, Channel1, Channel2, Channel3, Prefix, OutputFolder);        
    
%Nikon spinning disk confocal mode - TH/CS 2017
elseif strcmp(FileMode,'DSPIN')||strcmp(FileMode,'DND2')

    SPINDir=dir([Folder,filesep,'*.nd']);     %Nikon spinning disk
    ND2Dir=dir([Folder,filesep,'*.nd2']);    %Nikon point scanner .nd2 files
    
    % Excel-reading part of the LIF acquisition. In previous versions this comes after the bfopen stuff, but I don't think it
    % matters and the way I have it written, as the script loops through all relevant acquisitions it not only writes the FrameInfo struc,
    % but also writes the images. The image writing requires coatChannel and fiducialChannel
    % TH this strategy may also make it hard to assimilate this code in with HG's all powerful one-script-to-rule-them-all
    if strcmpi(ExperimentType,'1spot') || strcmp(ExperimentType,'2spot') || strcmp(ExperimentType,'2spot1color')
        %Figure out the different channels
        Channels={Channel1{1},Channel2{1}};
        %Coat protein channel
        coatChannel=find((~cellfun(@isempty,strfind(lower(Channels),'mcp')))|...
            (~cellfun(@isempty,strfind(lower(Channels),'pcp')))|...
            (~cellfun(@isempty,strfind(lower(Channels),'lambda'))));
        if length(coatChannel)>1
            error('Two coat proteins found. Should this be in 2spot2color mode?')
        elseif isempty(coatChannel)
            error('LIF Mode error: Channel name not recognized. Check MovieDatabase')
        end
        
        %Histone channel
        histoneChannel=find(~cellfun(@isempty,strfind(lower(Channels),'his')));
        fiducialChannel = histoneChannel;
        %Distinguish between not having histone, but having a dummy channel
        if isempty(fiducialChannel)
            if find(~cellfun(@isempty,strfind(lower(Channels),'dummy')))
                fiducialChannel=0;
            else
                fiducialChannel=0;
                display('Could not find a histone channel. Proceeding without it.')
            end
        end
        
    elseif strcmpi(ExperimentType,'2spot2color')       %2 spots, 2 colors
        load('ReferenceHist.mat')
        if (~isempty(strfind(Channel1{1},'mCherry')))|(~isempty(strfind(Channel2{1},'mCherry')))
            if (~isempty(strfind(Channel1{1},'mCherry')))
                fiducialChannel=1;
            elseif (~isempty(strfind(Channel2{1},'mCherry')))
                fiducialChannel=2;
            else
                error('mCherry channel not found. Cannot generate the fake nuclear image')
            end
        end
        
    elseif strcmpi(ExperimentType,'input')        %Protein input mode
        %This mode assumes that one channel corresponds to the input and
        %that there is no second channel or, at the most, there is a
        %histone channel
        if ~isempty(strfind(lower(Channel2{1}),'his'))
            fiducialChannel=2;
            inputProteinChannel=1;
        elseif ~isempty(strfind(lower(Channel1{1}),'his'))
            fiducialChannel=1;
            inputProteinChannel=2;
        else
            fiducialChannel=0;
            inputProteinChannel=1;
            warning('No histone channel found. Finding nuclei using the protein input channel.')
        end
        
    elseif strcmpi(ExperimentType, 'inputoutput')
        if (~isempty(strfind(lower(Channel2),'mcp')))&...
                ~isempty(strfind(lower(Channel2),'pcp'))
            coatChannel=2;
        elseif (~isempty(strfind(lower(Channel1),'mcp')))&...
                ~isempty(strfind(lower(Channel1),'pcp'))
            coatChannel=1;
        else
            error('No MCP or PCP channel detected. Check MovieDatabase')
        end
        
        if (~isempty(strfind(Channel1{1},'mCherry')))|(~isempty(strfind(Channel2{1},'mCherry')))
            if (~isempty(strfind(Channel1{1},'mCherry')))
                fiducialChannel=1;
                histoneChannel=1;
            elseif (~isempty(strfind(Channel2{1},'mCherry')))
                fiducialChannel=2;
                histoneChannel=2;
            else
                error('mCherry channel not found. Cannot generate the fake nuclear image')
            end
        end
    else
        error('Experiment type not recognized. Check MovieDatabase')
    end
    %%Thus endeth the excel file reading section.
    
    
    
    %Get the flat-field information
    FFDir = dir([SourcePath, filesep, Prefix(1:10), filesep,  'FF*.tif']);
    
    if length(FFDir)==1
        DSPINFF = bfopen([SourcePath, filesep, Prefix(1:10), filesep, FFDir(1).name]);
    elseif isempty(FFDir)
        display('Error: no flat field file found, proceeding with it.');
    else
        DSPINFF = uigetfile([SourcePath, filesep, Prefix(1:10)],'Select flatfield file');
    end
    
    if exist('DSPINFF')
        FF = DSPINFF{1}{1,1};
    
        %Find the channel with the highest counts (CS this will be the coat protein channel that we care about)
        for iFrame=1:size(DSPINFF{1},1)
            MaxValue(iFrame)=max(max(DSPINFF{1}{iFrame,1}));
        end
        [~,ChannelToUse]=max(MaxValue);
        imwrite(DSPINFF{1}{ChannelToUse,1},[OutputFolder,filesep,Prefix,'_FF.tif']);
    end
        
    if strcmp(FileMode,'DSPIN')  
        NSeries = length(SPINDir);
    else
        NSeries = length(ND2Dir);
    end
    
    FrameTimes = cell(NSeries,1);
    TimeDiff = zeros(NSeries,1);
    TimeVector = [];
    TimeDiff2 = zeros(NSeries,1);

    %Deal with ordering of files if more than 9 series in the
    %acquisition (so 1, 2, .... 10 .  rather than 1, 10, 2...)
    fileName = cell(1,length(SPINDir));

    for ii = 1:NSeries
        if strcmp(FileMode, 'DSPIN')
            fileName{ii} = SPINDir(ii).name;
        else
            fileName{ii} = D(ii).name;
        end
    end
    fileNameSort = natsort(fileName);          %for DSPIN data: need to have natsort function downloaded from matlab central fileexchange and in mRNADynamics folder
    
    %This first for loop is to load the data & metadata from each series and to generate the overall time vector
    for iSeries = 1:NSeries

        %Load the data
        DSPINdata(iSeries,:) = bfopen([Folder, filesep, char(fileNameSort(iSeries))]);
        DSPINimages{iSeries,:} = DSPINdata{iSeries,1}(:,1);
        DSPINmeta(iSeries) = DSPINdata{iSeries,4};                             % load the OME metadata
        DSPINmeta2(iSeries) = DSPINdata{iSeries,2};           % metadata useful for time vector

        NSlices = str2double(DSPINmeta(iSeries).getPixelsSizeZ(0));      % # of z slices
        NChannels = DSPINmeta(iSeries).getChannelCount(0);
        NPlanes(iSeries) = DSPINmeta(iSeries).getPlaneCount(0);                   % # planes (individual files in this series) for calculating #time points
        NTimePoints(iSeries) = NPlanes(iSeries)/NSlices/NChannels;

        
        %Generate the time vector. This has to be done slightly differently
        %for DSPIN and DND2 data.
        
        %Get the time-from-start-of-SERIES for each timepoint in the series
        FrameIncrement = NSlices.*NChannels;   %increment to step through to get the time-from-start for each timepoint
        TimeFromStart = zeros(1,NTimePoints(iSeries));
        
        for iTP = 2:NTimePoints(iSeries)
            try
                TimeFromStart(iTP) = str2double(DSPINmeta(iSeries).getPlaneDeltaT(0,(FrameIncrement*(iTP-1))).value);
            catch
                TimeFromStart(iTP) = str2double(DSPINmeta(iSeries).getPlaneDeltaT(0,(FrameIncrement*(iTP-1))));
            end
        end

        FrameTimes{iSeries} = TimeFromStart./60;       %time from start of iSeries that first z of channel1 at each timepoint is taken, in minutes

        %Get global start time of the series
        if strcmp(FileMode,'DSPIN')
            temp = DSPINmeta2(iSeries).get(['timestamp #', iIndex(1,1)]);   %to get global start time from .nd metadata
            try
                temp = datevec(temp, 'yyyymmdd HH:MM:SS.FFF');
            catch
                temp = datevec(temp, 'yyyy-mm-ddTHH:MM:SS.FFF');
            end
        else 
            temp = DSPINmeta2(iSeries).get('Global dTimeAbsolute');  %to get global start time from .nd2 file metadata
            temp = datetime(temp,'convertfrom','juliandate'); %gets the date right and the dt between acq'ns is const.
            temp = datevec(temp);
        end
        SeriesStartTime(iSeries,:) = temp;       

        %Find the difference in global time between it and the previous series
        %(except for the first series, which has no previous one) 
        if iSeries==1
            TimeDiff2 = 0;
        elseif iSeries>1
            TimeDiff(iSeries) = (etime(SeriesStartTime(iSeries,:), SeriesStartTime(iSeries-1,:)))./60;    %time diff (mins) between start of one series and the next
            TimeDiff2(iSeries) = TimeDiff(iSeries)+TimeDiff2(iSeries-1);        %converted into time since start of movie
        end

        %Add the difference in global time to the time-from-start-of-series in
        %each case
        temp2 = FrameTimes{iSeries} + TimeDiff2(iSeries);
        TimeVector = [TimeVector, temp2];     %This is the full time vector of time since the start of the movie for each timepoint

        clear TimeFromStart temp temp2        

    end
            
    %This second for loop is to write the images to new files in the PreprocessedData folder
    disp('Extracting Images')
    %Create a blank image

    BlankImage=uint16(zeros(size(DSPINdata{1}{1,1})));       
    tp = 1; %counter for timepoints to enable daisy-chaining of timepoint numbers across all series

    for iSeries = 1:NSeries

        disp(['Series ', num2str(iSeries)])

        for iTP = 1:NTimePoints(iSeries)             %for each timepoint in this series
            
            %For the coat protein channel, save  blank images at the beginning and end of the stack (i.e. slice 1 and slice NSlices+2)
            NameSuffix=['_ch',num2str(0),num2str(1)];       %CS 20171119: this is a quick fix
            NewName=[Prefix,'_',iIndex((iTP+(tp-1)),3),'_z',iIndex(1,2),NameSuffix,'.tif'];                % this business with (tp-1) is to daisy chain the timepoint numbering of the images across all series         
            imwrite(BlankImage,[OutputFolder,filesep,NewName]);
            NewName=[Prefix,'_',iIndex((iTP+(tp-1)),3),'_z',iIndex(NSlices+2,2),NameSuffix,'.tif'];
            imwrite(BlankImage,[OutputFolder,filesep,NewName]);

             %Copy the coat protein images into the intervening-numbered files
            n=1;        %Counter for slices

            if strcmp(FileMode,'DSPIN')  % .nd files list images as tp1: Ch1 z1:21, Ch2 z1:21 before moving to tp2
                firstImage = (iTP-1)*FrameIncrement+1;      %index of first z slice in the timepoint to find the file
                lastImage = (iTP-1)*FrameIncrement+NSlices; %index of last z slice in the timepoint to find the file
                for ii = firstImage:lastImage
                    NewName=[Prefix,'_',iIndex(iTP+(tp-1),3),'_z',iIndex(n+1,2),NameSuffix,'.tif'];      %n+1 to not overwrite initial blank image
                    imwrite(DSPINimages{iSeries}{ii},[OutputFolder,filesep,NewName]);
                    n = n+1;
                end
            else   %for DND2 mode: .nd2 files list images as tp1: z1 Ch1,Ch2 z2 Ch1,Ch2... before moving to tp2
                firstImage = (iTP-1)*FrameIncrement+1;
                lastImage = (iTP-1)*FrameIncrement+2*NSlices;
                for ii = firstImage:2:lastImage
                    NewName=[Prefix,'_',iIndex(iTP+(tp-1),3),'_z',iIndex(n+1,2),NameSuffix,'.tif'];      %n+1 to not overwrite initial blank image
                    imwrite(DSPINimages{iSeries}{ii},[OutputFolder,filesep,NewName]);
                    n = n+1;
                end
            end

            %Get the His slices for this timepoint
            m = 1;      %Counter for slices
            
            if strcmp(FileMode,'DSPIN')
                firstHisImage = (iTP-1)*FrameIncrement+(NSlices+1);
                lastHisImage = (iTP-1)*FrameIncrement+(2*NSlices);
                
                for iSlice = firstHisImage:lastHisImage                    %takes each slice of the current His channel and puts it in HisSlices
                    HisSlices(:,:,m)=DSPINimages{iSeries}{iSlice};
                    m = m+1;
                end
            else
                firstHisImage = (iTP-1)*FrameIncrement+2;
                lastHisImage = (iTP-1)*FrameIncrement+2*NSlices;
                for iSlice = firstHisImage:2:lastHisImage                    %takes each slice of the current His channel and puts it in HisSlices
                    HisSlices(:,:,m)=DSPINimages{iSeries}{iSlice};
                    m = m+1;
                end
            end

            %Make a projection of the stack of images for the timepoint
            Projection=max(HisSlices,[],3);            %CS: use max projection to improve segmentation unless encounter issue below.
            %Projection=median(HisSlices,3);           %median projection avoids issue of occasional bright reflective top/bottom slices.

            %There's some bug in BioFormats for Windows that loads some
            %histone frames as blank. I have no clue where this is coming
            %from, but we should send a report to them. In the meantime, if
            %this happens, we'll use the histone from the previous frame.
            if max(max(Projection))==0
                warning(['Could not load Histone channel in frame ',...
                    num2str(iTP),' of series ',num2str(iSeries),...
                    '. Proceeding by copying the histone image from the previous frame. Send a report of this bug to BioFormats!'])
                %Get the His slices for the previous timepoint
                m = 1;      %Counter for slices
                firstHisImage = (iTP-2)*FrameIncrement+(NSlices+1);
                lastHisImage = (iTP-2)*FrameIncrement+(2*NSlices);

                for iSlice = firstHisImage:lastHisImage                    %takes each slice of the current His channel and puts it in HisSlices
                    HisSlices(:,:,m)=DSPINimages{iSeries}{iSlice};
                    m = m+1;
                end
                
                Projection=max(HisSlices,[],3);            %CS: use max projection to improve segmentation unless encounter issue below.
            end
            
            %Write the projection to a file
            NewName = [Prefix,'-His_',iIndex(iTP+(tp-1),3),'.tif'];
            imwrite(uint16(Projection),[OutputFolder,filesep,NewName]);

        end

        tp = tp + NTimePoints(iSeries);

    end  
        
        %Save relevant information in FrameInfo structure. Atm the
        %naming conventions are slightly different. Fix later.
        for iFrame=1:length(TimeVector)
            FrameInfo(iFrame).LinesPerFrame = str2double(DSPINmeta(1).getPixelsSizeY(0));
            FrameInfo(iFrame).PixelsPerLine = str2double(DSPINmeta(1).getPixelsSizeX(0));
            FrameInfo(iFrame).NumberSlices = NSlices;
            FrameInfo(iFrame).FileMode = 'DSPIN';
            try
                FrameInfo(iFrame).PixelSize = str2num(DSPINmeta(1).getPixelsPhysicalSizeX(0).value);
                FrameInfo(iFrame).ZStep = str2double(DSPINmeta(1).getPixelsPhysicalSizeZ(0).value);
            catch
                FrameInfo(iFrame).PixelSize = str2num(DSPINmeta(1).getPixelsPhysicalSizeX(0));
                FrameInfo(iFrame).ZStep = str2double(DSPINmeta(1).getPixelsPhysicalSizeZ(0));
            end
            %HG: This should be in minutes. I made the appropriate change.
            FrameInfo(iFrame).Time = TimeVector(iFrame)*60;      %CS: this reports time from start of first acquisition
        end
        
        % Extract parameters needed later in addParticlePosition
        % Find position of the zoom image (movie) in the scope ref frame. Just use the first image file to get these parameters
        zPositionX = str2double(DSPINmeta(1).getPlanePositionY(0,0));      %TH these are in um, are are reversed for some reason
        zPositionY = str2double(DSPINmeta(1).getPlanePositionX(0,0));
        
        % Pixel size in um
        try
            PixelSizeZoom=str2double(DSPINmeta(1).getPixelsPhysicalSizeX(0).value); %may need to try .getValue instead
        catch
            PixelSizeZoom=str2double(DSPINmeta(1).getPixelsPhysicalSizeX(0));
        end
        
        % Size of the zoom images (movie) in pix
        zImageSizeXpix = str2double(DSPINmeta(1).getPixelsSizeX(0));
        zImageSizeYpix = str2double(DSPINmeta(1).getPixelsSizeY(0));
        
        %Save extra metadata info use later on
        mkdir([DropboxFolder,filesep,Prefix])
        save([DropboxFolder,filesep,Prefix,filesep,'bfMetaData.mat'],...
            'zPositionX','zPositionY','PixelSizeZoom','zImageSizeXpix','zImageSizeYpix')
  
end

%Skipping frames?
if ~isempty(SkipFrames)
    %Filter FrameInfo
    FrameFilter=ones(size(FrameInfo));
    FrameFilter(SkipFrames)=0;
    FrameInfo=FrameInfo(logical(FrameFilter));
    
    %Rename all Histone channel images
    %Find all the Histone files
    D=dir([OutputFolder,filesep,'*-His*.tif']);
    %Delete the skipped files
    for i=SkipFrames
        delete([OutputFolder,filesep,D(i).name]);
    end
    %Rename all remaining files
    D=dir([OutputFolder,filesep,'*-His*.tif']);
    for i=1:length(D)
        if ~strcmp([OutputFolder,filesep,D(i).name],...
                [OutputFolder,filesep,D(i).name(1:end-7),iIndex(i,3),'.tif'])
            movefile([OutputFolder,filesep,D(i).name],...
                [OutputFolder,filesep,D(i).name(1:end-7),iIndex(i,3),'.tif'])
        end
    end
    
    %Rename all coat protein channel images
    %Find all the coat protein files
    D=dir([OutputFolder,filesep,'*_z01.tif']);
    %Delete the skipped files
    for i=SkipFrames
        D2=dir([OutputFolder,filesep,D(i).name(1:end-6),'*.tif']);
        for j=1:length(D2)
            delete([OutputFolder,filesep,D2(j).name])
        end
    end
    %Rename all remaining files
    D=dir([OutputFolder,filesep,'*_z01.tif']);
    for i=1:length(D)
        if ~strcmp([OutputFolder,filesep,D(i).name],...
                    [OutputFolder,filesep,D(i).name(1:end-11),iIndex(i,3),'_z01.tif'])
            D2=dir([OutputFolder,filesep,D(i).name(1:end-6),'*.tif']);
            for j=1:length(D2)
                movefile([OutputFolder,filesep,D2(j).name],...
                    [OutputFolder,filesep,D2(j).name(1:end-11),iIndex(i,3),D2(j).name(end-7:end)])
            end
        end
    end
end


%Save the information about the various frames
mkdir([DropboxFolder,filesep,Prefix])
save([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'],...
    'FrameInfo')
