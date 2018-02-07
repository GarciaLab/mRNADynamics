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
    Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder...
    ] = readMovieDatabase(PrefixOverrideFlag);

%Determine whether we're dealing with 2-photon data from Princeton or LSM
%data. 2-photon data uses TIF files. In LSM mode multiple files will be
%combined into one.
DTIF=dir([Folder,filesep,'*.tif']);
DLSM=dir([Folder,filesep,'*.lsm']);     %Zeiss confocal, old LSM format
DLIF=dir([Folder,filesep,'*.lif']);     %Leica confocal
DCZI=dir([Folder,filesep,'*.czi']);     %Zeiss confocal, new CZI format
DLAT=dir([Folder,filesep,'*_Settings.txt']);
DSPIN=dir([Folder,filesep,'*.nd']);     %Nikon spinning disk
DND2=dir([Folder,filesep,'*.nd2']);    %Nikon point scanner .nd2 files

if ~isempty(DTIF) && isempty(DLSM) && isempty(DCZI) && isempty(DSPIN)
    if isempty(DLIF)
        if isempty(DLAT)
            disp('2-photon @ Princeton data mode')
            D=DTIF;
            FileMode='TIF';
        else
            disp('Lattice Light Sheet data mode')
            D=DTIF;
            FileMode='LAT';
        end
    else
        disp('LIF export mode')
        D=DTIF;
        FileMode='LIFExport';
    end
elseif isempty(DTIF) && ~isempty(DLSM)
    disp('LSM mode')
    D=DLSM;
    FileMode='LSM';
elseif isempty(DTIF) && ~isempty(DCZI)
    disp('LSM (CZI) mode')
    D=DCZI;
    FileMode='LSM';
elseif ~isempty(DSPIN)
    disp('Nikon spinning disk mode with .nd files')
    D=dir([Folder,filesep,'*.tif']);     %spinning disk with .nd output files
    if isempty(D)
        error('No TIF files found')
    end
    FileMode = 'DSPIN';
elseif ~isempty(DND2)
    disp('Nikon LSM Mode');
    D=dir([Folder,filesep,'*.nd2']);
    FileMode='DND2';               
else
    error('File type not recognized. For LIF files, were they exported to TIF?')
end

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

    %Get the structure with the acquisition information
    ImageInfo = imfinfo([Folder,filesep,D(1).name]);
    
    %Do we have a second channel for Histone?
    if strcmp(Channel2,'His-RFP')
        HisChannel=1;
    else
        HisChannel=0;
    end

    
    %Get the flat-field information
    %Figure out the zoom factor
    try
        Zoom=ExtractInformationField(ImageInfo(1),'state.acq.zoomFactor=');
    catch
        error('Are you trying to use LIF mode but don''t have the .lif in your folder?')
    end
    %Look for the file
    FFDir=dir([Folder,filesep,'..',filesep,'*FF',Zoom(1:end-1),'x*.*']);



    %If there's more than one match then ask for help
    if length(FFDir)==1
        FFFile=FFDir(1).name;
    elseif isempty(FFDir)
        disp('Warning, no flat field file found. Press any key to proceed without it');
        FFImage=ones(ImageInfo(1).Height,ImageInfo(1).Width);
        pause
    else
        FFFile=uigetfile([Folder,filesep,'..',filesep,'FF',Zoom(1),'x*.*'],'Select flatfield file');
    end

    %If it's there, copy the image to the folder.
    if ~isempty(FFDir)
        FFImage=imread([Folder,filesep,'..',filesep,FFFile],1);
        imwrite(FFImage,[OutputFolder,filesep,Prefix,'_FF.tif']);
    end


    % ES 2013-10-30: ScanImage versions handle averaging differently.
    ScanImageVersionS = ExtractInformationField(ImageInfo(1), 'state.software.version=');
    
    %Figure out how many repeats of each image were taken, if averaging was
    %used and how many slices per stack we have
    NRepeats=str2num(ExtractInformationField(ImageInfo(1),'state.acq.numberOfFrames='));
    if strcmp(ScanImageVersionS(1:end-1), '3') % Referring to ScanImage 3.5.1
        AverageFlag=str2num(ExtractInformationField(ImageInfo(1),'state.acq.averaging='));
    elseif strcmp(ScanImageVersionS(1:end-1), '3.8') % Referring to ScanImage 3.8
        NumAvgFramesSave=str2num(ExtractInformationField(ImageInfo(1),'state.acq.numAvgFramesSave='));
        if NumAvgFramesSave > 1
            AverageFlag = 1;
        elseif NumAvgFramesSave == 1
            AverageFlag = 0;
        end
    end

    if AverageFlag
        NImages = numel(ImageInfo);
    else
        NImages = numel(ImageInfo)/NRepeats;
    end


    %Different indexes for the for-loops depending if we have multiple channels

    if ~HisChannel
        %If there is no second channel this is easy
        jChannel1=1:NRepeats:numel(ImageInfo);
        kChannel1=1:NRepeats;
        lChannel1=2:NRepeats;
    else
        %If there is a second channel for histone we will still use the
        %alignment information from the first channel.
        jChannel1=1:NRepeats*2:numel(ImageInfo);
        kChannel1=1:2:NRepeats*2;
        lChannel1=2:2:NRepeats*2;
    end



    %If the images were not averaged we will want to align them and average
    %them ourselves.

    if ~AverageFlag
        h=waitbar(0,'Finding shifts for individual images');
        for i=1:length(D)     
            waitbar(i/length(D),h)


            for j=jChannel1
                l=1;
                for k=kChannel1
                   Image(:,:,l)=imread([Folder,filesep,D(i).name],j+k-1);

                    %Make the images smaller to allow for the shifts
                    ImageSmall(:,:,l)=Image(MaxShift+1:(end-MaxShift),...
                        MaxShift+1:(end-MaxShift),l);
                    l=l+1;
                end

                for k=2:NRepeats
                    for xShift=-(MaxShift-1)/2:(MaxShift-1)/2
                        for yShift=-(MaxShift-1)/2:(MaxShift-1)/2
                            %Create a shifted image. Notice that we always
                            %align with respect to the first image. This should
                            %be fine but could be improved to "track" the
                            %shift.
                            clear ShiftedImageSmall
                            ShiftedImageSmall(:,:)=...
                                Image((1+MaxShift+yShift):(end-MaxShift+yShift),...
                                (1+MaxShift+xShift):(end-MaxShift+xShift),k);

                                cc=corrcoef(double(ImageSmall(:,:,1)),...
                                    double(ShiftedImageSmall));

                                CorrMatrix(yShift+(MaxShift-1)/2+1,...
                                    xShift+(MaxShift-1)/2+1)=cc(1,2);
                        end
                    end
                    [RowAlign ColAlign]=find(CorrMatrix==max(max(CorrMatrix)));

                    xShiftImages(i,j,k-1)=ColAlign-(MaxShift-1)/2-1;
                    yShiftImages(i,j,k-1)=RowAlign-(MaxShift-1)/2-1;
                end
            end


        end
        close(h)
    end

    if AverageFlag
        h=waitbar(0,'Copying to FISH folder');
    else
        h=waitbar(0,'Aligning images and copying to FISH folder');
    end

    
    %FrameInfo was defined above for the Leica and Zeiss modes. Here,
    %I'm going to clear it and go with the original definition. I might
    %have to go back to ExtractImageInformation and change a few things
    %though.
    clear FrameInfo
    
    for i=1:length(D)
        Suffix{i}=[iIndex(i,3),'_z??'];
        ImageInfo = imfinfo([Folder,filesep,D(i).name]);
        waitbar(i/length(D),h)
        
        
        
        
        FrameInfo(i)=ExtractImageInformation(ImageInfo(1));

    %If AverageFlag we just need to separate the images. Otherwise we'll do
    %the alignment and averaging.
        if AverageFlag
            error('Change the code to leave a blank image at the beginning and end of the stack')


            l=1;
            for j=jChannel1
                %Image=uint16(double(imread([Folder,filesep,D(i).name],j))./FFImage);
                Image=uint16(double(imread([Folder,filesep,D(i).name],j)));

                NewName=[Prefix,'_',iIndex(i,3),'_z',...
                    iIndex(l,2),'.tif'];
                imwrite(Image,[OutputFolder,filesep,NewName]);

                %Also copy the histone channel if it's present
                if HisChannel
                    ImageHistone(:,:,l)=imread([Folder,filesep,D(i).name],j+1);
                end

                l=l+1;

            end
            NewName=[Prefix,'-His_',iIndex(i,3),'.tif'];

            %Get the maximum projection and cap the highest brightness pixles
            MaxHistoneImage=max(ImageHistone,3);
            MaxHistoneImage(MaxHistoneImage>MaxHistone)=MaxHistone;

            imwrite(MaxHistoneImage,[OutputFolder,filesep,NewName]);
        else
            l=1;

            %Use the size information of this image to calculate create a
            %blank image for the beginning and end of the stack.
            BlankImage=uint16(zeros(ImageInfo(1).Height,ImageInfo(1).Width));
            %Save the first image
            NewName=[Prefix,'_',iIndex(i,3),'_z',...
                iIndex(1,2),'.tif'];
            imwrite(BlankImage,[OutputFolder,filesep,NewName],...
                'Compression','none');
            %Save the last image
            NewName=[Prefix,'_',iIndex(i,3),'_z',...
                iIndex(length(jChannel1)+2,2),'.tif'];
            imwrite(BlankImage,[OutputFolder,filesep,NewName],...
                'Compression','none');



            for j=jChannel1
                %Load the first image. We're not shifting this one
                ImageForAvg(:,:,1)=double(imread([Folder,filesep,D(i).name],j));

                %Now shift and average the actual images
                m=2;
                for k=kChannel1(2:end)
                    ImageForAvg(:,:,m)=double(imread([Folder,filesep,D(i).name],j+k-1));
                    ImageForAvg(:,:,m)=ShiftImage(ImageForAvg(:,:,m),...
                        xShiftImages(i,j,m-1),yShiftImages(i,j,m-1));
                    m=m+1;
                end
                Image=uint16(mean(ImageForAvg,3));  
                NewName=[Prefix,'_',iIndex(i,3),'_z',...
                    iIndex(l+1,2),'.tif'];
                imwrite(Image,[OutputFolder,filesep,NewName]);

                %Align the histone channel
                if HisChannel
                    ImageForAvg(:,:,1)=imread([Folder,filesep,D(i).name],j+1);
                    m=2;
                    for k=kChannel1(2:end)
                        ImageForAvg(:,:,m)=imread([Folder,filesep,D(i).name],j+k);
                        ImageForAvg(:,:,m)=ShiftImage(ImageForAvg(:,:,m),...
                            xShiftImages(i,j,m-1),yShiftImages(i,j,m-1));
                        m=m+1;
                    end
                    ImageHistone(:,:,l)=uint16(mean(ImageForAvg,3));

                end



                l=l+1;
            end

            if HisChannel
                ImageHistoneMax=uint16(max(ImageHistone,[],3));

                %Cap the highest brightness pixles
                ImageHistoneMax(ImageHistoneMax>MaxHistone)=MaxHistone;

                NewName=[Prefix,'-His_',iIndex(i,3),'.tif'];
                imwrite(ImageHistoneMax,[OutputFolder,filesep,NewName]);
            end
        end
    end
    
    

    %Get the actual time corresponding to each frame in seconds and add it to
    %FrameInfo
    for i=1:length(FrameInfo)
        FrameInfo(i).Time=etime(datevec(FrameInfo(i).TimeString),datevec(FrameInfo(1).TimeString));
    end
    
    
    %Add the information about the mode
    for i=1:length(FrameInfo)
        FrameInfo(i).FileMode='TIF';
    end
    
    
    close(h)

elseif strcmp(FileMode, 'LAT')
      
    %Do we have a second channel for Histone?
    if strcmp(Channel2,'His-RFP')
        HisChannel=1;
    else
        HisChannel=0;
    end
    
    %Figure out the different channels
    Channels={Channel1{1},Channel2{1}};

    %Coat protein channel
    MCPChannel=find((~cellfun(@isempty,strfind(lower(Channels),'mcp')))|...
        (~cellfun(@isempty,strfind(lower(Channels),'pcp')))|...
        (~cellfun(@isempty,strfind(lower(Channels),'lambda'))));
    if length(MCPChannel)>1
        error('Two coat proteins found. Should this be in 2spot2color mode?')
    elseif length(MCPChannel)==0    
        error('LAT Mode error: Channel name not recognized. Check MovieDatabase.XLSX')
    end
        
     
    %Histone channel
    HisChannel=find(~cellfun(@isempty,strfind(lower(Channels),'his')));
    %Distinguish between not having histone, but having a dummy channel
    if isempty(HisChannel)
        if find(~cellfun(@isempty,strfind(lower(Channels),'dummy')))
            HisChannel=0;%find(~cellfun(@isempty,strfind(lower(Channels),'dummy')));
        else
            HisChannel=0;
            display('Could not find a histone channel. Proceeding without it.')
        end
    end

    %Load the data
    im_stack = {};
    his_stack = [];
    mcp_stack = {};
    for j = 1:length(DTIF)
        fname = [Folder, filesep, DTIF(j).name];
        info = imfinfo(fname);
        num_images = numel(info);
        if HisChannel
            for i = 1:num_images
                im_stack{j, i} = imread(fname, i, 'Info', info);
                if ~isempty(strfind(fname, 'CamA'))
                    his_stack{j,i} = imread(fname, i, 'Info', info);
                    his_array(j,i, :, :) = imread(fname, i, 'Info', info);
                elseif ~isempty(strfind(fname, 'CamB'))
                    mcp_stack{j-size(his_stack, 1),i} = imread(fname, i, 'Info', info);
                else
                    error('Something is wrong with your channels. Please doublecheck moviedatabase')
                end 
            end
        else
            mcp_stack = im_stack;
        end
    end

    %Extract the metadata for each series
    NSeries=1; %Will always be true for lattice mode.
    NSlices=size(mcp_stack,2);
    NPlanes=numel(mcp_stack);

    %Number of channels %TO DO:  when we start using histone
    if ~isempty(his_stack) && ~isempty(mcp_stack)
        NChannels=2;
    else 
        NChannels=1;
    end

    %Finally, use this information to determine the number of frames in
    %each series
    NFrames=size(mcp_stack,1);

    %Get rid of the last frame as it is always incomplete because
    %that's when we stopped it
    NFrames=NFrames-1;
    NPlanes = NPlanes - NSlices;

    %Generate FrameInfo
    FrameInfo=struct('LinesPerFrame',{},'PixelsPerLine',{},...
        'NumberSlices',{},'ZStep',{},'FileMode',{},...
        'PixelSize',{});

    %Extract time information from text metadata file

    metaID = fopen([Folder, filesep, DLAT.name]);
    metastring = fscanf(metaID,'%s');
    tok = strsplit(metastring,{'Cycle(s):','Cycle(Hz'});
    timestep = str2double(tok{2});

    Frame_Times = zeros(1,NFrames*NSlices);
    for i = 1:length(Frame_Times)
        Frame_Times(i) = i*timestep;
    end

    %Get the time stamp corresponding to the first slice of each
    %Z-stack
    m = 1;
    for j=1:NPlanes/NFrames:NPlanes           
        InitialStackTime(m)=Frame_Times(j);
        m = m+1;
    end

    for i=1:NFrames
        FrameInfo(i).PixelsPerLine=256; %to do: need to parse this (ROI line in the metadata text file)
        FrameInfo(i).LinesPerFrame=512; % "
        %FrameInfo(i).PixelSize=str2num(LIFMeta.getPixelsPhysicalSizeX(1));
        FrameInfo(i).ZStep = .5; %to do: need to parse from metadata (but only if the metadata is correct)
        FrameInfo(i).NumberSlices=min(NSlices);
        FrameInfo(i).FileMode='LAT';
        FrameInfo(i).Time=InitialStackTime(i);
        FrameInfo(i).PixelSize = .1; %should parse this from metadata if it gets included
    end

    %Copy the data
    h=waitbar(0,'Extracting Lattice images');
    %Create a blank image
    BlankImage=uint16(zeros(size(mcp_stack{1,1})));
    %Now do His-RFP
    m = 1;
    if HisChannel 
        for j = 1:NFrames
            ims = squeeze(his_array(j, :, :, :));
            if strcmp(ProjectionType,'medianprojection')
                Projection=squeeze(median(ims,1));
            else
                Projection=squeeze(max(ims,[],1));
            end
            imwrite(uint16(Projection),...
            [OutputFolder,filesep,Prefix,'-His_',iIndex(m,3),'.tif']);
        m = m + 1;
        end
    end
    m=1;        %Counter for number of frames
    for j=1:NFrames
        %First do the MCP channel
        %Save the blank images at the beginning and end of the
        %stack
        NewName=[Prefix,'_',iIndex(j,3),'_z',iIndex(1,2),'.tif'];
        imwrite(BlankImage,[OutputFolder,filesep,NewName]);
        NewName=[Prefix,'_',iIndex(j,3),'_z',iIndex(NSlices+2,2),'.tif'];
        imwrite(BlankImage,[OutputFolder,filesep,NewName]);
        %Copy the rest of the images
        z = 2;
        for s = 1:NSlices
            NewName=[Prefix,'_',iIndex(j,3),'_z',iIndex(z,2),'.tif'];
            imwrite(mcp_stack{j,s},[OutputFolder,filesep,NewName]);
            z = z + 1;
        end
        m=m+1;
    end
    close(h)
                  
%LSM mode
elseif strcmp(FileMode,'LSM')
    
%     warning('Still need to add the FF information') NL: I think this
%     warning is out-dated
    
    %What type of experiment do we have?
    if strcmp(ExperimentType,'1spot') || strcmp(ExperimentType,'2spot') || strcmp(ExperimentType,'2spot1color')
    
        %Figure out the different channels
        if ~isempty(strfind(Channel1{1},'MCP'))
            coatChannel=1;
        elseif  strfind(Channel1{1},'His')
            histoneChannel=1;
        else
            error('LSM Mode error: Channel name not recognized. Check MovieDatabase.XLSX')
        end

        if ~isempty(strfind(Channel2{1},'MCP'))
            coatChannel=2;
        elseif  strfind(Channel2{1},'His')
            histoneChannel=2;
        else
            error('LSM Mode error: Channel name not recognized. Check MovieDatabase.XLSX')
        end
        fiducialChannel=histoneChannel;
        NSeries=length(D);
        Frame_Times=[];     %Store the frame information
        
        %m=1;        %Counter for number of frames
        clear m
        
        h=waitbar(0,'Extracting LSM images');
        for LSMIndex=1:NSeries
            waitbar(LSMIndex/NSeries,h);
            %Load the file using BioFormats
            LSMImages=bfopen([Folder,filesep,D(LSMIndex).name]);
            %Extract the metadata for each series
            LSMMeta = LSMImages{:, 4};      %OME Metadata
            LSMMeta2 = LSMImages{:, 2};     %Original Metadata
            
            %Figure out the number of slices in each series
            NSlices(LSMIndex)=str2num(LSMMeta.getPixelsSizeZ(0));
            %Number of channels
            NChannels(LSMIndex)=LSMMeta.getChannelCount(0);
            %Total number of planes acquired
            NPlanes(LSMIndex)=LSMMeta.getPlaneCount(0);
            %Finally, use this information to determine the number of frames in
            %each series
            NFrames(LSMIndex)=NPlanes(LSMIndex)/NSlices(LSMIndex)/NChannels(LSMIndex);
            
            %Check that the acquisition wasn't stopped before the end of a
            %cycle. If that is the case, some images in the last frame will
            %be blank and we need to remove them.
            if sum(sum(LSMImages{1}{end,1}))==0
                %Reduce the number of frames by one
                NFrames(LSMIndex)=NFrames(LSMIndex)-1;
                %Reduce the number of planes by NChannels*NSlices
                NPlanes(LSMIndex)=NPlanes(LSMIndex)-NChannels(LSMIndex)*NSlices(LSMIndex);
            end
            
            %First, get the starting time. This is not accessible in the
            %OME format, so we need to pull it out from the original
            %Metadata
            if NFrames(LSMIndex)<10
                NDigits=1;
            elseif NFrames(LSMIndex)<100
                NDigits=2;
            elseif NFrames(LSMIndex)<1000
                NDigits=3;
            else
                error('This program cannot support more than 1000 frames. This can be easily fixed')
            end
            
            %Get the starting time of this acquisition
            %This is different if I have an LSM or CZI file
            if ~isempty(DLSM)
                StartingTime(LSMIndex) = LSMMeta2.get(['TimeStamp #',iIndex(1,NDigits)]);
            elseif ~isempty(DCZI)
                TimeStampString=LSMMeta2.get('Global Information|Image|T|StartTime #1');
                TimeStampStrings{LSMIndex}=TimeStampString;
                %Get the number of days since 1/1/0000
                TimeStamp=datetime(TimeStampString(1:19),'InputFormat','yyyy-MM-dd''T''HH:mm:ss');
                %Get the milliseconds, note that we're not using them!
                MilliSeconds=TimeStamp(21:end-1);
                %Convert the time to seconds. We're ignoring the seconds.
                StartingTime(LSMIndex)=datenum(TimeStamp)*24*60*60;
            else
                error('Wrong file format. The code should not have gotten this far.')
            end
            
            
            %There seems to be some ambiguity with how information is being
            %pulled out of the metadata. This might be due to BioFormats
            %versioning. We'll check which case we have to provide
            %compatibility for both types of outcome. Whether we get the
            %"value" field or not is determined by the ValueField flag.
            ValueField=1;
            try
                LSMMeta.getPlaneDeltaT(0,0).value
            catch
                ValueField=0;
            end
            
            
            
            %Now get the time for each frame. We start the timer at the first time point
            %of the first data series.
            for j=0:(NSlices(LSMIndex)*NChannels(LSMIndex)):(NPlanes(LSMIndex)-1)
                if ~isempty(LSMMeta.getPlaneDeltaT(0,j))
                    if ValueField                    
                        Frame_Times=[Frame_Times,...
                            str2num(LSMMeta.getPlaneDeltaT(0,j).value)+...
                            StartingTime(LSMIndex)-StartingTime(1)];
                    else
                        Frame_Times=[Frame_Times,...
                            str2num(LSMMeta.getPlaneDeltaT(0,j))+...
                            StartingTime(LSMIndex)-StartingTime(1)];
                    end
                else  %If there's only one frame the DeltaT will be empty
                    Frame_Times=[Frame_Times,...
                        StartingTime(LSMIndex)-StartingTime(1)];
                end
            end

            %Save the information in FrameInfo
            if LSMIndex==1
                FrameRange=1:NFrames(LSMIndex);
            else
                FrameRange=(1:NFrames(LSMIndex))+length(FrameInfo);
            end
                
            for i=FrameRange
                FrameInfo(i).LinesPerFrame=str2double(LSMMeta.getPixelsSizeY(0));
                FrameInfo(i).PixelsPerLine=str2double(LSMMeta.getPixelsSizeX(0));
                FrameInfo(i).NumberSlices=min(NSlices);
                FrameInfo(i).FileMode='LSMExport';
                if ValueField
                    FrameInfo(i).PixelSize=str2num(LSMMeta.getPixelsPhysicalSizeX(0).value);
                    FrameInfo(i).ZStep=str2double(LSMMeta.getPixelsPhysicalSizeZ(0).value);
                else
                    FrameInfo(i).PixelSize=str2num(LSMMeta.getPixelsPhysicalSizeX(0));
                    FrameInfo(i).ZStep=str2double(LSMMeta.getPixelsPhysicalSizeZ(0));
                end
                FrameInfo(i).Time=Frame_Times(i);        %In seconds
            end
        
            %Copy the data
            %Create a blank image
            BlankImage=uint16(zeros(size(LSMImages{1}{1,1})));
            %Redefine LSMImages as a cell. I think we need this so that the
            %parfor loop doesn't freak out.
            LSMImages=LSMImages{1};
            %Size of images for generating blank images inside the loop
            Rows=size(LSMImages{1,1},1);
            Columns=size(LSMImages{1,1},2);
            
            
             for j=1:length(FrameRange)%NFrames(LSMIndex) 
                %First do the coat protein channel
                %Save the blank images at the beginning and end of the
                %stack
                NameSuffix=['_ch',iIndex(coatChannel,2)];
                
                NewName=[Prefix,'_',iIndex(FrameRange(j),3),'_z',iIndex(1,2),NameSuffix,'.tif'];
                imwrite(BlankImage,[OutputFolder,filesep,NewName]);
                NewName=[Prefix,'_',iIndex(FrameRange(j),3),'_z',iIndex(min(NSlices)+2,2),NameSuffix,'.tif'];
                imwrite(BlankImage,[OutputFolder,filesep,NewName]);
                %Copy the rest of the images
                n=1;        %Counter for slices
                for k=((j-1)*NSlices(LSMIndex)*NChannels(LSMIndex)+1+(coatChannel-1)):...
                        NChannels:(j*NSlices(LSMIndex))*NChannels
                    if n<=min(NSlices)
                        NewName=[Prefix,'_',iIndex(FrameRange(j),3),'_z',iIndex(n+1,2),NameSuffix,'.tif'];
                        imwrite(LSMImages{k,1},[OutputFolder,filesep,NewName]);
                        n=n+1;
                    end
                end

                %Now do His-RFP
                HisSlices=zeros([size(LSMImages{1,1},1),...
                    size(LSMImages{1,1},2),NSlices(LSMIndex)]);
                n=1;
                for k=((j-1)*NSlices(LSMIndex)*NChannels(LSMIndex)+1+(fiducialChannel-1)):...
                        NChannels(LSMIndex):(j*NSlices(LSMIndex))*NChannels(LSMIndex)
                    HisSlices(:,:,n)=LSMImages{k,1};
                    n=n+1;
                end
                Projection=median(HisSlices,3);
                imwrite(uint16(Projection),...
                            [OutputFolder,filesep,Prefix,'-His_',iIndex(FrameRange(j),3),'.tif']);
                %m=m+1;
            end
            
            
        end
        close(h)

        %Find the FF information
        
        %The FF can be in the folder with the data or in the folder
        %corresponding to the day.
        D1=dir([Folder,filesep,'FF*.lsm']);
        D1=[D1,dir([Folder,filesep,'FF*.czi'])];
        D2=dir([Folder,filesep,'..',filesep,'FF*.lsm']);
        D2=[D2,dir([Folder,filesep,'..',filesep,'FF*.czi'])];

        FFPaths={};
        for i=1:length(D1)
            FFPaths{end+1}=[Folder,filesep,D1(i).name];
        end
        for i=1:length(D2)
            FFPaths{end+1}=[Folder,filesep,'..',filesep,D2(i).name];
        end

        %Go through the FF files and see which one matches the pixel size
        %and image pixel number
        FFToUse=[];
        for i=1:length(FFPaths)
            LSMFF=bfopen(FFPaths{i});
            LSMFFMeta = LSMFF{:, 4};
            if (LSMFFMeta.getPixelsPhysicalSizeX(0)==LSMFFMeta.getPixelsPhysicalSizeX(0))&...
                    (LSMFFMeta.getPixelsSizeY(0)==LSMFFMeta.getPixelsSizeY(0))&...
                    (LSMFFMeta.getPixelsSizeX(0)==LSMFFMeta.getPixelsSizeX(0))
                FFToUse=[FFToUse,i];
            end
        end

        if length(FFToUse)> 1
            warning('Too many flat field images match the pixel and image size size')
            FFToUse = uigetfile('Select which flat field image to use');
        elseif length(FFToUse)==0
            warning('No flat field image found')
            clear LIFFF
        else
            LSMFF=bfopen(FFPaths{FFToUse});

            %Find the channel with the highest counts
            for i=1:size(LSMFF{1},1)
                MaxValue(i)=max(max(LSMFF{1}{i,1}));
            end
            [Dummy,ChannelToUse]=max(MaxValue);
            imwrite(LSMFF{1}{ChannelToUse,1},...
                [OutputFolder,filesep,Prefix,'_FF.tif']);
        end
    else
        error('Experiment type not supported for LSM format. Ask HG.')
    end
    
%LIFExport mode
elseif strcmp(FileMode,'LIFExport')
    
    %Extract time information from xml files
        XMLFolder=Folder;
        SeriesFiles = dir([XMLFolder,filesep,'*Series*Properties.xml']);
        if isempty(SeriesFiles)
            XMLFolder=[Folder,filesep,'MetaData'];
            SeriesFiles = dir([XMLFolder,filesep,'*Series*Properties.xml']);
            if isempty(SeriesFiles)
                error('XML MetaFiles could not be found. Did they get exported using the LAS software?')
            end
        end
 
        %Load the file using BioFormats
        %Figure out which one is not the FF
        LIFIndex=find(cellfun(@isempty,strfind({DLIF.name},'FF')));
        %Load the data, this might cause problems with really large sets
        LIFImages=bfopen([Folder,filesep,DLIF(LIFIndex).name]);
        %Extract the metadata for each series
        LIFMeta = LIFImages{:, 4};
        NSeries=LIFMeta.getImageCount()-1;                
        %Figure out the number of slices in each series
        NSlices = [];
        for i=1:NSeries
            NSlices(i)=str2num(LIFMeta.getPixelsSizeZ(i-1));
        end

        %Number of planes per series
        NPlanes = [];
        for i=1:NSeries
            NPlanes(i)=LIFMeta.getPlaneCount(i-1);
        end

        %Number of channels
        NChannels=LIFMeta.getChannelCount(0);

        %Finally, use this information to determine the number of frames in
        %each series        
        NFrames=NPlanes./NSlices/NChannels;
        %Get rid of the last frame as it is always incomplete because
        %that's when we stopped it
        NFrames=NFrames-1;
        NPlanes = NPlanes - NSlices*NChannels;      
        Frame_Times = zeros(1,sum(NFrames.*NSlices));
%         Frame_Times = [];
        m=1;
        for i = 1:NSeries
            xDoc = xmlread([XMLFolder,filesep,SeriesFiles(i).name]);
            TimeStampList = xDoc.getElementsByTagName('TimeStamp');
            for k = 0:(NFrames(i)*NSlices(i)*NChannels)-1
                TimeStamp = TimeStampList.item(k);
                Date = char(TimeStamp.getAttribute('Date'));
                Time = char(TimeStamp.getAttribute('Time'));
                Milli = char(TimeStamp.getAttribute('MiliSeconds'));
                time_in_days = datenum(strcat(Date,'-',Time,'-',Milli),'dd/mm/yyyy-HH:MM:SS AM-FFF');
                Frame_Times(m)=time_in_days*86400;
                m=m+1;
            end
        end
        
        First_Time=Frame_Times(1);
        for i = 1:length(Frame_Times)
            if Frame_Times(i)==0
                Frame_Times(i) = 0;
            else
                Frame_Times(i) = Frame_Times(i)-First_Time;
            end
        end

        %Get the time stamp corresponding to the first slice of each
        %Z-stack
        m=1;
        for i=1:NSeries
            if i==1
                StartIndex=1;
            else
                StartIndex=sum(NPlanes(1:i-1))+1;
            end
            for j=StartIndex:(NSlices(i).*NChannels):sum(NPlanes(1:i))
                InitialStackTime(m)=Frame_Times(j);
                m=m+1;
            end
        end

        for i=1:sum(NFrames)
            FrameInfo(i).LinesPerFrame=str2double(LIFMeta.getPixelsSizeY(0));
            FrameInfo(i).PixelsPerLine=str2double(LIFMeta.getPixelsSizeX(0));
            FrameInfo(i).NumberSlices=min(NSlices);
            FrameInfo(i).FileMode='LIFExport';
            %This is to allow for backwards compatibility with BioFormats
            if ~isempty(str2num(LIFMeta.getPixelsPhysicalSizeX(0)))
                FrameInfo(i).PixelSize=str2num(LIFMeta.getPixelsPhysicalSizeX(0));
                FrameInfo(i).ZStep=str2double(LIFMeta.getPixelsPhysicalSizeZ(0));
            else
                FrameInfo(i).PixelSize=str2num(LIFMeta.getPixelsPhysicalSizeX(0).value);
                FrameInfo(i).ZStep=str2double(LIFMeta.getPixelsPhysicalSizeZ(0).value);
            end
            FrameInfo(i).Time=InitialStackTime(i);
        end
        %Find the flat field (FF) information
        
        %The flat field image can be in the folder with the data or in the folder
        %corresponding to the date.
        D1=dir([Folder,filesep,'FF*.lif']);
        D2=dir([Folder,filesep,'..',filesep,'FF*.lif']);
        FFPaths={};
        for i=1:length(D1)
            FFPaths{end+1}=[Folder,filesep,D1(i).name];
        end
        for i=1:length(D2)
            FFPaths{end+1}=[Folder,filesep,'..',filesep,D2(i).name];
        end
        
        %Go through the FF files and see which one matches the pixel size
        %and image pixel number
        FFToUse=[];
        for i=1:length(FFPaths)
            LIFFF=bfopen(FFPaths{i});
            LIFFFMeta = LIFFF{:, 4};
            
            try
                %Check whether the number of pixels and pixel sizes match
                if (str2num(LIFFFMeta.getPixelsPhysicalSizeX(0).value)==...
                        str2num(LIFMeta.getPixelsPhysicalSizeX(0).value))&...
                        (str2num(LIFFFMeta.getPixelsSizeY(0))==...
                        str2num(LIFMeta.getPixelsSizeY(0)))&...
                        (str2num(LIFFFMeta.getPixelsSizeX(0))==...
                        str2num(LIFMeta.getPixelsSizeX(0)))
                    FFToUse=[FFToUse,i];
                %Sometimes, the number of pixels is right, but there's a slight
                %difference in the zoom factor. If the difference is less than
                %1%, then include it anyway
                elseif ~(str2num(LIFFFMeta.getPixelsPhysicalSizeX(0).value)==...
                        str2num(LIFMeta.getPixelsPhysicalSizeX(0).value))&...
                        (str2num(LIFFFMeta.getPixelsSizeY(0))==...
                        str2num(LIFMeta.getPixelsSizeY(0)))&...
                        (str2num(LIFFFMeta.getPixelsSizeX(0))==...
                        str2num(LIFMeta.getPixelsSizeX(0)))
                    if abs(1-str2num(LIFFFMeta.getPixelsPhysicalSizeX(0).value)/...
                            str2num(LIFMeta.getPixelsPhysicalSizeX(0).value))<0.01
                        warning('Same image size found for data and flat field, but a zoom difference smaller than 1%. Using FF image anyway.')
                        FFToUse=[FFToUse,i];
                    end
                end
            catch
                if (str2num(LIFFFMeta.getPixelsPhysicalSizeX(0))==...
                        str2num(LIFMeta.getPixelsPhysicalSizeX(0)))&...
                        (str2num(LIFFFMeta.getPixelsSizeY(0))==...
                        str2num(LIFMeta.getPixelsSizeY(0)))&...
                        (str2num(LIFFFMeta.getPixelsSizeX(0))==...
                        str2num(LIFMeta.getPixelsSizeX(0)))
                    FFToUse=[FFToUse,i];
                %Sometimes, the number of pixels is right, but there's a slight
                %difference in the zoom factor. If the difference is less than
                %1%, then include it anyway
                elseif ~(str2num(LIFFFMeta.getPixelsPhysicalSizeX(0))==...
                        str2num(LIFMeta.getPixelsPhysicalSizeX(0)))&...
                        (str2num(LIFFFMeta.getPixelsSizeY(0))==...
                        str2num(LIFMeta.getPixelsSizeY(0)))&...
                        (str2num(LIFFFMeta.getPixelsSizeX(0))==...
                        str2num(LIFMeta.getPixelsSizeX(0)))
                    if abs(1-str2num(LIFFFMeta.getPixelsPhysicalSizeX(0))/...
                            str2num(LIFMeta.getPixelsPhysicalSizeX(0)))<0.01
                        warning('Same image size found for data and flat field, but a zoom difference smaller than 1%. Using FF image anyway.')
                        FFToUse=[FFToUse,i];
                    end
                end
            end
        end 
        
        if length(FFToUse)> 1
            warning('Too many flat field images match the pixel and image size size')
            [FFFile,FFPath] =...
                uigetfile([Folder,filesep,'*.lif'],'Select which flat field image to use');
            LIFFF=bfopen([FFPath,FFFile]);
        elseif isempty(FFToUse)
            warning('No flat field image found')
            clear LIFFF
        else
            LIFFF=bfopen(FFPaths{FFToUse});
        end

        %If a flatfield image was found, process it
        if exist('LIFFF')
            %Find the channel with the highest counts
            for i=1:size(LIFFF{1},1)
                MaxValue(i)=max(max(LIFFF{1}{i,1}));
            end
            [~,ChannelToUse]=max(MaxValue);
            imwrite(LIFFF{1}{ChannelToUse,1},...
                [OutputFolder,filesep,Prefix,'_FF.tif']);
        end
        
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
                error('LIF Mode error: Channel name not recognized. Check MovieDatabase.XLSX')
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
            % MCP-mCherry as a fake histone channel
            if (~isempty(strfind(Channel1{1},'mCherry')))||(~isempty(strfind(Channel2{1},'mCherry')))
                if (~isempty(strfind(Channel1{1},'mCherry')))
                    fiducialChannel=1;
                    histoneChannel=1;
                elseif (~isempty(strfind(Channel2{1},'NLSmCherry')))
                    fiducialChannel=2;
                    histoneChannel=2;
                else
                    warning('mCherry channel not found. Cannot generate the fake nuclear image');
                end
            end
            
        elseif strcmpi(ExperimentType,'2spot2color')       %2 spots, 2 colors
            load('ReferenceHist.mat')
            fiducialChannel=0;
            histoneChannel=0;

            if (~isempty(strfind(Channel1{1},'mCherry')))||(~isempty(strfind(Channel2{1},'mCherry')))
                if (~isempty(strfind(Channel1{1},'mCherry')))
                    fiducialChannel=1;
                    histoneChannel=1;
                elseif (~isempty(strfind(Channel2{1},'NLSmCherry')))
                    fiducialChannel=2;
                    histoneChannel=2;
                else
                    warning('mCherry channel not found. Cannot generate the fake nuclear image');
                end
            end
        
        elseif strcmpi(ExperimentType,'input')        %Protein input mode
            %This mode assumes that at least one channel corresponds to the input.
            %It also check whether the second channel is histone. If there is
            %no histone channel it creates a fake channel using one of the
            %inputs.
            
            %Parse the information from the different channels
            Channels={Channel1{1},Channel2{1}};
            
            %We have no coat protein here.
            coatChannel=0;
            
            %Histone channel.
            histoneChannel=find(~cellfun(@isempty,strfind(lower(Channels),'his')));
            if isempty(histoneChannel)
                histoneChannel=0;
            else
                fiducialChannel=histoneChannel;
            end

            %Input channels
            inputProteinChannel=~cellfun(@isempty,Channels);
            if histoneChannel
                inputProteinChannel(histoneChannel)=0;
            else
                %If there was no histone channel, we need to choose which
                %input channel to use as our fiducial channel. We'll use
                %the first channel for now. We can try to be smarted about
                %this later on.
                warning('No histone channel found. Finding nuclei using the protein input channel.')
                fiducialChannel=1;                
            end
            inputProteinChannel=find(inputProteinChannel);
            
            %Save the information about the number of channels in FrameInfo
            for i=1:length(FrameInfo)
                FrameInfo(i).NChInput=length(inputProteinChannel);
            end
                        
%             
%             
%             if ~isempty(strfind(lower(Channel2{1}),'his'))
%                 fiducialChannel=2;
%                 inputProteinChannel=1;
%                 coatChannel=0;
%                 histoneChannel=2;
%             elseif ~isempty(strfind(lower(Channel1{1}),'his'))
%                 fiducialChannel=1;
%                 inputProteinChannel=2;
%                 coatChannel=0;
%                 histoneChannel=1;
%             else
%                 inputProteinChannel=1;
%                 fiducialChannel=1;      %We're assuming we can use the protein channel for segmentation
%                 coatChannel=0;
%                 histoneChannel=1;
%                 warning('No histone channel found. Finding nuclei using the protein input channel.')
%             end
%             
        elseif strcmpi(ExperimentType, 'inputoutput')
            if (~isempty(strfind(lower(Channel2),'mcp')))&...
                    ~isempty(strfind(lower(Channel2),'pcp'))
                coatChannel=2;
            elseif (~isempty(strfind(lower(Channel1),'mcp')))&...
                    ~isempty(strfind(lower(Channel1),'pcp'))
                coatChannel=1;
            else
                error('No MCP or PCP channel detected. Check MovieDatabase.XLSX')
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
            error('Experiment type not recognized. Check MovieDatabase.xlsx')
        end
        %Copy the data
        h=waitbar(0,'Extracting LIFExport images');
        %Create a blank image
        BlankImage=uint16(zeros(size(LIFImages{1}{1,1})));
        m=1;        %Counter for number of frames
        %Load the reference histogram for the fake histone channel
        load('ReferenceHist.mat')
        for i=1:NSeries
            waitbar(i/NSeries,h)
            for j=1:NFrames(i) 
                for q=1:NChannels
                    if (strcmpi(ExperimentType,'1spot') ||...
                            strcmp(ExperimentType,'2spot') ||...
                            strcmp(ExperimentType,'2spot1color')) && ...
                            q==coatChannel
                        %Save the blank images at the beginning and end of the
                        %stack
                        NameSuffix=['_ch',iIndex(q,2)];
                        NewName=[Prefix,'_',iIndex(m,3),'_z',iIndex(1,2),NameSuffix,'.tif'];
                        imwrite(BlankImage,[OutputFolder,filesep,NewName]);
                        NewName=[Prefix,'_',iIndex(m,3),'_z',iIndex(min(NSlices)+2,2),NameSuffix,'.tif'];
                        imwrite(BlankImage,[OutputFolder,filesep,NewName]);
                        %Copy the rest of the images
                        n=1;        %Counter for slices
                        firstImage = (j-1)*NSlices(i)*NChannels+1+(q-1);
                        lastImage = j*NSlices(i)*NChannels;
                        for k=firstImage:NChannels:lastImage
                            if n<=min(NSlices)
                                NewName=[Prefix,'_',iIndex(m,3),'_z',iIndex(n+1,2),NameSuffix,'.tif'];
                                   imwrite(LIFImages{i}{k,1},[OutputFolder,filesep,NewName]);
                                n=n+1;
                            end
                        end
                    elseif strcmpi(ExperimentType,'2spot2color')
                        NameSuffix=['_ch',iIndex(q,2)];

                        %Save the blank images at the beginning and end of the
                        %stack
                        NewName=[Prefix,'_',iIndex(m,3),'_z',iIndex(1,2),NameSuffix,'.tif'];
                        imwrite(BlankImage,[OutputFolder,filesep,NewName]);
                        NewName=[Prefix,'_',iIndex(m,3),'_z',iIndex(min(NSlices)+2,2),NameSuffix,'.tif'];
                        imwrite(BlankImage,[OutputFolder,filesep,NewName]);
                        %Copy the rest of the images
                        n=1;        %Counter for slices
                        firstImage = (j-1)*NSlices(i)*NChannels+1+(q-1);
                        lastImage = j*NSlices(i)*NChannels;
                        for k=firstImage:NChannels:lastImage
                            if n<=min(NSlices)
                                NewName=[Prefix,'_',iIndex(m,3),'_z',iIndex(n+1,2),NameSuffix,'.tif'];
                                   imwrite(LIFImages{i}{k,1},[OutputFolder,filesep,NewName]);
                                n=n+1;
                            end
                        end
                    %input-output mode
                    elseif strcmpi(ExperimentType, 'inputoutput')
                        %are we dealing with the coat channel?
                        if q==coatChannel
                            %Save the blank images at the beginning and end of the
                            %stack
                            NameSuffix=['_ch',iIndex(q,2)];
                            NewName=[Prefix,'_',iIndex(m,3),'_z',iIndex(1,2),NameSuffix,'.tif'];
                            imwrite(BlankImage,[OutputFolder,filesep,NewName]);
                            NewName=[Prefix,'_',iIndex(m,3),'_z',iIndex(min(NSlices)+2,2),NameSuffix,'.tif'];
                            imwrite(BlankImage,[OutputFolder,filesep,NewName]);
                            %Copy the rest of the images
                            n=1;        %Counter for slices
                            firstImage = (j-1)*NSlices(i)*NChannels+1+(q-1);
                            lastImage = j*NSlices(i)*NChannels;

                            TempNameSuffix = ['_ch',iIndex(q,2)];
                            for k=firstImage:NChannels:lastImage
                                if n<=min(NSlices)
                                    NewName=[Prefix,'_',iIndex(m,3),'_z',iIndex(n+1,2),TempNameSuffix,'.tif'];
                                       imwrite(LIFImages{i}{k,1},[OutputFolder,filesep,NewName]);
                                    n=n+1;
                                end
                            end
                            
                        %This is for the input channel    
                        else
                            %Save the blank images at the beginning and end of the
                            %stack
                            NameSuffix=['_ch',iIndex(q,2)];
                            NewName=[Prefix,'_',iIndex(m,3),'_z',iIndex(1,2),NameSuffix,'.tif'];
                            imwrite(BlankImage,[OutputFolder,filesep,NewName]);
                            NewName=[Prefix,'_',iIndex(m,3),'_z',iIndex(min(NSlices)+2,2),NameSuffix,'.tif'];
                            imwrite(BlankImage,[OutputFolder,filesep,NewName]);
                            %Copy the rest of the images
                            n=1;        %Counter for slices
                            firstImage = (j-1)*NSlices(i)*NChannels+1+(q-1);
                            lastImage = j*NSlices(i)*NChannels;

                            TempNameSuffix = ['_ch',iIndex(q,2)];
                            for k=firstImage:NChannels:lastImage
                                if n<=min(NSlices)
                                    NewName=[Prefix,'_',iIndex(m,3),'_z',iIndex(n+1,2),TempNameSuffix,'.tif'];
                                       imwrite(LIFImages{i}{k,1},[OutputFolder,filesep,NewName]);
                                    n=n+1;
                                end
                            end
                        end
                           
                    elseif strcmpi(ExperimentType, 'input')&&sum(q==inputProteinChannel)
                        %Are we dealing with one or two channels?
                        if length(inputProteinChannel)==1
                            NameSuffix=['_ch',iIndex(q,2)];
                        else
                            NameSuffix=['_ch',iIndex(q,2)];
                        end
                        %Save the blank images at the beginning and end of the
                        %stack
                        NewName=[Prefix,'_',iIndex(m,3),'_z',iIndex(1,2),NameSuffix,'.tif'];
                        imwrite(BlankImage,[OutputFolder,filesep,NewName]);
                        NewName=[Prefix,'_',iIndex(m,3),'_z',iIndex(min(NSlices)+2,2),NameSuffix,'.tif'];
                        imwrite(BlankImage,[OutputFolder,filesep,NewName]);
                        %Copy the rest of the images
                        n=1;        %Counter for slices
                        firstImage = (j-1)*NSlices(i)*NChannels+1+(q-1);
                        lastImage = j*NSlices(i)*NChannels;
                        for k=firstImage:NChannels:lastImage
                            if n<=min(NSlices)
                                NewName=[Prefix,'_',iIndex(m,3),'_z',iIndex(n+1,2),NameSuffix,'.tif'];
                                   imwrite(LIFImages{i}{k,1},[OutputFolder,filesep,NewName]);
                                n=n+1;
                            end
                        end             
                    end
                end
                %Now copy nuclear tracking images
                if fiducialChannel
                    HisSlices=zeros([size(LIFImages{i}{1,1},1),size(LIFImages{i}{1,1},2),NSlices(i)]);
                    otherSlices = HisSlices;
                    n=1;
                    firstImage = (j-1)*NSlices(i)*NChannels+1+(fiducialChannel-1);
                    lastImage = j*NSlices(i)*NChannels;
                    for k = firstImage:NChannels:lastImage
                        HisSlices(:,:,n)=LIFImages{i}{k,1};
                        otherSlices(:,:,n) = LIFImages{i}{k+1,1};
                        n=n+1;
                    end
                    if histoneChannel
                        if strcmp(ProjectionType,'medianprojection')
                            Projection=median(HisSlices,3);
                        else
                            Projection=max(HisSlices,[],3);
                        end
                        
                        %Think about the case when there is no His channel,
                        %and it is inputoutput mode or 1spot mode or 2spot2color.
                        %(MCP-mCherry)
                        if (isempty(strfind(Channel1{1}, 'His')))&&(isempty(strfind(Channel2{1}, 'His')))
                            if strcmpi(ExperimentType, 'inputoutput')|strcmpi(ExperimentType, '1spot')|strcmpi(ExperimentType,'2spot2color')|strcmpi(ExperimentType,'input')
                                if (~isempty(strfind(Channel1{1}, 'NLS')))|(~isempty(strfind(Channel2{1}, 'NLS')))
                                    %don't invert with NLS-MCP-mCherry
                                else
                                    %We don't want to use all slices. Only the center ones
                                    StackCenter=round((min(NSlices)-1)/2);
                                    StackRange=StackCenter-1:StackCenter+1;
                                    if strcmp(ProjectionType,'medianprojection')
                                        Projection=median(HisSlices(:,:,StackRange),[],3);
                                    else
                                        Projection=max(HisSlices(:,:,StackRange),[],3);
                                    end
                                    %invert images to make nuclei bright
                                    Projection=imcomplement(Projection);                            
                                end
                                Projection=histeq(mat2gray(Projection),ReferenceHist);
                                Projection = Projection*10000;
                            end
                        end
                    else 
                        %We don't want to use all slices. Only the center ones
                        StackCenter=round((min(NSlices)-1)/2);
                        StackRange=StackCenter-1:StackCenter+1;
                        if strcmp(ProjectionType,'medianprojection')
                            Projection=median(HisSlices(:,:,StackRange),3);
                            otherProjection=median(otherSlices(:,:,StackRange),3);
                        else
                            Projection=max(HisSlices(:,:,StackRange),[],3);
                            otherProjection=max(otherSlices(:,:,StackRange),[],3);
                        end
                        Projection = Projection + otherProjection;
                    end
 
                    imwrite(uint16(Projection),...
                    [OutputFolder,filesep,Prefix,'-His_',iIndex(m,3),'.tif']);

                    
                end
            m=m+1;
            end
        end
        close(h)
        
elseif strcmp(FileMode, 'LAT')
    [Output, FrameInfo] = ExportDataForFISH_Lattice(Prefix, D, Folder, OutputFolder, Channel1, Channel2, TAGOnly, ImageInfo); 
    
    
%Nikon spinning disk confocal mode - TH/CS 2017
elseif strcmp(FileMode,'DSPIN')||strcmp(FileMode,'DND2')
    
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
            error('LIF Mode error: Channel name not recognized. Check MovieDatabase.XLSX')
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
            error('No MCP or PCP channel detected. Check MovieDatabase.XLSX')
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
        error('Experiment type not recognized. Check MovieDatabase.xlsx')
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
        NSeries = length(DSPIN);
    else
        NSeries = length(DND2);
    end
    
    FrameTimes = cell(NSeries,1);
    TimeDiff = zeros(NSeries,1);
    TimeVector = [];
    TimeDiff2 = zeros(NSeries,1);

    %Deal with ordering of files if more than 9 series in the
    %acquisition (so 1, 2, .... 10 .  rather than 1, 10, 2...)
    fileName = cell(1,length(DSPIN));

    for ii = 1:NSeries
        if strcmp(FileMode, 'DSPIN')
            fileName{ii} = DSPIN(ii).name;
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
