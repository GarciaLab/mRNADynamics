function Prefix=ExportDataForFISH(varargin)

%This function grabs individual stacks taken by ScanImage and splits them in
%multiple channels so that it can be analyzed by the FISH code.
%It adds a blank image at the beginning and the end so that the code
%doesn't discard columns that peak at the edges of the Z-stack. It also
%generates a TAG file with lots of channels. One for each time point.

%Options:
%TAGOnly: Generates only the tag file


%Note:
%Flatfield image: We assume there is no background of fluorescence / dark
%current. I should check this. Maybe do a dark count reading for every time
%I take data?


%Where do things go?
%
%1) The original raw data: D:\Hernan\LivemRNA\Analysis
%2) The data processed such that it can be analyzed by the FISH code:
%   D:\Hernan\FISHDrosophila\Data
%3) The data analyzed by the FISH code: D:\Hernan\FISHDrosophila\Analysis
%4) The resulting structures from the particle tracking:
%   D:\Hernan\Dropbox\LivemRNAData
%
%The idea of (4) being in Dropbox is that I don't need to be synchronizing
%the part related to the manual analysis.


%Parameters:
NIndices=3;     %Number of indices ScanImage used to save the files
MaxShift=9;     %Maximum shift in pixels corresponding to image shift and
                %alignment
MaxHistone=1000;    %Maximum intensity for the histone channel. Anything above
                    %this will be capped.
                

%Look at parameters
PrefixOverrideFlag = 0;
TAGOnly=0;
SkipFrames=[];
k=1;
while k<length(varargin)
    if strcmp(lower(varargin{k}),'tagonly')
            TAGOnly=1;
    elseif strcmp(lower(varargin{k}),'skipframes')
        SkipFrames=varargin{k+1};
        k=k+1;
        warning('SkipFrame mode.')
    else
        Prefix = varargin{k};
        PrefixOverrideFlag = 1;
    end
    k=k+1;
end

                
%Figure out the initial folders. We'll update the Drobpox one later on in the code.
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,SchnitzcellsFolder]=...
    DetermineLocalFolders;
                    
%Get the folder with the data
if ~PrefixOverrideFlag
    Folder=uigetdir(SourcePath,'Select folder with data');
else
    HyphenPositionR = find(Prefix == '-');
    DateFolderS = Prefix(1 : HyphenPositionR(3)-1);
    LineFolderS = Prefix(HyphenPositionR(3)+1 : end);
    Folder = [SourcePath, filesep, DateFolderS, filesep, LineFolderS];
end

%Determine whether we're dealing with 2-photon data from Princeton or LSM
%data. 2-photon data uses TIF files. In LSM mode multiple files will be
%combined into one.
DTIF=dir([Folder,filesep,'*.tif']);
DLSM=dir([Folder,filesep,'*.lsm']);
DLIF=dir([Folder,filesep,'*.lif']);

if (length(DTIF)>0)&(length(DLSM)==0)
    if length(DLIF)==0
        display('2-photon @ Princeton data mode')
        D=DTIF;
        FileMode='TIF';
    else
        display('LIF export mode')
        D=DTIF;
        FileMode='LIFExport';
    end
elseif (length(DTIF)==0)&(length(DLSM)>0)
    display('LSM mode')
    D=DLSM;
    FileMode='LSM';
else
    error('File type not recognized')
end


%Get the information from the last two folders in the structure
if ~PrefixOverrideFlag
    SlashPositions=strfind(Folder,filesep);
    
    Prefix=[Folder((SlashPositions(end-1)+1):(SlashPositions(end)-1)),'-',...
        Folder((SlashPositions(end)+1):(end))];
end

%Figure out what type of experiment we have
[XLSNum,XLSTxt]=xlsread([DropboxFolder,filesep,'MovieDatabase.xlsx']);
DataFolderColumn=find(strcmp(XLSTxt(1,:),'DataFolder'));
ExperimentTypeColumn=find(strcmp(XLSTxt(1,:),'ExperimentType'));
Channel1Column=find(strcmp(XLSTxt(1,:),'Channel1'));
Channel2Column=find(strcmp(XLSTxt(1,:),'Channel2'));

% Convert the prefix into the string used in the XLS file
Dashes = strfind(Prefix, '-');
PrefixRow = find(strcmp(XLSTxt(:, DataFolderColumn),...
    [Prefix(1:Dashes(3)-1), filesep, Prefix(Dashes(3)+1:end)]));
% ES 2013-10-06: Removing the hard-coding for selecting the date string in
% 'Prefix'. This is because I tend to put a letter after the date:
% '2013-10-06A', for instance, instead of '2013-10-06'. This allows me to
% use different flat field images for multiple movies acquired in one day,
% which is necessary because they might be imaged at different angles
% (something that my microscope supports).
    
%PrefixRow=find(strcmp(XLSTxt(:,DataFolderColumn),[Prefix(1:10),'\',Prefix(12:end)])|...
%    strcmp(XLSTxt(:,DataFolderColumn),[Prefix(1:10),'/',Prefix(12:end)]));

ExperimentType=XLSTxt(PrefixRow,ExperimentTypeColumn);
Channel1=XLSTxt(PrefixRow,Channel1Column);
Channel2=XLSTxt(PrefixRow,Channel2Column);

%Now that we know all of this, reload the folders taking into account the
%new Dropbox one
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,SchnitzcellsFolder]=...
    DetermineLocalFolders(Prefix);

%Create the output folder
OutputFolder=[FISHPath,filesep,'Data',filesep,...
    Prefix];
mkdir(OutputFolder)


%Get the structure with the acquisition information
ImageInfo = imfinfo([Folder,filesep,D(1).name]);

%Do we have a second channel for Histone?
if strcmp(Channel2,'His-RFP')
    HisChannel=1;
else
    HisChannel=0;
end


if strcmp(FileMode,'TIF')

    %Get the flat-field information
    %Figure out the zoom factor
    Zoom=ExtractInformationField(ImageInfo(1),'state.acq.zoomFactor=');
    %Look for the file
    FFDir=dir([Folder,filesep,'..',filesep,'*FF',Zoom(1:end-1),'x*.*']);



    %If there's more than one match then ask for help
    if length(FFDir)==1
        FFFile=FFDir(1).name;
    elseif isempty(FFDir)
        display('Warning, no flat field file found. Press any key to proceed without it');
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
        NImagesForTAG=NImages;
    else
        %If there is a second channel for histone we will still use the
        %alignment information from the first channel.
        jChannel1=1:NRepeats*2:numel(ImageInfo);
        kChannel1=1:2:NRepeats*2;
        lChannel1=2:2:NRepeats*2;
        NImagesForTAG=NImages/2;
    end



    %If the images were not averaged we will want to align them and average
    %them ourselves.

    %Check that we didn't just want to generate the TAG file
    if ~TAGOnly
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
    end




    if AverageFlag
        h=waitbar(0,'Copying to FISH folder');
    else
        h=waitbar(0,'Aligning images and copying to FISH folder');
    end

    for i=1:length(D)
        Suffix{i}=[iIndex(i,3),'_z??'];
        ImageInfo = imfinfo([Folder,filesep,D(i).name]);
        waitbar(i/length(D),h)
        FrameInfo(i)=ExtractImageInformation(ImageInfo(1));

        %Check that we don't just want to calculate the TAG file
        if ~TAGOnly
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
    end
    
    %Add the information about the mode
    for i=1:length(FrameInfo)
        FrameInfo(i).FileMode='TIF';
    end
    
    
    close(h)
    
    %TAG file information
    Output{1}=['id ',Prefix,'_'];
    Output{2}='';
    Output{3}='1';
    Output{4}=['frames ',num2str(length(FrameInfo)),':1:',num2str(NImagesForTAG+2)];
    Output{5}=['suffix ???_z??'];
    Output{6}=['flat FF'];
    
%LSM mode
elseif strcmp(FileMode,'LSM')
    
    warning('Still need to add the FF information')
    
    %Parse all the time information from the files
    for i=1:length(D)
        parselsmwithtime([Folder,filesep,D(i).name])
    end
    
    %Load the information
    for i=1:length(D)
        LSMData(i)=load([Folder,filesep,D(i).name(1:end-3),'mat']);
    end
    
    %Generate FrameInfo
    FrameInfo=struct('LinesPerFrame',{},'PixelsPerLine',{},...
        'NumberSlices',{},'ZStep',{},'FileMode',{},...
        'PixelSize',{});
    %Some that I'm not assigning: ZoomFactor, Rotation
    for i=1:length(D)
        for j=1:LSMData(i).Datas.LSM_info.DimensionTime
            FrameInfo(end+1).LinesPerFrame=LSMData(i).Datas.LSM_info.DimensionY;
            FrameInfo(end).PixelsPerLine=LSMData(i).Datas.LSM_info.DimensionX;
            FrameInfo(end).NumberSlices=LSMData(i).Datas.LSM_info.DimensionZ;
            FrameInfo(end).FileMode='LSM';
            FrameInfo(end).PixelSize=LSMData(i).Datas.LSM_info.VoxelSizeX;
            FrameInfo(end).ZStep=LSMData(i).Datas.LSM_info.VoxelSizeZ;
            FrameInfo(end).Time=eval(['LSMData(i).Datas.LSM_info.Timeinfo.Posi',num2str(j)]);
        end
    end
    
    %Copy the data
    m=1;
    NSlices=FrameInfo(m).NumberSlices;
    h=waitbar(0,'Extracting LSM images');
    for i=1:length(D)
        %Get this frame
        for j=1:LSMData(i).Datas.LSM_info.DimensionTime
            waitbar(m/length(FrameInfo),h)
            Im = LoadLsmToMat([Folder,filesep,D(i).name(1:end-4)],j);
            
            %Go through the slices
            HisRFP=zeros(size(Im.Slice1,1),size(Im.Slice1,2),NSlices);
            %Save a black image before and after the stack
            imwrite(uint16(zeros(size(Im.Slice1,1),size(Im.Slice1,2))),...
                    [OutputFolder,filesep,Prefix,'_',iIndex(m,3),'_z',...
                    iIndex(1,2),'.tif']);
            imwrite(uint16(zeros(size(Im.Slice1,1),size(Im.Slice1,2))),...
                [OutputFolder,filesep,Prefix,'_',iIndex(m,3),'_z',...
                iIndex(NSlices+2,2),'.tif']);
            
            for k=1:NSlices
                %Save the MCP-GFP information
                imwrite(eval(['Im.Slice',num2str(k),'(:,:,1)']),...
                    [OutputFolder,filesep,Prefix,'_',iIndex(m,3),'_z',...
                    iIndex(k+1,2),'.tif'])
                
                %Get the Histone-RFP image
                HisRFP(:,:,k)=eval(['Im.Slice',num2str(k),'(:,:,2)']);
            end
            %Project and save the Histone-RFP image
            HisRFP=max(HisRFP,[],3);
            imwrite(uint16(HisRFP),...
                [OutputFolder,filesep,Prefix,'-His_',iIndex(m,3),'.tif'])
            m=m+1;
        end
    end
    close(h)
    
   %TAG file information
    Output{1}=['id ',Prefix,'_'];
    Output{2}='';
    Output{3}='1';
    Output{4}=['frames ',num2str(length(FrameInfo)),':1:',num2str(NSlices+2)];
    Output{5}=['suffix ???_z??'];
    %Output{6}=['flat FF'];
    
%LIFExport mode
elseif strcmp(FileMode,'LIFExport')
    
    warning('Still need to add the FF information')
    warning('All information about the images is being entered manually')
    

    %Determine the number of slices
    NSlices=length(dir([Folder,filesep,D(1).name(1:end-11),'*ch00.tif']));
    
    SeriesFiles = dir([Folder,filesep,'Series*Properties.xml']);
    NSeries = length(SeriesFiles);
    % Number of LUT images
    N_LUT = length(dir([Folder,filesep,'Series*ch*LUT.tif']));
    
    %Determine the number of frames (subtract LUT), divide by both channels)
    NFrames=(length(D)-N_LUT)/(2*NSlices);
    
    %Generate FrameInfo
    FrameInfo=struct('LinesPerFrame',{},'PixelsPerLine',{},...
        'NumberSlices',{},'ZStep',{},'FileMode',{},...
        'PixelSize',{});
    %Some that I'm not assigning: ZoomFactor, Rotation
    
    %Extract time information from xml files
    Frame_Times = zeros(1,NFrames*NSlices);
    for i = 1:length(SeriesFiles)
        xDoc = xmlread([Folder,filesep,SeriesFiles(i).name]);
        TimeStampList = xDoc.getElementsByTagName('TimeStamp');
        if i == 1
            start_index = 0;
        else
            start_index = ceil((nnz(Frame_Times)+1)/42)*42+1;
        end
        for k = 0:TimeStampList.getLength-1
            TimeStamp = TimeStampList.item(k);
            Date = char(TimeStamp.getAttribute('Date'));
            Time = char(TimeStamp.getAttribute('Time'));
            Milli = char(TimeStamp.getAttribute('MiliSeconds'));
            time_in_days = datenum(strcat(Date,'-',Time,'-',Milli),'dd/mm/yyyy-HH:MM:SS AM-FFF');
            Frame_Times(start_index+k+1)=time_in_days*86400;
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
    
    for i=1:NFrames
        FrameInfo(i).LinesPerFrame=nan;
        FrameInfo(i).PixelsPerLine=nan;
        FrameInfo(i).NumberSlices=NSlices;
        FrameInfo(i).FileMode='LIFExport';
        FrameInfo(i).PixelSize=0.22;
        FrameInfo(i).ZStep=0.5;
        % EDIT TIMESTAMP FROM XML RATHER THAN GUESSING 40
        FrameInfo(i).Time=Frame_Times(1+NSlices*i);
    end
    
    %Copy the data
    slice_num=1;
    % absolute frame
    m=1;
    % frame without blanks
    n=1;
    %records which frames were deleted for being blank
    blank_frameindex = zeros(1,NFrames);
    blank_frames = zeros(1,NFrames);
    h=waitbar(0,'Extracting LIFExport images');
    
    
    OriginalFrameCounter=1;
    for i = 1:length(D) 
        if ~isempty(findstr(D(i).name,'_ch00'))
            waitbar(i/length(D),h)
            Image=uint16(double(imread([Folder,filesep,D(i).name])));
            if max(Image(:))==0  
                blank_frames(m) = blank_frames(m)+1;
                blank_frameindex(m) = m;
                if slice_num == NSlices
                    m=m+1;
                    slice_num = 1;
                else
                    slice_num = slice_num+1;
                end
            elseif (max(Image(:))~=0)&&(blank_frames(m)==0)
                NewName=[Prefix,'_',iIndex(n,3),'_z',iIndex(slice_num,2),'.tif'];
                imwrite(Image,[OutputFolder,filesep,NewName]);
                if slice_num == NSlices
                    m=m+1;
                    n=n+1;
                    slice_num = 1;
                else
                    slice_num = slice_num+1;
                end
            end
        end
    end
    close(h)
    blank_frameindex(blank_frameindex==0) = [];
    
    %Remove FrameInfo information for blank frames
    FrameInfo(blank_frameindex)=[];
    
    %Do the histone channel, channel 01
    slice_num=1;
    % absolute frame
    m=1;
    % frame without blanks
    n=1;
    h=waitbar(0,'Histone images');
    for i = 1:length(D) 
        if ~isempty(findstr(D(i).name,'_ch01'))
            waitbar(i/length(D),h)
            Hist(:,:,slice_num)=uint16(double(imread([Folder,filesep,D(i).name])));
            if blank_frames(m)==0
                if slice_num == NSlices
                    HisRFP=max(Hist,[],3);
                    imwrite(uint16(HisRFP),...
                    [OutputFolder,filesep,Prefix,'-His_',iIndex(n,3),'.tif'])
                    m=m+1;
                    n=n+1;
                    slice_num = 1;
                else
                    slice_num = slice_num+1;
                end
            else
                if slice_num == NSlices
                    m=m+1;
                    slice_num = 1;
                else
                    slice_num = slice_num+1;
                end
            end
        end
    end
    close(h)
    
    % TAG File Information
    Output{1}=['id ',Prefix,'_'];
    Output{2}='';
    Output{3}='1';
    Output{4}=['frames ',num2str(length(FrameInfo)),':1:',num2str(NSlices)];
    Output{5}=['suffix ???_z??'];
    %Output{6}=['flat FF'];
    
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
    
    %Rename all MCP channel images
    %Find all the MCP files
    D=dir([OutputFolder,filesep,'*_z01.tif']);
    %Delete the skipped files
    for i=SkipFrames
        for j=1:NSlices
            delete([OutputFolder,filesep,D(i).name(1:end-6),iIndex(j,2),'.tif'])
        end
    end
    %Rename all remaining files
    D=dir([OutputFolder,filesep,'*_z01.tif']);
    for i=1:length(D)
        if ~strcmp([OutputFolder,filesep,D(i).name],...
                    [OutputFolder,filesep,D(i).name(1:end-11),iIndex(i,3),'_z01.tif'])
            for j=1:NSlices
                movefile([OutputFolder,filesep,D(i).name(1:end-6),iIndex(j,2),'.tif'],...
                    [OutputFolder,filesep,D(i).name(1:end-11),iIndex(i,3),'_z',iIndex(j,2),'.tif'])
            end
        end
    end
    
    %Redefine the output for the FISH code
    Output{4}=['frames ',num2str(length(FrameInfo)),':1:',num2str(NSlices)];
end



%Create the TAG file for the FISH analysis code
fid = fopen([OutputFolder,filesep,Prefix,'.tag'], 'wt');

for i=1:length(Output)
    fprintf(fid, '%s \n', Output{i});
end


fclose(fid);


%Save the information about the various frames
mkdir([DropboxFolder,filesep,Prefix])
save([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'],...
    'FrameInfo')

