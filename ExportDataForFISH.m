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
 
        LIFDir=dir([Folder,filesep,'*.lif']);     %Leica confocal
        %Load the file using BioFormats
        %Figure out which one is not the FF
        LIFIndex=find(cellfun(@isempty,strfind({LIFDir.name},'FF')));
        %Load the data, this might cause problems with really large sets
        LIFImages=bfopen([Folder,filesep,LIFDir(LIFIndex).name]);
        %Extract the metadata for each series
        LIFMeta = LIFImages{:, 4};
        % NSeries=LIFMeta.getImageCount(); %AR 2/4/2018 Not sure why this subtracts one, but it causes an error when there's only one series.
        % if NSeries == 0
        % NSeries = LIFMeta.getImageCount();
        % end
        NSeries = LIFMeta.getImageCount();
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
        % Frame_Times = [];
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
        
        if strcmpi(ExperimentType,'1spot') || strcmpi(ExperimentType,'2spot') ||...
                strcmpi(ExperimentType,'2spot1color') || strcmpi(ExperimentType,'inputoutput')
            %Figure out the different channels
            if ~isempty(Channel3)
                Channels={Channel1{1},Channel2{1},Channel3{1}};
            else
                Channels={Channel1{1},Channel2{1}};
            end
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
            % MCP-mCherry as a fake histone channel in case there's no
            % His-iRFP (Last edited : 3/28/2018, YJK)
            if isempty(fiducialChannel)&&...
                    ((~isempty(strfind(Channel1{1},'mCherry')))||(~isempty(strfind(Channel2{1},'mCherry'))))
                if (~isempty(strfind(Channel1{1},'mCherry')))
                    fiducialChannel=1;
                    histoneChannel=1;
                elseif (~isempty(strfind(Channel2{1},'mCherry')))
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
                elseif (~isempty(strfind(Channel2{1},'mCherry')))
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
    %       Comment out 'inputoutput' mode since it's now same as '1spot' mode
    %       (By YJK on 3/28/2018)
    %         elseif strcmpi(ExperimentType, 'inputoutput')
    %             if (~isempty(strfind(lower(Channel2),'mcp')))&...
    %                     ~isempty(strfind(lower(Channel2),'pcp'))
    %                 coatChannel=2;
    %             elseif (~isempty(strfind(lower(Channel1),'mcp')))&...
    %                     ~isempty(strfind(lower(Channel1),'pcp'))
    %                 coatChannel=1;
    %             else
    %                 error('No MCP or PCP channel detected. Check MovieDatabase')
    %             end
    %             if (~isempty(strfind(Channel1{1},'mCherry')))|(~isempty(strfind(Channel2{1},'mCherry')))
    %                 if (~isempty(strfind(Channel1{1},'mCherry')))
    %                     fiducialChannel=1;
    %                     histoneChannel=1;
    %                 elseif (~isempty(strfind(Channel2{1},'mCherry')))
    %                     fiducialChannel=2;
    %                     histoneChannel=2;
    %                 else
    %                     error('mCherry channel not found. Cannot generate the fake nuclear image')
    %                 end
    %             end
        else
            error('Experiment type not recognized. Check MovieDatabase')
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
                        
                        %YJK : Think about the case when there is no His channel,
                        %and it is inputoutput mode or 1spot mode or 2spot2color.
                        %We can use (MCP-mCherry) either inverted or raw
                        %images to make fake histone images.
                        if (isempty(strfind(Channel1{1}, 'His')))&&(isempty(strfind(Channel2{1}, 'His')))&&(isempty(Channel3))
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
