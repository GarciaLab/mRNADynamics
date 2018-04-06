function FrameInfo = processZeissConfocalLSMData(Folder, D, ExperimentType, Channel1, Channel2, Prefix, OutputFolder)
    %warning('Still need to add the FF information') NL: I think this
    %warning is out-dated
    
    %What type of experiment do we have?
    if strcmp(ExperimentType,'1spot') || strcmp(ExperimentType,'2spot') || strcmp(ExperimentType,'2spot1color')
    
        %Figure out the different channels
        if ~isempty(strfind(Channel1{1},'MCP'))
            coatChannel=1;
        elseif  strfind(Channel1{1},'His')
            histoneChannel=1;
        else
            error('LSM Mode error: Channel name not recognized. Check MovieDatabase')
        end

        if ~isempty(strfind(Channel2{1},'MCP'))
            coatChannel=2;
        elseif  strfind(Channel2{1},'His')
            histoneChannel=2;
        else
            error('LSM Mode error: Channel name not recognized. Check MovieDatabase')
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
            LSMDir=dir([Folder,filesep,'*.lsm']);     %Zeiss confocal, old LSM format
            CZIDir=dir([Folder,filesep,'*.czi']);     %Zeiss confocal, new CZI format
            if ~isempty(LSMDir)
                StartingTime(LSMIndex) = LSMMeta2.get(['TimeStamp #',iIndex(1,NDigits)]);
            elseif ~isempty(CZIDir)
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
end