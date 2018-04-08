function FrameInfo = process2PhotonPrincetonData(Folder, D, FrameInfo, Channel2, MaxShift, MaxHistone, OutputFolder)    
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
end