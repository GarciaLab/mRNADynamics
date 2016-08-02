% Lattice file mode 
function [Output, FrameInfo] = ExportDataForFISH_Lattice(Prefix, D, Folder, OutputFolder, Channel1, Channel2, TAGOnly, ImageInfo)

    %Get the flat-field information
    %Figure out the zoom factor
%    disp 'TO DO: something about the zoom. Right now, we do not have whole-embryo images anyways.'
%    Zoom=ExtractInformationField(ImageInfo(1),'state.acq.zoomFactor=');
    %Look for the file
    FFDir=dir([Folder,filesep,'..',filesep,'*FF_ill*']);
    FFdarkDir=dir([Folder,filesep,'..',filesep,'*FF_dark*']);
    SmplDarkDir=dir([Folder,filesep,'..',filesep,'*sample_dark*']);
    %If there's more than one match then ask for help
    % load Flat Field
    if length(FFDir)==1
        FFFile=FFDir(1).name;
        FFImage = imread([Folder,filesep,'..',filesep,FFDir(1).name]);
    elseif isempty(FFDir)
        display('Warning, no flat field file found. Press any key to proceed without it');
        FFImage=ones(ImageInfo(1).Height,ImageInfo(1).Width);
        pause
    else
        FFFile=uigetfile([Folder,filesep,'..',filesep,'*.*'],'Select flatfield file');
        FFImage = imread([FFFile]);
    end
    % load dark count of Flat Field
    if length(FFdarkDir)==1
        FFdarkFile=FFdarkDir(1).name;
        FFdarkImage = imread([Folder,filesep,'..',filesep,FFdarkDir(1).name]);
    elseif isempty(FFdarkDir)
        display('Warning, no flat field dark image file found. Press any key to proceed without it');
        FFdarkImage=ones(ImageInfo(1).Height,ImageInfo(1).Width);
        pause
    else
        FFdarkFile=uigetfile([Folder,filesep,'..',filesep,'*.*'],'Select flatfield darkcount image file');
        FFdarkImage = imread([FFdarkFile]);
    end
    % load sample dark count
    if length(SmplDarkDir)==1
        SmplDarkFile=SmplDarkDir(1).name;
        SmplDarkImage = imread([Folder,filesep,'..',filesep,SmplDarkDir(1).name]);
    elseif isempty(SmplDarkDir)
        display('Warning, no sample dark image file found. Press any key to proceed without it');
        SmplDarkImage=ones(ImageInfo(1).Height,ImageInfo(1).Width);
        pause
    else
        SmplDarkFile=uigetfile([Folder,filesep,'..',filesep,'x*.*'],'Select sample darkcount image file');
        SmplDarkImage = imread([SmplDarkFile]);
    end
    % load dark count of sample images
    filtStd=30;         %This came from the FISH code.
%    FFImage=FFImage - uint16(mean(FFdarkImage(:))); % subtract dark counts.
    FFImage=FFImage - mean(FFdarkImage(:)); % subtract dark counts.
    FFImage=imfilter(FFImage,fspecial('gaussian',2*filtStd,filtStd),'symmetric');
    FFImage=imdivide(double(FFImage),double(max(FFImage(:))));
    % suppress data in regions where flatfield is less than flat field
    % threshold
    FFthreshold = .75;
    FFMask = zeros(size(FFImage));
    FFMask(FFImage<FFthreshold) = 1; % Mask is 1 where Image needs to be cut out.
    SmplDarkMean = uint16(mean(SmplDarkImage(:)));

    %dont copy into folder because of deskewing -> different resulting
    %image sizes. FF correction can only be made before deskewing
%     %If it's there, copy the image to the folder.
%     if ~isempty(FFDir)
%         FFImage=imread([Folder,filesep,'..',filesep,FFFile],1);
%         imwrite(FFImage,[OutputFolder,filesep,Prefix,'_FF.tif']);
%     end



    % find all the z-stacks corresponding to the two channels
    % TO DO: generalize to more channels
%     Channel0images = strcmp(D,'ch0');%dir([Folder,filesep,'*ch0*']);
%     Channel1images = strcmp(D,'ch1');%dir([Folder,filesep,'*ch1*']);
%     % check number of images and set them to the same number if different
%     if length(Channel0images) ~= length(Channel1images)
%         disp 'Warning: different number of time points for various channels. cutting off at the shortest.'
%         if length(Channel0images) < length(Channel1images)
%             Channel1images = Channel1images(1:length(Channel0images));
%         else
%             Channel0images = Channel0images(1:length(Channel1images));
%         end
%     end

    % HERE WAS the code to correlate and find best image shifts
    % TO DO: maybe insert it again to improve on deskewing. but not
    % necessary now


    % initialize a buffer image that contains one plane of the data
    BufferImageNoDeskew = imread([Folder,filesep,D(1).name],1);
    ImageInfo = imfinfo([Folder,filesep,D(1).name]);
    Cutoff = strfind(ImageInfo(1).Filename,'_ch');
    LatticeTxt = [ImageInfo(1).Filename(1:Cutoff), 'Settings.txt'];
    MyFrameInfo=Lattice_ExtractImageInformation(ImageInfo(1), LatticeTxt);
    % for deskewing later
    zStep = MyFrameInfo.ZStep;
    xStep = 0.104; % nm TODO: have the software generate a better metadata that includes this and read it
    angle = 31.8*pi/180; %angle between stage and bessel beam lattice
    shiftX = zStep/xStep *cos(angle);
    BufferImageDeskew = uint16(zeros(length(BufferImageNoDeskew(:,1)), length(BufferImageNoDeskew(1,:))+round(shiftX*length(ImageInfo))));
    HisImage = BufferImageDeskew;
    
    IsHis = 0;
    nrHis=0;

    h=waitbar(0,'Copying to FISH folder');

    for i=1:length(D)
        Suffix{i}=[iIndex(i,3),'_z??'];
        ImageInfo = imfinfo([Folder,filesep,D(i).name]);
        % get settings text file name and extract the image information
        Cutoff = strfind(ImageInfo(1).Filename,'_ch');
        LatticeTxt = [ImageInfo(1).Filename(1:Cutoff), 'Settings.txt'];
        MyFrameInfo=Lattice_ExtractImageInformation(ImageInfo(1), LatticeTxt);

        waitbar(i/length(D),h)

        % check if this is a histone channel or not for later processing
        if (~isempty(strfind(ImageInfo(1).Filename,'_ch0')) && strcmp(Channel1,'His-RFP')) ...
                || (~isempty(strfind(ImageInfo(1).Filename,'_ch1')) && strcmp(Channel2,'His-RFP'))
            IsHis = 1;
        else
            IsHis = 0;
        end
        %Check that we don't just want to calculate the TAG file
        if ~TAGOnly
            % load the image, do a flat field correction, then deskew it
            % and do a maximum projection or save it directly
            
            % deskew
            for j=1:length(ImageInfo)
                BufferImageNoDeskew = imread([ImageInfo(1).Filename],j)-SmplDarkMean; % To do: have one dark image per channel.
                % FLat Field correction in here!
                %BufferImageNoDeskew = imdivide(BufferImageNoDeskew,FFImage);
                BufferImageNoDeskew = uint16(imdivide(double(BufferImageNoDeskew),FFImage));
                %BufferImageNoDeskew = BufferImageNoDeskew-250;
                %disp 'bad background correction. improve.'
                % add background in regions that appeared by deskewing.
                % otherwise the DoG threshold needs to be set higher
                % doesn't need to be exact, only roughly
                meanbkg = sort([BufferImageNoDeskew(:,1)', BufferImageNoDeskew(:,end)', BufferImageNoDeskew(1,:), BufferImageNoDeskew(end,:)]);
                meanbkg = mean(meanbkg(1:floor(end/2)));
                BufferImageNoDeskew(FFMask>0)=meanbkg;
                % do the deskew
                BufferImageDeskew(:)=0;
                startx = round(shiftX*(j-1))+1;
                endx = startx+length(BufferImageNoDeskew(1,:))-1;
                xnotused = 1:length(BufferImageDeskew(1,:));
                xnotused(startx:endx)=[];
                BufferImageDeskew(:,startx:endx) = BufferImageNoDeskew;
                BufferImageDeskew(:,xnotused)=meanbkg;
                if IsHis % then only max-project and save after the loop
                    if j==1
                        HisImage(:)=0;
                    end
                    HisImage = max(BufferImageDeskew, HisImage);
                else % not His: save every frame.
                    NewName=[Prefix,'_',iIndex(Lattice_FindFromFilename(ImageInfo(1).Filename ,'stack', '')+1,3),'_z',...
                        iIndex(j+1,2),'.tif'];
                    imwrite(BufferImageDeskew,[OutputFolder,filesep,NewName]);
                    if j==1
                        % Add a blank image on top and bottom to get around border
                        % effects. Then loop through the stack and save the images.
                        BufferImageDeskew(:)=0;
                        NewName=[Prefix,'_',iIndex(Lattice_FindFromFilename(ImageInfo(1).Filename ,'stack', '')+1,3),'_z',...
                                    iIndex(1,2),'.tif'];
                        imwrite(BufferImageDeskew,[OutputFolder,filesep,NewName]);
                        NewName=[Prefix,'_',iIndex(Lattice_FindFromFilename(ImageInfo(1).Filename ,'stack', '')+1,3),'_z',...
                                    iIndex(length(ImageInfo)+2,2),'.tif'];
                        imwrite(BufferImageDeskew,[OutputFolder,filesep,NewName]);
                    end
                end
            end
            MyFrameInfo.ZStep = zStep * sin(angle);
            
            
            % if this is a histone stack, max-project it.
            if IsHis
                NewName=[Prefix,'-His_',iIndex(Lattice_FindFromFilename(ImageInfo(1).Filename, 'stack' ,'')+1,3),'.tif'];
                imwrite(HisImage,[OutputFolder,filesep,NewName]);
            end
        end
        if IsHis
            nrHis=nrHis+1;
            FrameInfo(nrHis) = MyFrameInfo; % this is just for compability with existing code
        end
    end
    
    %Add the information about the mode
    for i=1:length(FrameInfo)
        FrameInfo(i).FileMode='LAT';
    end
    
    NImagesForTAG = length(FrameInfo);
    
    close(h)
    
    % write spatial and temporal resolution into dataStructure and save in
    % the Analysis/Analsysis folder. probably, this is not the long-term
    % solution? I guess we are supposed to restrict ourselves to
    % Analysis/Raw Data here..
    disp('warning: saving dataStrucuture into Analysis folder from ExportDataForFISH_Lattice')
    % get time resolution
    if length(FrameInfo)>1
        deltaT = str2num(FrameInfo(2).TimeString)-str2num(FrameInfo(1).TimeString); % in s
    else
        deltaT = 0; % does not make sense with one frame
    end
    AcquisitionSettings.time_resolution = deltaT;
    AcquisitionSettings.space_resolution = FrameInfo(1).PixelSize;
%    mkdir([OutputFolder,filesep,'..',filesep,'..',filesep,'Analysis',filesep,Prefix,'_'])
    save([OutputFolder,filesep,Prefix,'-AcquisitionSettings.mat'], 'AcquisitionSettings')
    
    %TAG file information
    Output{1}=['id ',Prefix,'_'];
    Output{2}='';
    Output{3}='1';
    Output{4}=['frames ',num2str(length(FrameInfo)),':1:',num2str(length(ImageInfo)+2)];% #frames:1:#slices
    Output{5}=['suffix ???_z??'];
%    Output{6}=['flat FF'];