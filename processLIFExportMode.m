function FrameInfo = processLIFExportMode(Folder, ExperimentType, ProjectionType, Channel1, Channel2, Channel3, Prefix, OutputFolder)
  %Extract time information from xml files
  XMLFolder = Folder;
  SeriesFiles = dir([XMLFolder, filesep, '*Series*Properties.xml']);
  if isempty(SeriesFiles)
      XMLFolder = [Folder, filesep, 'MetaData'];
      SeriesFiles = dir([XMLFolder, filesep, '*Series*Properties.xml']);
      if isempty(SeriesFiles)
          error('XML MetaFiles could not be found. Did they get exported using the LAS software?')
      end
  end

  %Leica confocal
  LIFDir = dir([Folder,filesep,'*.lif']);

  %Load the file using BioFormats
  %Figure out which one is not the FF
  LIFIndex = find(cellfun(@isempty, strfind({LIFDir.name}, 'FF')));

  %Load the data, this might cause problems with really large sets
  LIFImages = bfopen([Folder, filesep, LIFDir(LIFIndex).name]);

  %Extract the metadata for each series
  LIFMeta = LIFImages{:, 4};

  % NSeries=LIFMeta.getImageCount(); %AR 2/4/2018 Not sure why this subtracts one, but it causes an error when there's only one series.
  % if NSeries == 0
  % NSeries = LIFMeta.getImageCount();
  % end
  NSeries = LIFMeta.getImageCount();
  
  %Figure out the number of slices in each series
  NSlices = [];
  for i = 1:NSeries
      NSlices(i) = str2num(LIFMeta.getPixelsSizeZ(i-1));
  end

  %Number of planes per series
  NPlanes = [];
  for i = 1:NSeries
      NPlanes(i) = LIFMeta.getPlaneCount(i-1);
  end

  %Number of channels
  NChannels = LIFMeta.getChannelCount(0);

  %Finally, use this information to determine the number of frames in each series        
  NFrames = NPlanes./NSlices/NChannels;

  %Get rid of the last frame as it is always incomplete because that's when we stopped it
  NFrames = NFrames - 1;
  NPlanes = NPlanes - NSlices * NChannels;      
  Frame_Times = zeros(1, sum(NFrames.*NSlices));
  % Frame_Times = [];

  m = 1;
  for i = 1:NSeries
    xDoc = xmlread([XMLFolder,filesep,SeriesFiles(i).name]);
    TimeStampList = xDoc.getElementsByTagName('TimeStamp');
    for k = 0:(NFrames(i) * NSlices(i) * NChannels) - 1
        TimeStamp = TimeStampList.item(k);
        Date = char(TimeStamp.getAttribute('Date'));
        Time = char(TimeStamp.getAttribute('Time'));
        Milli = char(TimeStamp.getAttribute('MiliSeconds'));
        time_in_days = datenum(strcat(Date, '-', Time, '-', Milli), 'dd/mm/yyyy-HH:MM:SS AM-FFF');
        Frame_Times(m) = time_in_days*86400;
        m=m+1;
    end
  end
  
  First_Time=Frame_Times(1);
  for i = 1:length(Frame_Times)
    if Frame_Times(i) == 0
      Frame_Times(i) = 0;
    else
      Frame_Times(i) = Frame_Times(i) - First_Time;
    end
  end

  %Get the time stamp corresponding to the first slice of each Z-stack
  m = 1;
  for i = 1:NSeries
    if i == 1
      StartIndex = 1;
    else
      StartIndex = sum(NPlanes(1:i-1)) + 1;
    end
    for j = StartIndex:(NSlices(i).*NChannels):sum(NPlanes(1:i))
      InitialStackTime(m) = Frame_Times(j);
      m = m + 1;
    end
  end

  for i = 1:sum(NFrames)
    FrameInfo(i).LinesPerFrame = str2double(LIFMeta.getPixelsSizeY(0));
    FrameInfo(i).PixelsPerLine = str2double(LIFMeta.getPixelsSizeX(0));
    FrameInfo(i).NumberSlices = min(NSlices);
    FrameInfo(i).FileMode = 'LIFExport';
    %This is to allow for backwards compatibility with BioFormats
    if ~isempty(str2num(LIFMeta.getPixelsPhysicalSizeX(0)))
        FrameInfo(i).PixelSize = str2num(LIFMeta.getPixelsPhysicalSizeX(0));
        FrameInfo(i).ZStep = str2double(LIFMeta.getPixelsPhysicalSizeZ(0));
    else
        FrameInfo(i).PixelSize = str2num(LIFMeta.getPixelsPhysicalSizeX(0).value);
        FrameInfo(i).ZStep = str2double(LIFMeta.getPixelsPhysicalSizeZ(0).value);
    end
    FrameInfo(i).Time = InitialStackTime(i);
  end
  %Find the flat field (FF) information
  
  %The flat field image can be in the folder with the data or in the folder
  %corresponding to the date.
  D1 = dir([Folder, filesep, 'FF*.lif']);
  D2 = dir([Folder, filesep, '..', filesep, 'FF*.lif']);
  FFPaths = {};
  for i = 1:length(D1)
    FFPaths{end+1}=[Folder, filesep, D1(i).name];
  end
  for i = 1:length(D2)
    FFPaths{end+1} = [Folder, filesep, '..', filesep,D2(i).name];
  end
  
  %Go through the FF files and see which one matches the pixel size
  %and image pixel number
  FFToUse = [];
  for i = 1:length(FFPaths)
    LIFFF = bfopen(FFPaths{i});
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
                FFToUse = [FFToUse, i];
            end
        end
    catch
      if (str2num(LIFFFMeta.getPixelsPhysicalSizeX(0))==...
              str2num(LIFMeta.getPixelsPhysicalSizeX(0)))&...
              (str2num(LIFFFMeta.getPixelsSizeY(0))==...
              str2num(LIFMeta.getPixelsSizeY(0)))&...
              (str2num(LIFFFMeta.getPixelsSizeX(0))==...
              str2num(LIFMeta.getPixelsSizeX(0)))
          FFToUse = [FFToUse, i];
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
              FFToUse = [FFToUse, i];
          end
      end
    end
  end 
  
  if length(FFToUse)> 1
    warning('Too many flat field images match the pixel and image size size')
    [FFFile,FFPath] =...
        uigetfile([Folder, filesep, '*.lif'], 'Select which flat field image to use');
    LIFFF = bfopen([FFPath, FFFile]);
  elseif isempty(FFToUse)
    warning('No flat field image found')
    clear LIFFF
  else
    LIFFF = bfopen(FFPaths{FFToUse});
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
    if (fiducialChannel==0)&&...
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
    Channels = {Channel1{1}, Channel2{1}};
    
    %We have no coat protein here.
    coatChannel = 0;
    
    %Histone channel.
    histoneChannel = find(~cellfun(@isempty,strfind(lower(Channels),'his')));
    if isempty(histoneChannel)
      histoneChannel = 0;
    else
      fiducialChannel = histoneChannel;
    end

    %Input channels
    inputProteinChannel=~cellfun(@isempty,Channels);
    if histoneChannel
        inputProteinChannel(histoneChannel) = 0;
    else
      %If there was no histone channel, we need to choose which
      %input channel to use as our fiducial channel. We'll use
      %the first channel for now. We can try to be smarted about
      %this later on.
      warning('No histone channel found. Finding nuclei using the protein input channel.')
      fiducialChannel = 1;                
    end
    inputProteinChannel = find(inputProteinChannel);
    
    %Save the information about the number of channels in FrameInfo
    for i = 1:length(FrameInfo)
      FrameInfo(i).NChInput = length(inputProteinChannel);
    end
  else
    error('Experiment type not recognized. Check MovieDatabase')
  end

  %Copy the data
  h = waitbar(0, 'Extracting LIFExport images');
  %Create a blank image
  BlankImage = uint16(zeros(size(LIFImages{1}{1,1})));
  m = 1;        %Counter for number of frames
  %Load the reference histogram for the fake histone channel
  load('ReferenceHist.mat')
  for i = 1:NSeries
    waitbar(i/NSeries,h)
    for j=1:NFrames(i) 
        for q=1:NChannels
            if (strcmpi(ExperimentType,'1spot') ||...
                    strcmp(ExperimentType,'2spot') ||...
                    strcmp(ExperimentType,'2spot1color')) && ...
                    q == coatChannel
              %Save the blank images at the beginning and end of the stack
              NameSuffix = ['_ch',iIndex(q,2)];
              NewName = [Prefix, '_', iIndex(m,3), '_z', iIndex(1,2), NameSuffix, '.tif'];
              imwrite(BlankImage, [OutputFolder, filesep, NewName]);
              NewName = [Prefix, '_', iIndex(m,3), '_z', iIndex(min(NSlices)+2, 2), NameSuffix, '.tif'];
              imwrite(BlankImage, [OutputFolder, filesep, NewName]);
              %Copy the rest of the images
              n = 1;        %Counter for slices
              firstImage = (j-1) * NSlices(i) * NChannels + 1 + (q - 1);
              lastImage = j * NSlices(i) * NChannels;
              for k = firstImage:NChannels:lastImage
                if n <= min(NSlices)
                    NewName = [Prefix, '_', iIndex(m,3), '_z', iIndex(n + 1, 2), NameSuffix, '.tif'];
                       imwrite(LIFImages{i}{k,1}, [OutputFolder, filesep, NewName]);
                    n = n + 1;
                end
              end
            elseif strcmpi(ExperimentType,'2spot2color')
              NameSuffix=['_ch',iIndex(q,2)];

              %Save the blank images at the beginning and end of the
              %stack
              NewName = [Prefix, '_', iIndex(m,3), '_z', iIndex(1,2), NameSuffix, '.tif'];
              imwrite(BlankImage, [OutputFolder, filesep, NewName]);
              NewName = [Prefix, '_', iIndex(m,3), '_z', iIndex(min(NSlices) + 2, 2), NameSuffix, '.tif'];
              imwrite(BlankImage, [OutputFolder, filesep, NewName]);
              %Copy the rest of the images
              n = 1;        %Counter for slices
              firstImage = (j-1) * NSlices(i) * NChannels + 1 + (q - 1);
              lastImage = j * NSlices(i) * NChannels;
              for k = firstImage:NChannels:lastImage
                if n <= min(NSlices)
                  NewName = [Prefix, '_', iIndex(m,3), '_z', iIndex(n + 1, 2), NameSuffix, '.tif'];
                     imwrite(LIFImages{i}{k,1}, [OutputFolder, filesep, NewName]);
                  n=n+1;
                end
              end
            %input-output mode
            elseif strcmpi(ExperimentType, 'inputoutput')
                %are we dealing with the coat channel?
                if q == coatChannel
                    %Save the blank images at the beginning and end of the stack
                    NameSuffix = ['_ch', iIndex(q, 2)];
                    NewName = [Prefix, '_', iIndex(m, 3), '_z', iIndex(1, 2), NameSuffix, '.tif'];
                    imwrite(BlankImage, [OutputFolder, filesep, NewName]);
                    NewName = [Prefix, '_', iIndex(m, 3), '_z', iIndex(min(NSlices) + 2, 2), NameSuffix, '.tif'];
                    imwrite(BlankImage, [OutputFolder, filesep, NewName]);
                    %Copy the rest of the images
                    n = 1;        %Counter for slices
                    firstImage = (j-1) * NSlices(i) * NChannels + 1 + (q - 1);
                    lastImage = j * NSlices(i) * NChannels;

                    TempNameSuffix = ['_ch', iIndex(q, 2)];
                    for k = firstImage:NChannels:lastImage
                      if n <= min(NSlices)
                          NewName = [Prefix, '_', iIndex(m, 3), '_z', iIndex(n + 1, 2), TempNameSuffix, '.tif'];
                             imwrite(LIFImages{i}{k,1}, [OutputFolder, filesep, NewName]);
                          n = n + 1;
                      end
                    end
                    
                %This is for the input channel    
                else
                  %Save the blank images at the beginning and end of the stack
                  NameSuffix = ['_ch', iIndex(q, 2)];
                  NewName = [Prefix, '_', iIndex(m, 3), '_z', iIndex(1, 2), NameSuffix, '.tif'];
                  imwrite(BlankImage, [OutputFolder, filesep, NewName]);
                  NewName = [Prefix, '_', iIndex(m, 3), '_z', iIndex(min(NSlices) + 2, 2), NameSuffix, '.tif'];
                  imwrite(BlankImage, [OutputFolder, filesep, NewName]);
                  %Copy the rest of the images
                  n = 1;        %Counter for slices
                  firstImage = (j-1) * NSlices(i) * NChannels + 1 + (q - 1);
                  lastImage = j * NSlices(i) * NChannels;

                  TempNameSuffix = ['_ch', iIndex(q, 2)];
                  for k = firstImage:NChannels:lastImage
                    if n <= min(NSlices)
                        NewName = [Prefix, '_', iIndex(m, 3), '_z', iIndex(n + 1, 2), TempNameSuffix, '.tif'];
                           imwrite(LIFImages{i}{k,1}, [OutputFolder, filesep, NewName]);
                        n = n + 1;
                    end
                  end
                end
                   
            elseif strcmpi(ExperimentType, 'input') && sum(q == inputProteinChannel)
                %Are we dealing with one or two channels?
                if length(inputProteinChannel) == 1
                  NameSuffix = ['_ch', iIndex(q, 2)];
                else
                  NameSuffix=['_ch', iIndex(q, 2)];
                end
                %Save the blank images at the beginning and end of the stack
                NewName = [Prefix, '_', iIndex(m, 3), '_z', iIndex(1, 2), NameSuffix, '.tif'];
                imwrite(BlankImage, [OutputFolder, filesep, NewName]);
                NewName = [Prefix, '_', iIndex(m,3), '_z', iIndex(min(NSlices) + 2, 2), NameSuffix, '.tif'];
                imwrite(BlankImage, [OutputFolder, filesep, NewName]);
                %Copy the rest of the images
                n = 1;        %Counter for slices
                firstImage = (j-1) * NSlices(i) * NChannels + 1 + (q - 1);
                lastImage = j * NSlices(i) * NChannels;
                for k = firstImage:NChannels:lastImage
                  if n <= min(NSlices)
                    NewName = [Prefix, '_', iIndex(m, 3), '_z', iIndex(n + 1, 2), NameSuffix, '.tif'];
                       imwrite(LIFImages{i}{k,1}, [OutputFolder, filesep, NewName]);
                    n=n+1;
                  end
                end             
            end
        end
        %Now copy nuclear tracking images
        if fiducialChannel
          HisSlices = zeros([size(LIFImages{i}{1,1},1), size(LIFImages{i}{1,1},2), NSlices(i)]);
          otherSlices = HisSlices;
          n = 1;
          firstImage = (j-1) * NSlices(i) * NChannels + 1 + (fiducialChannel - 1);
          lastImage = j * NSlices(i) * NChannels;
          for k = firstImage:NChannels:lastImage
              HisSlices(:,:,n) = LIFImages{i}{k,1};
              otherSlices(:,:,n) = LIFImages{i}{k+1,1};
              n = n + 1;
          end
          if histoneChannel
            if strcmp(ProjectionType,'medianprojection')
              Projection = median(HisSlices, 3);
            else
              Projection = max(HisSlices, [], 3);
            end
            
            %YJK : Think about the case when there is no His channel,
            %and it is inputoutput mode, 1spot mode or 2spot2color.
            %We can use (MCP-mCherry) either inverted or raw
            %images to make fake histone images.
            
            if (isempty(strfind(Channel1{1}, 'His')))&&(isempty(strfind(Channel2{1}, 'His')))&&(isempty(strfind(Channel3{1}, 'His')))
                if strcmpi(ExperimentType, 'inputoutput')|strcmpi(ExperimentType, '1spot')|strcmpi(ExperimentType,'2spot2color')|strcmpi(ExperimentType,'input')

                if (~isempty(strfind(Channel1{1}, 'NLS')))|(~isempty(strfind(Channel2{1}, 'NLS')))
                  %don't invert with NLS-MCP-mCherry
                else
                  %We don't want to use all slices. Only the center ones
                  StackCenter = round((min(NSlices) - 1) / 2);
                  StackRange = StackCenter - 1:StackCenter + 1;
                  if strcmp(ProjectionType, 'medianprojection')
                      Projection = median(HisSlices(:,:,StackRange), [], 3);
                  else
                      Projection = max(HisSlices(:,:,StackRange), [], 3);
                  end
                  %invert images to make nuclei bright
                  Projection = imcomplement(Projection);
                end
                Projection = histeq(mat2gray(Projection), ReferenceHist);
                Projection = Projection * 10000;
              end
            end
          else 
              %We don't want to use all slices. Only the center ones
              StackCenter = round((min(NSlices) - 1) / 2);
              StackRange = StackCenter - 1:StackCenter + 1;
              if strcmp(ProjectionType, 'medianprojection')
                Projection = median(HisSlices(:,:,StackRange), 3);
                otherProjection= median(otherSlices(:,:,StackRange), 3);
              else
                Projection = max(HisSlices(:,:,StackRange), [], 3);
                otherProjection = max(otherSlices(:,:,StackRange), [], 3);
              end
              Projection = Projection + otherProjection;
          end

          imwrite(uint16(Projection),...
          [OutputFolder, filesep, Prefix, '-His_', iIndex(m, 3), '.tif']);
        end
    m = m + 1;
    end
  end
  close(h)
end
