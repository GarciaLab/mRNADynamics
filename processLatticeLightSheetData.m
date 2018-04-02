 function FrameInfo = processLatticeLightSheetData(Folder, D, Channel1, Channel2, ProjectionType, Prefix, OutputFolder)    
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
        error('LAT Mode error: Channel name not recognized. Check MovieDatabase')
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
    
    %Figure out the different channels
    %MT, 2018-02-11: Yes, this is redundant, but it's a temporary thing -
    %FIX LATER
    coatChannel = 0;
    histoneChannel = 0;
    if ~isempty(strfind(Channel1{1},'MCP'))
        coatChannel=1;
    elseif  strfind(Channel1{1},'His')
        histoneChannel=1;
    else
        error('LAT Mode error: Channel name not recognized. Check MovieDatabase.XLSX')
    end

    if ~isempty(strfind(Channel2{1},'MCP'))
        coatChannel=2;
    elseif  strfind(Channel2{1},'His')
        histoneChannel=2;
    else
        error('LAT Mode error: Channel name not recognized. Check MovieDatabase.XLSX')
    end
    
    %Load the data
    im_stack = {};
    his_stack = [];
    mcp_stack = {};
    his_channel = 0;        %Keep track of whether data is his or mcp
    mcp_channel = 0;
    warning('While analyzing Lattice data, assuming Channel1 field in MovieDatabase.xlsx references CamA, and Channel2 references CamB.')
    h=waitbar(0,'Separating CamA/CamB into His and MCP Channels');
    for j = 1:length(D)
        waitbar(j/length(D),h)
        fname = [Folder, filesep, D(j).name];
        if contains(fname, 'CamA')
            if histoneChannel == 1
                his_channel = 1;
                mcp_channel = 0;
            else
                his_channel = 0;
                mcp_channel = 1;
            end
        elseif contains(fname,'CamB')
            if histoneChannel == 2
                his_channel = 1;
                mcp_channel = 0;
            else
                his_channel = 0;
                mcp_channel = 1;
            end
        end   
        info = imfinfo(fname);
        num_images = numel(info);
        for i = 1:num_images
            im_stack{j, i} = imread(fname, i, 'Info', info);
            if his_channel && ~mcp_channel
                his_stack{j,i} = imread(fname, i, 'Info', info);
                his_array(j,i, :, :) = imread(fname, i, 'Info', info);
            elseif ~his_channel && mcp_channel
                mcp_stack{j-size(his_stack, 1),i} = imread(fname, i, 'Info', info);
            else
                error('Something is wrong with your channels. Please doublecheck moviedatabase')
            end
        end
    end
    close(h)
    
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
    LATDir=dir([Folder,filesep,'*_Settings.txt']);
    metaID = fopen([Folder, filesep, LATDir.name]);
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
        NameSuffix=['_ch',iIndex(coatChannel,2)];
        
        NewName=[Prefix,'_',iIndex(j,3),'_z',iIndex(1,2),NameSuffix,'.tif'];
        imwrite(BlankImage,[OutputFolder,filesep,NewName]);
        NewName=[Prefix,'_',iIndex(j,3),'_z',iIndex(NSlices+2,2),NameSuffix,'.tif'];
        imwrite(BlankImage,[OutputFolder,filesep,NewName]);
        %Copy the rest of the images
        z = 2;
        for s = 1:NSlices
            NewName=[Prefix,'_',iIndex(j,3),'_z',iIndex(z,2),NameSuffix,'.tif'];
            imwrite(mcp_stack{j,s},[OutputFolder,filesep,NewName]);
            z = z + 1;
        end
        m=m+1;
    end
    close(h)
end
