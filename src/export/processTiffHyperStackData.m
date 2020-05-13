function FrameInfo = processTiffHyperStackData(Folder, D, FrameInfo, Channel1, Channel2, OutputFolder, Prefix)

%Get the structure with the acquisition information
ImageInfo = imfinfo([Folder,filesep,D(1).name]);
% examine file names to determine channel IDS
files = {D.name};
if length(files) > 1
    error('Code does not currently support multiple channels')
end

channels = {char(Channel1) char(Channel2)};
% find mcrp channel 
MCPChannel = find(contains(channels,'MCP'));
HisChannel = find(contains(channels,'His-RFP'));

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

%Fill in info for as many FrameInfo fields as possible using
%"ImageDescriptionField"
ImageDescriptionString = ImageInfo(1).ImageDescription;

% number of frames
FrameStart = strfind(ImageDescriptionString,'frames=')+7;
FrameStop = strfind(ImageDescriptionString,'hyperstack')-1;
NFrames = str2num(ImageDescriptionString(FrameStart:FrameStop));

% number of Z slices
SliceStart = strfind(ImageDescriptionString,'slices=')+7;
SliceStop = strfind(ImageDescriptionString,'frames')-1;
NSlices = str2num(ImageDescriptionString(SliceStart:SliceStop));

% image dims 
LinesPerFrame = ImageInfo(1).Height;
PixelsPerLine = ImageInfo(1).Width;

% pixel dimensions (dummy vals for now)
PixelSize = 0.1;
ZStep = 0.5;

for f = 1:NFrames
    FrameInfo(f).Time = f; % keep in frame units for now
    FrameInfo(f).NumberSlices = NSlices;
    FrameInfo(f).PixelsPerLine = PixelsPerLine;
    FrameInfo(f).LinesPerFrame = LinesPerFrame;
    FrameInfo(f).PixelSize = PixelSize;
    FrameInfo(f).ZStep = ZStep;   
    FrameInfo(f).FileMode = 'HyperTiff';   
end    


% read in raw data
disp('Loading raw image data...')
MCPFilter = contains(files,['ch' iIndex(MCPChannel,2)]);
RawMCPData = bfopen([Folder,filesep,D(MCPFilter).name]);

% now export images slice by slice
% disp('Exporting image slices...')
h=waitbar(0,'Exporting individual Tiff slices');
for f = 1:NFrames
           
    waitbar(f/NFrames,h)
    
    %Use the size information of this image to calculate create a
    %blank image for the beginning and end of the stack.
    BlankImage=uint16(zeros(FrameInfo(f).LinesPerFrame,FrameInfo(f).PixelsPerLine));
    %Save the first image
    NewName=[Prefix,'_',iIndex(f,3),'_z',  iIndex(1,2),'_ch' iIndex(MCPChannel,2), '.tif'];
    imwrite(BlankImage,[OutputFolder,filesep,NewName],'Compression','none');
    %Save the last image
    NewName=[Prefix,'_',iIndex(f,3),'_z', iIndex(NSlices+2,2),'_ch' iIndex(MCPChannel,2), '.tif'];
    imwrite(BlankImage,[OutputFolder,filesep,NewName],'Compression','none');                
        
    parfor z=1:NSlices
        
        % extract image
        ImagePlane =  uint16(RawMCPData{1}{(f-1)*NSlices + z,1});

        % write to file
        NewName=[Prefix,'_',iIndex(f,3),'_z', iIndex(z+1,2),'_ch' iIndex(MCPChannel,2) '.tif'];
        imwrite(ImagePlane,[OutputFolder,filesep,NewName]);
      
    end

end

% %Get the actual time corresponding to each frame in seconds and add it to
% %FrameInfo
% for i=1:length(FrameInfo)
%     FrameInfo(i).Time=etime(datevec(FrameInfo(i).TimeString),datevec(FrameInfo(1).TimeString));
% end
% 


close(h)
end