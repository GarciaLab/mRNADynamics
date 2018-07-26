function FrameInfo=ExtractImageInformation(ImageInfo, LatticeTxt)

%This function returns a structure with the relevant information related to
%the acquisition parameters

%Imaging conditions
FrameInfo.LinesPerFrame=ImageInfo.Height;
FrameInfo.PixelsPerLine=ImageInfo.Width;
FrameInfo.ZoomFactor=NaN;
FrameInfo.Rotation=NaN;

FrameInfo.PixelSize=0.104; % µm

%Time
FrameInfo.TimeString=num2str(Lattice_FindFromFilename(ImageInfo.Filename, '', 'msec')/1000); % in seconds

% %Position
% FrameInfo.XPosition=str2num(ExtractInformationField(ImageInfo,'state.motor.absXPosition='));
% FrameInfo.YPosition=str2num(ExtractInformationField(ImageInfo,'state.motor.absYPosition='));
% FrameInfo.ZPosition=str2num(ExtractInformationField(ImageInfo,'state.motor.absZPosition='));

%Stack
settingstxt = fileread(LatticeTxt);

findstr = 'Z PZT Offset, Interval (um), # of Pixels for Excitation (1)';
pos = strfind(settingstxt, findstr);
cutstr = settingstxt(pos-4:pos-1);
FrameInfo.NumberSlices=str2num(cutstr);

findstr = 'S PZT Offset, Interval (um), # of Pixels for Excitation (1)';
pos = strfind(settingstxt, findstr);
cutstr = settingstxt(pos-8:pos-5);
FrameInfo.ZStep=str2num(cutstr);
disp 'Warning: Z information retrieved really badly from the settings file. improve.'




