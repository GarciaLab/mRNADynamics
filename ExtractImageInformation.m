function FrameInfo=ExtractImageInformation(ImageInfo)

%This function returns a structure with the relevant information related to
%the acquisition parameters

%Imaging conditions
FrameInfo.LinesPerFrame=str2num(ExtractInformationField(ImageInfo,'state.acq.linesPerFrame='));
FrameInfo.PixelsPerLine=str2num(ExtractInformationField(ImageInfo,'state.acq.pixelsPerLine='));
FrameInfo.ZoomFactor=str2num(ExtractInformationField(ImageInfo,'state.acq.zoomFactor='));
FrameInfo.Rotation=str2num(ExtractInformationField(ImageInfo,'state.acq.scanRotation='));
FrameInfo.ScanAmplitudeX=str2num(ExtractInformationField(ImageInfo,'state.acq.scanAmplitudeX='));
FrameInfo.ScanAmplitudeY=str2num(ExtractInformationField(ImageInfo,'state.acq.scanAmplitudeY='));



%Time
FrameInfo.TimeString=ExtractInformationField(ImageInfo,'state.internal.triggerTimeString=');
FrameInfo.TimeString=FrameInfo.TimeString(2:end-2);

%Position
FrameInfo.XPosition=str2num(ExtractInformationField(ImageInfo,'state.motor.absXPosition='));
FrameInfo.YPosition=str2num(ExtractInformationField(ImageInfo,'state.motor.absYPosition='));
FrameInfo.ZPosition=str2num(ExtractInformationField(ImageInfo,'state.motor.absZPosition='));

%Stack
FrameInfo.NumberSlices=str2num(ExtractInformationField(ImageInfo,'state.acq.numberOfZSlices='));
FrameInfo.ZStep=str2num(ExtractInformationField(ImageInfo,'state.acq.zStepSize='));



