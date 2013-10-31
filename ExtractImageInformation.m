function FrameInfo=ExtractImageInformation(ImageInfo)

%This function returns a structure with the relevant information related to
%the acquisition parameters

%Imaging conditions
FrameInfo.LinesPerFrame=str2num(ExtractInformationField(ImageInfo,'state.acq.linesPerFrame='));
FrameInfo.PixelsPerLine=str2num(ExtractInformationField(ImageInfo,'state.acq.pixelsPerLine='));
FrameInfo.ZoomFactor=str2num(ExtractInformationField(ImageInfo,'state.acq.zoomFactor='));
FrameInfo.Rotation=str2num(ExtractInformationField(ImageInfo,'state.acq.scanRotation='));

% ES 2013-10-30: Compatibility with ScanImage 3.8
ScanImageVersionS = ExtractInformationField(ImageInfo, 'state.software.version=');
if strcmp(ScanImageVersionS(1:end-1), '3') % Referring to ScanImage 3.5.1
    FrameInfo.ScanAmplitudeX=str2num(ExtractInformationField(ImageInfo,'state.acq.scanAmplitudeX='));
    FrameInfo.ScanAmplitudeY=str2num(ExtractInformationField(ImageInfo,'state.acq.scanAmplitudeY='));
elseif strcmp(ScanImageVersionS(1:end-1), '3.8') % Referring to ScanImage 3.8
    FrameInfo.ScanAngleMultiplierFast=str2num(ExtractInformationField(ImageInfo,'state.acq.scanAngleMultiplierFast='));
    FrameInfo.ScanAngleMultiplierSlow=str2num(ExtractInformationField(ImageInfo,'state.acq.scanAngleMultiplierSlow='));
end

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



