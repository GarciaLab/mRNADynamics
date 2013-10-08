function output=ExtractInformationField(ImageInfo,Field)

%This function extracts a specific field form the ScanImage ImageInfo

%Find where the carriage returns are
cr=strfind(ImageInfo.ImageDescription,13);

PosField=strfind(ImageInfo.ImageDescription,Field);

%Find the closes carriage return to our field
[dummy,MinIndex]=min(abs(cr-PosField-length(Field)));

output=ImageInfo.ImageDescription(PosField+length(Field):cr(MinIndex));