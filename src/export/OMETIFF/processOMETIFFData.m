function FrameInfo = processOMETIFFData(rawDataFolder, OMETIFFFile, FrameInfo, Channel1, Channel2, Prefix, OutputFolder)
  disp('Processing OME-TIFF Data');
  disp(['RawData folder: ', rawDataFolder]);
  
  disp(['Channel 1: ', Channel1{1}]);
  disp(['Channel 2: ', Channel2{1}]);
  disp(['Output folder: ', OutputFolder]);

  omeTiffCompanionFile = [OMETIFFFile.folder, filesep, OMETIFFFile.name];
  disp(['OME-TIFF companion file: ', omeTiffCompanionFile]);
  
  % OME-TIFF XML companion file has this structure:
  %
  % <OME UUID="urn:uuid:b8cda420-61b6-4458-b5eb-e9eeb7fb697a"
  %    xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06"
  %    xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
  %    xsi:schemaLocation="http://www.openmicroscopy.org/Schemas/OME/2016-06
  %    http://www.openmicroscopy.org/Schemas/OME/2016-06/ome.xsd">
  %
  % <Image ID="Image:1" Name="Center_Single plane">                                                                                 
  %
  %     <Pixels ID="Pixels:1" SizeT="1"
  %         SizeX="2048" SizeY="2048" SizeZ="101" SizeC="2"
  %         PhysicalSizeX="0.173333333333333" PhysicalSizeY="0.173333333333333" PhysicalSizeZ="0.5"
  %         PhysicalSizeXUnit="µm" PhysicalSizeYUnit="µm" PhysicalSizeZUnit="µm"
  %         TimeIncrement="30" TimeIncrementUnit="s" DimensionOrder="XYZCT" Type="uint16">
  %
  %         <Channel ID="Channel:1" Name="Channel 1" Color="16711935" SamplesPerPixel="1"/>
  %         <Channel ID="Channel:2" Name="Channel 2" Color="-16711681" SamplesPerPixel="1"/>
  %
  %         <TiffData FirstT="0" FirstC="0" FirstZ="0" IFD="0" PlaneCount="101">
  %             <UUID FileName="t0001_Channel 1.tif">urn:uuid:66deddee-cf18-4928-83ea-8e9b95237f01</UUID>
  %         </TiffData>
  %         <TiffData FirstT="0" FirstC="1" FirstZ="0" IFD="0" PlaneCount="101">
  %             <UUID FileName="t0001_Channel 2.tif">urn:uuid:71017614-b2b2-4c79-bfbc-4a9ae06bc895</UUID>
  %         </TiffData>
  %     </Pixels>
  % </Image>
  % </OME>
  xml = xml2struct(omeTiffCompanionFile);
  
  disp('OME Pixel attributes:');
  pixelAttributes = xml.OME.Image.Pixels.Attributes;
  
  NChannels = str2double(pixelAttributes.SizeC);
  
  % TO-DO: get the series amount from...? Is it necessary?
  % NSeries = 1;
  
  NFrames = str2double(pixelAttributes.SizeT);
  
  % for FrameInfo
  LinesPerFrame = str2double(pixelAttributes.SizeX);
  PixelsPerLine = str2double(pixelAttributes.SizeY);
  NumberSlices = str2double(pixelAttributes.SizeZ);
  FileMode = 'OMETIFF'; %what if we fool the code here and use 'LIFExport' since it's similar?
  PixelSize = str2double(pixelAttributes.PhysicalSizeX);
  ZStep = str2double(pixelAttributes.PhysicalSizeZ);
  
  TimeIncrement = str2double(pixelAttributes.TimeIncrement);
  
  % zPosition = ? %not clear how to get it, for example in LIFExport it's
  % enclosed in a try/catch and only works for leica
  
  tifData = xml.OME.Image.Pixels.TiffData;
  tifDataCount = length(tifData);
  disp(['TIF files count: ', num2str(tifDataCount)]);
 
  % debug level - for now
  % loci.common.DebugTools.setRootLevel('INFO'); % commented out, for some
  % reason is not finding this bioformats class on the server.

  fprintf('Reading file: %s\n', omeTiffCompanionFile);
  TIFImages = bfopen(omeTiffCompanionFile);

  %Extract the metadata for each series
  TIFMeta = TIFImages{:, 4};

  for i = 1:sum(NFrames)
      
    FrameInfo(i).LinesPerFrame = LinesPerFrame;
    FrameInfo(i).PixelsPerLine = PixelsPerLine;
    FrameInfo(i).NumberSlices = NumberSlices;
    FrameInfo(i).FileMode = FileMode;
    FrameInfo(i).PixelSize = PixelSize;
    FrameInfo(i).ZStep = ZStep;
    FrameInfo(i).Time = (i - 1) * TimeIncrement; 
  end
  
  writeOMETifToOutputFolder(TIFImages, Prefix, NFrames, NChannels, NumberSlices, OutputFolder);
  
  disp('Finished processing OME-TIFF data');

end
