function LIFExportMode_flatFieldImage(LIFMeta, Folder, OutputFolder, Prefix, PreferredFileForTest)
  %The flat field image can be in the folder with the data or in the folder corresponding to the date.
  D1 = dir([Folder, filesep, 'FF*.lif']);
  D2 = dir([Folder, filesep, '..', filesep, 'FF*.lif']);
  FFPaths = {};
  for i = 1:length(D1)
    FFPaths{end+1}=[Folder, filesep, D1(i).name];
  end
  for i = 1:length(D2)
    FFPaths{end+1} = [Folder, filesep, '..', filesep,D2(i).name];
  end

  %Go through the FF files and see which one matches the pixel size and image pixel number
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
    FilePath = [];
    if (~empty(PreferredFileNameForTest)) 
      FilePath = [Folder, filesep, PreferredFileNameForTest];
      disp(['Too many flat field images, using file name specified for testing', FilePath])
    else
      warning('Too many flat field images match the pixel and image size size')
      [FFFile,FFPath] =...
          uigetfile([Folder, filesep, '*.lif'], 'Select which flat field image to use');
      FilePath = [FFPath, FFFile];
    end;
    LIFFF = bfopen(FilePath);
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
end