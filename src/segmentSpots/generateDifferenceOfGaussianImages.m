% Generates difference of Gaussian images
function [sigmas] = generateDifferenceOfGaussianImages(DogOutputFolder, pixelSize, customFilter, nCh, ExperimentType, coatChannel, numFrames, displayFigures, zSize, PreProcPath, Prefix, filterType, highPrecision)

  filterSize = round(2000 / pixelSize);
  sigmas = determineSigmas(customFilter, pixelSize);

  for channelIndex = 1:nCh
    % h=waitbar(0,'Generating DoG images');
    % (MT, 2018-02-11) Added support for lattice imaging, maybe 
    % temporary - FIX LATER

    if strcmpi(ExperimentType, 'inputoutput') || strcmpi(ExperimentType, 'lattice')
      nameSuffix= ['_ch', iIndex(coatChannel, 2)];
    else
      nameSuffix = ['_ch', iIndex(channelIndex, 2)];
    end
        
    for current_frame = 1:numFrames
      % waitbar(current_frame/numFrames,h);
      
      if displayFigures
      
        for zIndex = 1:zSize
          generateDoGs(DogOutputFolder, PreProcPath, Prefix, current_frame, nameSuffix, filterType, sigmas, filterSize,...
            highPrecision, zIndex, displayFigures);
        end
      
      else 
      
        parfor zIndex = 1:zSize   
          generateDoGs(DogOutputFolder, PreProcPath, Prefix, current_frame, nameSuffix, filterType, sigmas, filterSize,...
            highPrecision, zIndex, displayFigures);
        end
      
      end
    
    end

    % close(h);
  end 
    
end

function generateDoGs(DogOutputFolder, PreProcPath, Prefix, current_frame, nameSuffix, filterType, sigmas, filterSize, highPrecision, zIndex, displayFigures) 
  fileName = [PreProcPath, filesep, Prefix, filesep, Prefix, '_', iIndex(current_frame, 3), '_z',...
  iIndex(zIndex, 2), nameSuffix, '.tif'];

  im = double(imread(fileName));

  if strcmp(filterType, 'Difference_of_Gaussian')
    dog = filterImage(im, filterType, sigmas, filterSize);
    
    if highPrecision
      dog = (dog + 100) * 100;
    end
  
  else
    dog = filterImage(im, filterType, sigmas, []) + 100;
  end
  
  dog = padarray(dog(filterSize:end - filterSize - 1, filterSize:end - filterSize - 1), [filterSize, filterSize]);
  dog_name = ['DOG_', Prefix, '_', iIndex(current_frame, 3), '_z', iIndex(zIndex, 2), nameSuffix, '.tif'];
  dog_full_path = [DogOutputFolder, filesep, dog_name];
  imwrite(uint16(dog), dog_full_path)
  
  if displayFigures
    imshow(dog, [median(dog(:)), max(dog(:))]);
  end

end

function sigmas = determineSigmas(customFilter, pixelSize)
  %If customSigma is desired, prompts for the sigma1 and sigma2 values

  if ~customFilter

    %Initialize Difference of Gaussian filter parameters. filterSize >> sigma2
    %> sigma1
    sigma1 = pixelSize / pixelSize; %width of narrower Gaussian
    sigma2 = round(42000 / pixelSize); % width of wider Gaussian. AR 1/10/18: what is this number.

  end 

  sigmas = {sigma1, sigma2};

end 