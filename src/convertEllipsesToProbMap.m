function convertEllipsesToProbMap(Prefix,varargin)
% convertEllipsesToProbMap(Prefix,[Options])
%
% DESCRIPTION
% This script is for the case when we have a nice segmentation of nuclei,
% but bad tracking. So, we will use the Ellipses.mat to generate a fake
% probability map, which can be fed into the Tr2D for a better tracking.
%
% ARGUMENTS
% Prefix: Prefix of the data set to analyze
% [Any other parameters you have]
%
% OPTIONS
%
% CONTROLS
%
% OUTPUT
% Probability map : This should be fed into the Tr2D, thus should have the
% same name, file type as files saved in Tr2D directory.
% classificationImage0.tif 
%
% Author (contact): Yang Joon Kim (yjkim90@berkeley.edu)
% Created: 03/01/2019
% Last Updated: 03/01/2019
%
% Documented by: Yang Joon Kim (yjkim90@berkeley.edu)

%% Load the dataset
% We need nuclei segmentation(Ellipses.mat), pixel size, radius, etc.

  % Get the actual folder now that we have the Prefix
  [~, ~, DropboxFolder, ~, PreProcPath] = DetermineLocalFolders(Prefix);

  % What type of experiment are we dealing with? Get this out of MovieDatabase
  [~, ExperimentType, ~, ~, ~, ~, Channel1, Channel2, ~, ~, ~, ~, ~, ...
      nc9, nc10, nc11, nc12, nc13, nc14, ~] = getExperimentDataFromMovieDatabase(Prefix, DropboxFolder);

  % Set the destination folder (PreProcessedFolder)
  OutputFolder = [PreProcPath, filesep, Prefix];

  % Load the information about this image
  % Check if we have FrameInfo otherwise try to get the information straight
  % from the file.

  %SearchRadius = ceil(SearchRadiusMicrons / PixelSize);
    load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat']);
    nPixels_y = FrameInfo.LinesPerFrame; % number of pixels in y axis
    nPixels_x = FrameInfo.PixelsPerLine; % number of pixels in x axis
    nFrames = length(FrameInfo);
    ConversionUnit = FrameInfo.PixelSize; % um/pixel
    
  % Load the nuclei segmentation information
    load([DropboxFolder, filesep, Prefix, filesep, 'Ellipses.mat'], 'Ellipses');
    Nuclei = Ellipses;
    
    
    


%% Make a binary image of nuclei/cytoplasm
% Define the Probability map
ProbMap = zeros(nPixels_y, nPixels_x, nFrames);
    
    % for loop over frames
    for i=1:length(Nuclei)
        % Initialize the probability map for each frame
        probmap_temp = zeros(nPixels_y, nPixels_x);
        % At each frame, extract the info about nuclei position, and
        % radius. Then convert the zeros in the Prob map to ones.
        Nuclei_temp_currentFrame = Nuclei{i,1};
            for j=1:length(Nuclei_temp_currentFrame)
                % Move through between nuclei (segmented) and get the info
                % about their x, y positions, and also radius. 
                % Note1 : I'm making assumption that the nuclei are usually
                % circles, not ellipsoid, thus, I can just grab the third
                % value for the radius.
                % Note2 : The indexing is from the left top (Origin), thus
                % I need to convert that to actual position in the
                % prob.map.
                
                xPos = Nuclei_temp_currentFrame(j,1); % Center's x position
                yPos = Nuclei_temp_currentFrame(j,2); % Center's y position
                % Radius Adjustment
                % For many cases, the diameter is pretty large value, thus
                % I'll divide it with 2, for now.
                adjustment_factor_diameter = 3;
                Diameter = ceil((Nuclei_temp_currentFrame(j,3)-1)./ConversionUnit) /adjustment_factor_diameter ; % um divided by um/pixel, thus it's now pixels unit.
                
                % Make the nuclei region (pixels') prob map to be 1.
                for k = 1:nPixels_x
                    for l = 1:nPixels_y
                       distance = (k-xPos).^2 + (l-yPos).^2 ;
                       if distance<(Diameter/2)^2
                            probmap_temp(l,k) = 1;
                       end
                    end
                end
            end
         ProbMap(:,:,i) = probmap_temp;
    end
    
 %% Check the Probability map
% for i=1:nFrames
%     imshow(ProbMap(:,:,i),[])
%     pause
% end
%% Save
imwrite(ProbMap(:,:,1),[OutputFolder,filesep,...
                'tr2dProject',filesep,'segmentation',filesep,'weka',filesep,'classificationImage0.tif']); 
for i=2:nFrames
    imwrite(ProbMap(:,:,i),[OutputFolder,filesep,...
                'tr2dProject',filesep,'segmentation',filesep,'weka',filesep,'classificationImage0.tif'],'WriteMode','append'); 
end
end