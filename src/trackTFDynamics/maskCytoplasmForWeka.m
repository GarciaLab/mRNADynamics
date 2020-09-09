function maskCytoplasmForWeka(Prefix, varargin)
%
% function maskCytoplasmForWeka(Prefix, varagin)
%
% DESCRIPTION
% This function takes the PreProcessed TIFs for the specified experiment 
% and applies the nuclear mask from the max projected histone channel to
% them to generate masked TIFs that can be fed into Weka like the normal
% PreProcessed TIFs.
% This function is useful if your movies have high levels of cytoplasmic
% fluorescence, especially very bright aggregates, which might make it
% difficult to train a classifier without excessive false positives.
%
% Note: This function only works with the new 3D TIFF stacks that are
%       output by the current version of the pipeline
%       ***IT IS NOT BACKWARDS COMPATIBLE***
%
%
% INPUT ARGUMENTS
% Prefix: prefix for this experiment
% 
%
% OPTIONS
% 'radiusScale', radiusScale: Multiplicative scaling factor by which the 
%                             radius of the nuclear mask is increased to 
%                             ensure the whole nucleus is captured. 
%                             Default radiusScale = 1.3
%
% 'includeChannels', channelsToMask: Array containing the channel number(s)
%                                    that you want to mask, e.g. for a
%                                    dataset with 3 channels channelsToMask
%                                    = [1,2] will mask the first and second
%                                    channels.
%                                    Defaults to all channels with
%                                    PreProcessed images.
%
% 'maskNormalizedImages': Masks previously normalized PreProcessed images, 
%                         which are stored in a subfolder called
%                         'PreProcessed\normalizedImages'
%
%
% OUTPUT
% normalizedFolder: path to the folder where the normalized movie frames
%                   (saved as 3D, xyz TIFF stacks, same as the input movie)
%                   are saved
% TIFF files containing a bleaching corrected movie
%
%
% Author (contact): Meghan Turner (meghan_turner@berkeley.edu)
% Created: 08/12/2020
% Last Updated: N/A
%

% Set user input option defaults
radiusScale = 1.3;
channelsToMask = [];
maskNormIm = false;
% Determine if user set non-default options
determineMaskCytoplasmOptions;

% Get needed info for this experiment
liveExperiment = LiveExperiment(Prefix);

nCh = numel(liveExperiment.spotChannels);
channels = liveExperiment.Channels;
channelNames = cell(1,nCh);

xDim = liveExperiment.xDim;
yDim = liveExperiment.yDim;
zDim = liveExperiment.zDim;
nFrames = liveExperiment.nFrames;

dropboxFolder = liveExperiment.userResultsFolder;
resultsFolder = liveExperiment.resultsFolder;
preProcFolder = liveExperiment.preFolder;

% Set folder paths for normalized vs unnormalized starting images
if maskNormIm
    maskPreProcFolder = [preProcFolder filesep 'normMaskImages'];
    preProcFolder = [liveExperiment.preFolder filesep 'normalizedImages'];
    maskSuffix = '_normMask';
    mkdir(maskPreProcFolder);
elseif ~maskNormIm
    maskPreProcFolder = [preProcFolder filesep 'maskedImages'];
    maskSuffix = '_mask';
    mkdir(maskPreProcFolder);
end

Ellipses = getEllipses(liveExperiment);

% Adjust channels to mask based on user input
if isempty(channelsToMask)
    nCh = numel(channels);
    channelsToMask = 1:nCh;
else
    nCh = numel(channelsToMask);
end

% Loop over all channels that need to be masked
rawImDir = cell(1,nCh);
for channel = channelsToMask
    channelSearchString = ['*ch0' num2str(channel) '.tif'];
    rawImDir{channel} = dir([preProcFolder filesep channelSearchString]);
    
    h = waitbar(1/nFrames, ['Masking raw images for Channel 0' num2str(channel) '...']);
    % Mask each frame in the movie
    for currFrame = 1:nFrames
        waitbar(currFrame/nFrames,h);
        %Get the raw image
        currImPath = [rawImDir{1,channel}(currFrame).folder, filesep, rawImDir{1,channel}(currFrame).name];
        imStack = loadTiffStack(currImPath); %this is the xDim x yDim x zDim image matrix

        % Make the cytoplasmic mask using the max-projected His stack
        ellipseFrame = Ellipses{currFrame};
        nuclearMask = makeNuclearMask(ellipseFrame, [yDim xDim], radiusScale);
        nuclearMask = nuclearMask >= 1;  %overlapping nuclei are annotated with a value of 2, but I want a binarized mask   
%         figure(1)
%         imshow(nuclearMask,[])
        
        % Create the new file inside the maskedImages folder
        nameSuffix = ['_ch',iIndex(channel,2)];
        imageName = [Prefix, '_', iIndex(currFrame,3), maskSuffix, nameSuffix, '.tif'];
        
        % Apply the z-projected nuclear mask to each z slice (this might be
        % a poor mask for the ends of the z stack -- maybe fix later)
        paddedZDim = zDim + 2;
        for currZSlice = 1:paddedZDim
            imSlice = imStack(:,:,currZSlice);
            maskedSlice = double(imSlice) .* nuclearMask;
            %save this masked z slice
            if currZSlice == 1
                imwrite(uint16(maskedSlice), [maskPreProcFolder, filesep, imageName]);
            else
                imwrite(uint16(maskedSlice),[maskPreProcFolder, filesep, imageName], 'WriteMode', 'append');
            end
        end
    end
    close(h)
    
end


%% Nested function to process user input options for the parent function
function determineMaskCytoplasmOptions
    for i = 1:numel(varargin)
        if ischar(varargin{i})
            % radiusScale option
            if strcmpi(varargin{i},'radiusScale')
                if (i+1 <= numel(varargin)) & isnumeric(varargin{i+1})
                    radiusScale = varargin{i+1};
                elseif ((i+1 > numel(varargin)) | ~isnumeric(varargin{i+1}))
                    error(['Input option ''' varargin{i}, ''' must be followed by a numeric variable.'])
                end
            % includeChannels option
            elseif strcmpi(varargin{i},'includeChannels')
                if (i+1 <= numel(varargin)) & isnumeric(varargin{i+1})
                    channelsToMask = varargin{i+1};
                elseif ((i+1 > numel(varargin)) | ~isnumeric(varargin{i+1}))
                    error(['Input option ''' varargin{i}, ''' must be followed by a numeric array containing the channel(s) to include.'])
                end
            % maskNormalizedImages option
            elseif strcmpi(varargin{i},'maskNormalizedImages')
                maskNormIm = true;
            % Notify user of invalid options
            else
                error([varargin{i}, ' is not a valid input option.'])
            end
        end
    end
end


end


