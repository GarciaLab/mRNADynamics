
function generateNuclearTiffStacks(Prefix)

% get basic project info 
liveExperiment = LiveExperiment(Prefix);

% make output director
nucleusStackDir = [liveExperiment.preFolder filesep 'nucleusStacks' filesep];
mkdir(nucleusStackDir);
nucleusProbDir = [liveExperiment.procFolder filesep 'nucleusProbabilityMaps' filesep];
mkdir(nucleusProbDir);

% load reference histogram for generating fake histone channel
% load('ReferenceHist.mat', 'ReferenceHist');
% Ellipses = getEllipses(liveExperiment);

% get nucleus channels
NuclearChannels = liveExperiment.nuclearChannels;

invertedChannelFlags = contains(liveExperiment.Channels(NuclearChannels),'inverted');
inputChannelFlags = ismember(NuclearChannels,liveExperiment.inputChannels);

% load frame information
FrameInfo = getFrameInfo(liveExperiment);

% get approximate time resolution
Tres = median(diff([FrameInfo.Time]));
frameSigma = 20/Tres;
timeIncrement = ceil(1 * 90 / Tres);

% calculate size of filter to use for spatial trend correction
PixelSizeXY = FrameInfo(1).PixelSize;
fieldSigma = 15 ./ PixelSizeXY;

% identify frames where stage position likely occured
seriesVec = [FrameInfo.Series];
seriesStartFrames = find([1 diff(seriesVec)]);

% determine subset of frames to use for segmentation
framesToUse = sort([seriesStartFrames,seriesStartFrames-1 length(FrameInfo)]);
framesToUseRaw = framesToUse(2:end);
framesToUse = [];

for i = 1:length(framesToUseRaw)-1
    framesToUse = [framesToUse framesToUseRaw(i):timeIncrement:framesToUseRaw(i+1)-1];
end
framesToUse(end+1) = length(FrameInfo);  
   
parfor frameIndex = 1:length(framesToUse)
    
    imStackFrame = zeros(FrameInfo(1).LinesPerFrame,FrameInfo(1).PixelsPerLine,FrameInfo(1).NumberSlices+2,length(NuclearChannels));
%     metric_vec = NaN(1,length(NuclearChannels));
    for ch = 1:length(NuclearChannels)
      
        % use sliding window averaging to filter out shot noise
        CurrentFrame = framesToUse(frameIndex);
        CurrentSeries = seriesVec(CurrentFrame);
        framesToLoad = CurrentFrame-ceil(frameSigma):CurrentFrame+ceil(frameSigma);
        framesToLoad = framesToLoad(framesToLoad>=1&framesToLoad<=length(FrameInfo));        
        framesToLoad = framesToLoad(seriesVec(framesToLoad)==CurrentSeries);
        frameWeights = exp(-.5*((framesToLoad - CurrentFrame)/frameSigma).^2);
        imStackInit = zeros(FrameInfo(1).LinesPerFrame,FrameInfo(1).PixelsPerLine,FrameInfo(1).NumberSlices+2);
        for i = 1:length(framesToLoad)
            imStackInit = imStackInit + frameWeights(i)*double(getMovieFrame(liveExperiment, framesToLoad(i), NuclearChannels(ch)));
        end
        imStackTemp = imStackInit(:,:,2:end-1,:)./sum(frameWeights);
        
        % suppress bright pixels
        px98 = prctile(imStackTemp(:),98);
        imStackTemp(imStackTemp>px98)=px98;  
        
        % invert if appropriate
        if invertedChannelFlags(ch)
          imStackTemp = px98 - imStackTemp;%mat2gray(imcomplement(imStackTemp));
%           imStackTemp = histeq(mat2gray(imStackTemp), ReferenceHist);
        end
        
        if inputChannelFlags(ch)
          % account for spatial heterogeneity in nuclei 
          refStack = imgaussfilt(imStackTemp,fieldSigma);
          imStackTemp = imStackTemp./ refStack;
        end
        
        imStackTemp = imStackTemp / std(imStackTemp(:));
        imStackTemp = imStackTemp - nanmean(imStackTemp(:)); %NL: this is likely futile               
        imStackFrame(:,:,2:end-1,ch) = (imStackTemp);%./ refStack .* normStack;                
        
        % use DoG to estimate noise
%         [thresh,~] = multithresh(imgaussfilt(imStackTemp(:,:,5),5));
%         imLB = bwlabel(imgaussfilt(imStackTemp(:,:,5),5)>=thresh);
%         metric_vec(ch) = 1 + abs(length(Ellipses{framesToUse(frameIndex)})-length(unique(imLB)))/length(Ellipses{framesToUse(frameIndex)});
    end
         
    % combie into a single  stack to use for segmentation
    if any(invertedChannelFlags)
        imStack = mat2gray(max(imStackFrame,[],4));
    else
        imStack = mat2gray(mean(imStackFrame,4));
    end
    
    % write to output directory    
    outName = [Prefix '_' sprintf('%03d',framesToUse(frameIndex)) '_ch00.tif'];
    imwrite(imStack(:,:,1),[nucleusStackDir outName]);   
    for i = 2:size(imStack,3)
        imwrite(imStack(:,:,i),[nucleusStackDir outName],'WriteMode','append');   
    end
    
end
%      

