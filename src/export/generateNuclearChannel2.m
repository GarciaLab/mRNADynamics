% This will check all channels, then grab channels with ":Nuclear" in the
% string, then add the Projected images from multiple channels.
% For MovieDatabase Channels, we should put ":Nuclear" or
% "invertedNuclear" for those channels to be recognized for histone
% channel generation.
%
function Projection = generateNuclearChannel2(ProjectionType, Channels, ReferenceHist, movieMat, frame)

yDim = size(movieMat, 1);
xDim = size(movieMat, 2);
NSlices = size(movieMat, 3);

Channel1=Channels{1}; Channel2 = Channels{2}; Channel3 = Channels{3};

% Check how many channels have ":Nuclear" in the MovieDatabase.csv
NuclearChannels = [contains(Channel1, 'Nuclear', 'IgnoreCase', true),...
    contains(Channel2, 'Nuclear', 'IgnoreCase', true),...
    contains(Channel3, 'Nuclear', 'IgnoreCase', true)];
nNuclearChannels = sum(NuclearChannels);

InvertedChannels = [contains(Channel1, 'inverted', 'IgnoreCase', true), ...
    contains(Channel2, 'inverted', 'IgnoreCase', true), ...
    contains(Channel3, 'inverted', 'IgnoreCase', true)];
AllChannels = {Channel1,Channel2,Channel3}; %Cell array of channel labels

ProjectionTemp = zeros(yDim, xDim, nNuclearChannels);

if nNuclearChannels ~= 0
    
    for ChannelIndex = 1:nNuclearChannels
        
        
        % Find the corresponding channel
        AllNuclearChannels = find(~ NuclearChannels == 0, nNuclearChannels);
        nuclearChannel = AllNuclearChannels(ChannelIndex);
        
        % For all 'nuclear' channels, generate HisSlices, and do projection
        
        HisSlices = movieMat(:, :, :, frame, nuclearChannel);
        
        ProjectionTemp(:, :, ChannelIndex) = calculateProjection(ProjectionType, NSlices, HisSlices);
        
        if InvertedChannels(nuclearChannel) == 1
            ProjectionTemp(:, :, ChannelIndex) =...
                imcomplement(ProjectionTemp(:, :, ChannelIndex));
        end
        
        ProjectionTemp(:, :, ChannelIndex) = 256*...
            histeq(mat2gray(ProjectionTemp(:, :, ChannelIndex)), ReferenceHist);
        
        % Check if we are excluding this frame from this nuclear channel
        excludeFrames = 0;
        
        if and(contains(AllChannels{nuclearChannel},'['),contains(AllChannels{nuclearChannel},']'))
            excludeFrames = 1;
        end
        
%         if excludeFrames
%             framesToExcludeAll = extractBetween(AllChannels{nuclearChannel},'[',']');
%             framesToExcludeAll = framesToExcludeAll{1};
%             
%             %Split at commas to separate into multiple sets of exclusion frames
%             framesToExcludeAll = strsplit(framesToExcludeAll,',');
%             
%             %Save both the series and frames to exclude per series
%             firstFrame = [];
%             lastFrame = [];
%             seriesToExclude = [];
%             for i = 1:length(framesToExcludeAll)
%                 seriesAndFrames = strsplit(framesToExcludeAll{i},'_');
%                 seriesToExclude(end+1) = str2double(seriesAndFrames{1});
%                 framesToExclude = seriesAndFrames{2};
%                 colIndex = strfind(framesToExclude,':');
%                 firstFrame(end+1) = str2double(framesToExclude(1:(colIndex-1)));
%                 lastFrame(end+1) = str2double(framesToExclude((colIndex+1):end));
%             end
%             
%             %If this frame is within the excluded series and frames, throw it out
%             
%             for f = 1:length(seriesToExclude)
%                 frameWindow = firstFrame(f):lastFrame(f);
%                 if and(~isempty(find(ismember(framesIndex,frameWindow), 1)),...
%                         seriesIndex == seriesToExclude(f))
%                     ProjectionTemp(:, :, ChannelIndex) = nan(size(ProjectionTemp,1),size(ProjectionTemp,2));
%                     disp(['Excluding frame ',num2str(framesIndex),' of series ',...
%                         num2str(seriesIndex),' from nuclear projection of channel ',...
%                         num2str(nuclearChannel)]);
%                 end
%             end
%         end
    end
    
    % Getting average of all Projections
    if nNuclearChannels > 1
        Projection = nanmean(ProjectionTemp, 3);
    elseif nNuclearChannels == 1
        Projection = ProjectionTemp;
    end
   

Projection = uint8(Projection);


% imwrite(Projection, [OutputFolder, filesep, Prefix, '-His_', iIndex(numberOfFrames, 3), '.tif']);

end



