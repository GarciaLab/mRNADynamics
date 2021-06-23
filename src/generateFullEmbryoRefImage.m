% Edited from processFiducialChannel.m
% by Yang Joon Kim (yjkim90@berkeley.edu) on 10/23/2018
% This will check all channels, then grab channels with ":Nuclear" in the
% string, then add the Projected images from multiple channels.
% For MovieDatabase Channels, we should put ":Nuclear" or
% "invertedNuclear" for those channels to be recognized for histone
% channel generation.
function Projection = generateFullEmbryoRefImage(...
                        fullEmbryoImages, NSlices,...
                        NChannels,...
                        Channels)

Channel1=Channels{1}; Channel2 = Channels{2}; Channel3 = Channels{3};

% load reference histogram
load('ReferenceHist.mat', 'ReferenceHist');    

% calculate slice indices 
lowerSlice = round(NSlices * .33);
upperSlice = round(NSlices * .66);

% Check how many channels have ":Nuclear" in the MovieDatabase.csv
NuclearChannels = [contains(Channel1, 'Nuclear', 'IgnoreCase', true),...
    contains(Channel2, 'Nuclear', 'IgnoreCase', true),...
    contains(Channel3, 'Nuclear', 'IgnoreCase', true)];
nNuclearChannels = sum(NuclearChannels);

InvertedChannels = [contains(Channel1, 'inverted', 'IgnoreCase', true), ...
    contains(Channel2, 'inverted', 'IgnoreCase', true), ...
    contains(Channel3, 'inverted', 'IgnoreCase', true)];

Projection = [];
ProjectionTemp = zeros(size(fullEmbryoImages{1,1},1),size(fullEmbryoImages{1,1},2),nNuclearChannels);

if nNuclearChannels ~= 0
    
    for ChannelIndex = 1:nNuclearChannels
                
        % Find the corresponding channel
        AllNuclearChannels = find(~ NuclearChannels == 0, nNuclearChannels);
        nuclearChannel = AllNuclearChannels(ChannelIndex);
        
        
        imageList = nuclearChannel:NChannels:size(fullEmbryoImages,1);
        nuclearStack = zeros(size(fullEmbryoImages{1,1},1),size(fullEmbryoImages{1,1},2),NSlices);
        iter = 1;
        for imagesIndex = imageList
            nuclearStack(:, :, iter) = fullEmbryoImages{imagesIndex,1};
            iter = iter+1;
        end
        
        % calculate projection
        
        ProjectionTemp(:, :, ChannelIndex) = sum(nuclearStack(:,:,lowerSlice:upperSlice),3);
        
        % we need to generate an embryo mask
        im_temp = ProjectionTemp(:, :, ChannelIndex);
        ft_size = round(size(im_temp,1)/50);
        im_temp = imgaussfilt(im_temp,ft_size);
        thresh = multithresh(im_temp);
        im_mask = im_temp>thresh;
        
        % Think about "invertedNuclear", for example, MCP-mCherry, then
        % invert the ProjectionTemp using imcomplement
        if InvertedChannels(nuclearChannel) == 1               
            % generate complement
            cp_temp = ProjectionTemp(:, :, ChannelIndex);
            cp_temp(~im_mask) = max(cp_temp(:));
            cp_temp = imcomplement(cp_temp);
            % set outside to 0            
            ProjectionTemp(:, :, ChannelIndex) = cp_temp;
        end
        % normalize
        ProjectionTemp(:, :, ChannelIndex) = mat2gray(ProjectionTemp(:, :, ChannelIndex));
%         ProjectionTemp(:, :, ChannelIndex) = ProjectionTemp(:, :, ChannelIndex) * 255;        
        
    end
    
    % Getting average of all Projections    
    Projection = nanmean(ProjectionTemp, 3);    


end



