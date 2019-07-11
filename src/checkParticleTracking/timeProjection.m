function desiredProjection = timeProjection(Prefix,currentChannel,FrameInfo, DropboxFolder,PreProcPath, varargin)
% timeProjection(Prefix, currentChannel, [Options])
%
% DESCRIPTION
% This function can calculate the max time projection of the movie.
% By default it will make the max time projection using the max z
% projections. There is an option to only return either the max or median z
% projection of all the frames.
%
% ARGUEMENTS
% Prefix: Prefix of the data set to analyze
% currentChannel: The current channel of the spots. It is assumed that this
%                 is a number value.
%
% OPTIONS
% 'medianTimeProjection': Make the max time projection usin the median z projections
% 'nc', NC: Choose which nc will be inluded in the timeProjection where NC
%           must be a number or a vector of two numbers (the range of nc that
%           will be included). 
%           (Ex. NC = [12 14] will include nc12, nc13 and nc14)
% 'onlyMaxZProjections': Only makes max Z projections for all frames
% 'onlyMedianZProjections': Only makes median Z projections for all frames
% 
% OUTPUT
% By default, this will return the max time projection made from the 
% max z projections. It can also make the max time projection using median 
% z projections. It can also return only the max or median projection of
% all frames.
%
% Author (contact): Emma Luu (emma_luu@berkeley.edu)
% Created: 06/16/2017
% Last Updated: 11/12/2017
% 
% Documented by: Emma Luu (emma_luu@berkeley.edu)


%% Checking Varargin
h = waitbar(0,'Checking Varargin');
useMedian = 0; % Indicates whether to use median z projection for time projection
ncRange = 0; % Indicates whether to only do projection for a certain nc range
maxZOnly = 0; % Indicates whether to only make max Z projections
medianZOnly = 0; % Indicates whether to only make median Z projections

if length(varargin)
    for i=1:length(varargin)
        if strcmpi(varargin{i}, 'medianTimeProjection')
            useMedian = 1;
        elseif strcmpi(varargin{i}, 'onlyMaxZProjections')
            maxZOnly = 1;
        elseif strcmpi(varargin{i}, 'onlyMedianZProjections')
            medianZOnly = 1;
        end
        
        if strcmpi(varargin{i},'nc') % checking for the desired nc range
            ncRange = 1;
            if length((varargin{i+1})) == 2
                startNC = ['nc' num2str(varargin{i+1}(1))];
                endNC = ['nc' num2str(varargin{i+1}(2) +1)];% Not including the next nc
            else
                startNC = ['nc' num2str(varargin{i+1})]; 
                endNC = ['nc' num2str(varargin{i+1} + 1)]; % Not including the next nc
            end
        end
    end
end

%% Information about the Folders and Frames 
% Define and make them here instead of passing these from CheckParticleTracking
waitbar(0,h,'Getting Relevant Folder Information');

totalFrames = length(FrameInfo);
zSlices = FrameInfo(1).NumberSlices + 2; %Note: Blank slides are included

if totalFrames < 1E3
    NDigits = 3;
elseif totalFrames < 1E4
    NDigits = 4;
end

if ncRange
    
    %Getting nc range frame information
    %determine ncStarFrame dinamically from 'nc?' arguments in varargin above
    ncStartFrame = eval(startNC);
    
    if ncStartFrame == 0
        ncStartFrame = 1;
    end
    
    %(6/22) EL: Not sure how to better handle this sitatuion
    if strcmp('nc15',endNC) 
        ncEndFrame = totalFrames;
    else
        % Find the first frame of the next nc and make the previous frame 
        % the last frame looked at
        %determine ncEndFrame dinamically from 'nc?' arguments in varargin above
        ncEndFrame = eval(endNC); 
    end
    frameRange = ncStartFrame : ncEndFrame; % desired frame range
else
    frameRange = 1:totalFrames; % doing projections using all frames
end

%% Z Projections
% The variables below will store the max and median Z projections 
waitbar(0,h,'Starting Z Projections');
framesCompleted = 0;
if useMedian || medianZOnly
    desiredZProjs = zeros(FrameInfo(1).LinesPerFrame, FrameInfo(1).PixelsPerLine, zSlices);
    for CurrentFrame = frameRange
        desiredZProjs(:,:,CurrentFrame) = ...
            zProjections(Prefix, currentChannel, CurrentFrame, zSlices, NDigits,DropboxFolder,PreProcPath, FrameInfo, 'median');
        framesCompleted = 1 + framesCompleted;
        waitbar(framesCompleted/length(frameRange),h,'Making Z Median Projections');
    end
else
    im = zeros(FrameInfo(1).LinesPerFrame, FrameInfo(1).PixelsPerLine, 'uint16');
    desiredProjection = im;
    for CurrentFrame = frameRange
        for currentZ = 2:zSlices-1
            im = imread([PreProcPath,filesep,Prefix,filesep,...
                    Prefix,'_',iIndex(currentFrame,nDigits),'_z',iIndex(currentZ,2),'_ch0', num2str(currentChannel) ,'.tif']);
            desiredProjection = max(desiredProjection, im);
        end
        framesCompleted = 1 + framesCompleted;
%         plot(desiredZProjs(:,:,CurrentFrame))
%         title(['Frame: ' num2str(CurrentFrame)])
        waitbar(framesCompleted/length(frameRange),h,'Making Z Max Projections');
    end
end
%% Time Projections
% Taking the max and median with respect to the time axis (3) 

% if maxZOnly || medianZOnly
%     desiredProjection = desiredZProjs;
% else
%     desiredProjection = max(desiredZProjs,[],3);
% end
% figure()
% imshow(desiredTimeProjection,[0 80])
close(h)
end