function desiredTimeProjection = timeProjection(Prefix,varargin)
% timeProjection(Prefix, [Options])
%
% DESCRIPTION
% This function will calculate the max time projection of the movie.
% By default it will make the max time projection using the max z
% projections.
%
% ARGUEMENTS
% Prefix: Prefix of the data set to analyze
% DropboxFolder: Dropbox Folder the data set is in
%
% OPTIONS
% 'median': Make the max time projection usin the median z projections
% 'nc', NC: Choose which nc will be inluded in the timeProjection where NC
%           must be a number or a vector of two numbers (the range of nc that
%           will be included). 
%           (Ex. NC = [12 14] will include nc12, nc13 and nc14)
% OUTPUT
% This will return the max time projection made from either the 
% max z projections or median z projections.
%
% Author (contact): Emma Luu (emma_luu@berkeley.edu)
% Created: 06/16/2017
% Last Updated: 06/22/2017
%
% Documented by: Emma Luu (emma_luu@berkeley.edu)


%% Checking Varargin
h = waitbar(0,'Checking Varargin');
useMedian = 0;
ncRange = 0;
if length(varargin)
    
    for i=1:length(varargin)
        if strcmpi(varargin{i}, 'median')
            useMedian = 1;
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
[~,~,DropboxFolder,~,~]= DetermineLocalFolders(Prefix);
[~,~,DefaultDropboxFolder,~,~]= DetermineLocalFolders;
DataFolder=[DropboxFolder,filesep,Prefix];
FilePrefix=[DataFolder(length(DropboxFolder)+2:end),'_'];

if exist([DataFolder, filesep, 'FrameInfo.mat'])
    load([DataFolder, filesep, 'FrameInfo.mat']);
else
    error('noFrameInfo.mat found')
end

totalFrames = length(FrameInfo);
zSlices = FrameInfo(1).NumberSlices + 2; %Note: Blank slides are included

if totalFrames < 1E3
    NDigits = 3;
elseif totalFrames < 1E4
    NDigits = 4;
end

if ncRange
    %Find the corresponding entry in the XLS file    
    [~,~,XLSRaw]=xlsread([DefaultDropboxFolder,filesep,'MovieDatabase.xlsx']);
    DataFolderColumn=find(strcmp(XLSRaw(1,:),'DataFolder'));
    Dashes=findstr(FilePrefix,'-');

    if (~isempty(findstr(FilePrefix,'Bcd')))&(isempty(findstr(FilePrefix,'BcdE1')))&...
            (isempty(findstr(FilePrefix(1:end-1),'NoBcd')))&...
            (isempty(findstr(FilePrefix(1:end-1),'Bcd1x')))
        warning('This step in CheckParticleTracking will most likely have to be modified to work')
        XLSEntry=find(strcmp(XLSRaw(:,DataFolderColumn),...
            [Date,'\BcdGFP-HisRFP']));
    else
        XLSEntry=find(strcmp(XLSRaw(:,DataFolderColumn),...
            [FilePrefix(1:Dashes(3)-1),'\',FilePrefix(Dashes(3)+1:end-1)]));
        if isempty(XLSEntry)
            XLSEntry=find(strcmp(XLSRaw(:,DataFolderColumn),...
                [FilePrefix(1:Dashes(3)-1),'/',FilePrefix(Dashes(3)+1:end-1)]));
            if isempty(XLSEntry)
                error('Could not find data set in MovieDatabase.XLSX. Check if it is defined there.')
            end
        end
    end
    
    %Getting nc range frame information
    ncStartColumn=find(strcmp(XLSRaw(1,:),startNC));
    ncStartFrame = cell2mat(XLSRaw(XLSEntry,ncStartColumn));
    
    if ncStartFrame == 0
        ncStartFrame = 1;
    end
    
    %(6/22) EL: Not sure how to better handle this sitatuion
    if strcmp('nc15',endNC) 
        ncEndFrame = totalFrames;
    else
        % Find the first frame of the next nc and make the previous frame 
        % the last frame looked at
        ncEndColumn=find(strcmp(XLSRaw(1,:),endNC));
        ncEndFrame = cell2mat(XLSRaw(XLSEntry,ncEndColumn))-1; 
    end
    frameRange = ncStartFrame : ncEndFrame; % desired frame range
else
    frameRange = 1:totalFrames; % doing projections using all frames
end

%% Z Projections
% The variables below will store the max and median Z projections 
waitbar(0,h,'Starting Z Projections');
desiredZProjs = zeros(FrameInfo(1).LinesPerFrame, FrameInfo(1).PixelsPerLine, zSlices);
framesCompleted = 0;
if useMedian
    for CurrentFrame = frameRange
        [~,desiredZProjs(:,:,CurrentFrame)] = ...
            zProjections(Prefix, CurrentFrame, zSlices, NDigits);
        framesCompleted = 1 + framesCompleted;
        waitbar(framesCompleted/length(frameRange),h,'Making Z Median Projections');
    end
else
    for CurrentFrame = frameRange
        [desiredZProjs(:,:,CurrentFrame),~]= ...
            zProjections(Prefix, CurrentFrame, zSlices, NDigits);
        framesCompleted = 1 + framesCompleted;
        waitbar(framesCompleted/length(frameRange),h,'Making Z Max Projections');
    end
end
%% Time Projections
% Taking the max and median with respect to the time axis (3) 
desiredTimeProjection = max(desiredZProjs,[],3);
% figure(3)
% imshow(desiredTimeProjection,[0 80]);
close(h)
end