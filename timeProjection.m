function [maxTProj, medianTProj] = timeProjection(Prefix,varargin)
% This function will return the max and median time projection of the movie.
% It will make TotalFrames, ZSlixes and NDigits 
%(which are all variables from CheckParticleTracking). 

%Parameters:
%Prefix: Prefix of the data set to analyze
%justnc13 : Only look at particles that show up in nc13 

%% Checking Varargin
justNC13 = 0;
if length(varargin)
    for i=1:length(varargin)
        if strcmpi(varargin{i}, 'justnc13')
            justNC13 = 1;
        end
    end
end

%% Information about the Folders and Frames
[~,~,DropboxFolder,~,~]=...
    DetermineLocalFolders(Prefix);
DataFolder=[DropboxFolder,filesep,Prefix];
FilePrefix=[DataFolder(length(DropboxFolder)+2:end),'_'];

if exist([DataFolder, filesep, 'FrameInfo.mat'])
    load([DataFolder, filesep, 'FrameInfo.mat'])
else
    error('noFrameInfo.mat found')
end 

TotalFrames = length(FrameInfo);
ZSlices = FrameInfo(1).NumberSlices + 2; %Note: Blank slides are included

if TotalFrames < 1E3
    NDigits = 3;
elseif TotalFrames < 1E4
    NDigits = 4;
end

if justNC13
    %Find the corresponding entry in the XLS file    
    [~,~,XLSRaw]=xlsread([DropboxFolder,filesep,'MovieDatabase.xlsx']);
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
    
    %Getting nc13 frame information
    nc13Column=find(strcmp(XLSRaw(1,:),'nc13'));
    nc14Column=find(strcmp(XLSRaw(1,:),'nc14'));
    nc13Frame = cell2mat(XLSRaw(XLSEntry,nc13Column));
    nc14Frame = cell2mat(XLSRaw(XLSEntry,nc14Column));
    frameRange = nc13Frame : nc14Frame; %Only doing projections for nc13
else
    frameRange = 1:TotalFrames; % doing projections using all frames
end

%% Z Projections
% The variables below will store the max and median Z projections 
maxZProjs = [];
medianZProjs = [];
for CurrentFrame = frameRange
    [maxZProjs(:,:,CurrentFrame), medianZProjs(:,:,CurrentFrame)]= ...
        zProjections(Prefix, CurrentFrame, ZSlices, NDigits);
%     figure(1)
%     imshow(maxZProjs(:,:,CurrentFrame),[0 80],'Border','Tight')
%     figure(2)
%     imshow(medianZProjs(:,:,CurrentFrame),[],'Border','Tight')
%     drawnow
%     disp(CurrentFrame)
end

%% Time Projections
% Taking the max and median with respect to the time axis (3) 
maxTProj = max(maxZProjs,[],3);
medianTProj = max(medianZProjs,3);
% figure(3)
% imshow(maxTProj,[0 80])
end