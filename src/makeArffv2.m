function arff = makeArff(Prefix,varargin)

% fCell = evalin('base', 'fCell'); %for debugging the arffmaking
% makeArff(Prefix)
%
% DESCRIPTION
% Takes ground truth spot/pixel data curated from CheckParticleTracking and
% makes a .arff file compatible with Weka classification.
%
% ARGUMENTS
% 'Prefix': Useful for getting pixel information from FrameInfo not
%           present in the 'groundTruth' structure.
%
% OPTIONS
% 'numFrames': Optionally only take pixels from the first numFrames frames
% of the movie. Useful for debugging. 
%               
% OUTPUT
% 'arff': An arff-like matrix that can be fed into Weka after being saved as
% .arff or .csv
%
% Author (contact): Armando Reimer (areimer@berkeley.edu)
% Created: 7/14/2017
% Last Updated: 7/14/2017
%
% Documented by: Armando Reimer (areimer@berkeley.edu)

[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath] = ...
    DetermineLocalFolders(Prefix);
frameInfo = load([DropboxFolder,filesep, Prefix, filesep, 'FrameInfo.mat']);
frameInfo = frameInfo.FrameInfo;
rows = frameInfo(1).LinesPerFrame;
cols = frameInfo(1).PixelsPerLine;
zSlices = frameInfo(1).NumberSlices;
stackSize = rows*cols*zSlices;
Spots = load([DropboxFolder,filesep,Prefix,filesep,'Spots.mat']);
numFrames = length(Spots.Spots);

%Read in optional parameters
for i=1:length(varargin)
    if strcmp(varargin{i},'numFrames')
        numFrames = varargin{i+1};
    end
end

%Reading in the filters from a file that happens to contain all of them.
in_path = 'C:\Users\ArmandoReimer\Desktop\GermLineCloneP2P All filters data.arff';
fin = fileread(in_path);
filters = regexp(fin,'(?<=@attribute )\S*(?= )', 'match');
filters = filters(2:end-1);
numFilters = length(filters);
arfftemp = zeros(stackSize*(numFilters+2)*numFrames,numFilters+2);
n = 0;
for i = 1:numFrames
    imPath= [PreProcPath,filesep,Prefix,filesep,'stacks', filesep, iIndex(i,3),'.tif'];
    truthPath = [FISHPath,filesep,Prefix,'_', filesep, 'binary_masks', filesep, 'binary_stack_',Prefix,'_',iIndex(i,3),'.tiff'];
    im = zeros(rows,cols,zSlices);
    groundTruth = zeros(rows,cols,zSlices);
    
    for j = 2:zSlices+1 %we don't want the top and bottom slices
        im(:,:,j-1) = imread(imPath, j);
        groundTruth(:,:,j-1) = imread(truthPath, j);
    end
    
    arfftemp(n+1:stackSize, 1) = im(:)';
    n = n + stackSize;
    arfftemp(n+1:stackSize, 2) = groundTruth(:)';
    n = n + stackSize;
    
%   parfor o = 1:numFilters %just a reminder this can be parfored
    for o = 1:numFilters
        filterName = filters{o};
        filterType = regexp(filterName,'\D*(?=_\d)', 'match');
        filterType = filterType{1};
        sigmas = regexp(filterName,'(?<=_)\d(?=.)','match'); 
        if ~strcmp(filterName, 'original') && ~strcmp(filterName, 'class')
            f = filterImage(im, filterType, sigmas);
            arfftemp(n:stackSize, o+2) = f(:)'; %note for armando- this needs to be filled in so i can eliminate fcell entirely.
            %i still need to figure out how to define
            %lastnonzeroentry and filterIndex
            n = n + stackSize;
        end
    end
end

csvwrite('C:\Users\ArmandoReimer\Desktop\arffout.arff', arfftemp);
end
