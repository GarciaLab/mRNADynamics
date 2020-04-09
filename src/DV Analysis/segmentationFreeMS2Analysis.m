function [h, ax] = segmentationFreeMS2Analysis(Prefix, nBins, varargin)

%color is the bin color accepted as a string within one element cell.
%use an empty cell as default. accepted colors- 'red', 'yellow', 'cyan', 'magenta',
%'lightBlue'

scale = 1;
color = {'red'}; %just using red as the default
optionalResults = '';
ax = [];
nuclearMaskFlag = false;

%area for improvement: segment the nuclei and only include dog values from
%within nuclei 
for i = 1:length(varargin)
    if strcmpi(varargin{i}, 'color')
        color = varargin{i+1};
    elseif strcmpi(varargin{i}, 'scale')
        scale = varargin{i+1};
    elseif strcmpi(varargin{i}, 'ax')
        ax = varargin{i+1};
    elseif strcmpi(varargin{i}, 'nuclearMask')
        nuclearMaskFlag = true;
    end
end

if isempty(ax)
    fig = figure();
    ax = axes(figure);
end

[~,ProcPath,DropboxFolder] = readMovieDatabase(Prefix, optionalResults);

if nuclearMaskFlag
    load([DropboxFolder, filesep, Prefix, filesep, 'Ellipses.mat'], 'Ellipses');
end

load([DropboxFolder, filesep, Prefix, filesep, 'FrameInfo.mat']);


nFrames = length(FrameInfo);
nSlices = FrameInfo(1).NumberSlices;

dogDirPath = [ProcPath,filesep,Prefix,'_\dogs'];
dogDir = dir([dogDirPath, filesep,'DOG*']);
dogDir= {dogDir.name};
numIm = length(dogDir);

if contains(dogDir{1}, 'mat')
    saveType = 'mat';
elseif contains(dogDir{1}, 'tif')
    saveType = 'tif';
end

vals = [];

% try
%     parpool(20);
% end

for i = 1:numIm
    
    current_frame = floor(i/nSlices)+1;
    
    if strcmpi(saveType, 'tif')
        dog = imread([dogDirPath, filesep, dogDir{i}]);
    elseif strcmpi(saveType, 'mat')
        dog = load([dogDirPath, filesep, dogDir{i}]);
        dog = dog.plane;
    end
    
    
    if nuclearMaskFlag
         nuclearMask = ones(size(dog, 1), size(dog, 2));
        if current_frame <= nFrames
            ellipsesFrame = Ellipses{current_frame};
            nuclearMask = makeNuclearMask(ellipsesFrame, [size(dog,1), size(dog,2)], 'radScale', 1);
        end
    
        dog = dog.*nuclearMask;
    end
    
    vals(i) = log10(max(dog(:))+1);
    
    %imshow(dog, []);
    
end

vals(vals==0) = NaN;
h = histogram(ax, vals,nBins,'Normalization','pdf', 'facealpha', .6);
set(ax,'YScale','log');
xlabel(ax, 'log(max DoG intensity + 1) (au)');
ylabel(ax, 'frequency');
title(Prefix, 'Interpreter', 'none');
standardizeFigure(ax, [], color{1}, 'fontSize', 14);
hold on


end