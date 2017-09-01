function arff = makeArff(Prefix)
% makeArff(groundTruth)
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
% None.
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
Spots = load([DropboxFolder,filesep,Prefix,filesep,'Spots.mat']);
numFrames = length(Spots.Spots);

%Reading in the filters from a file that happens to contain all of them.
in_path = 'C:\Users\ArmandoReimer\Desktop\GermLineCloneP2P All filters data.arff';
fin = fileread(in_path);
filters = regexp(fin,'(?<=@attribute )\S*(?= )', 'match');
filters = filters(2:end-1);
numFilters = length(filters);

for i = 1:numFrames
    imPath= [PreProcPath,filesep,Prefix,filesep,'stacks', filesep, iIndex(i,3),'.tif'];
    truthPath = [FISHPath,filesep,Prefix,'_', filesep, 'binary_masks', filesep, 'binary_stack_',Prefix,'_',iIndex(i,3),'.tiff'];
    im = zeros(rows,cols,zSlices);
    groundTruth = zeros(rows,cols,zSlices);
    for j = 1:zSlices
        im(:,:,j) = imread(imPath, j);
        groundTruth(:,:,j) = imread(truthPath, j);
    end
    if i==1
        fCell = [{'original'; im},{'class'; groundTruth}];
        for o = 1:numFilters
            filterName = filters{o};
            filterType = regexp(filterName,'\D*(?=_\d)', 'match');
            filterType = filterType{1};
            sigmas = regexp(filterName,'(?<=_)\d(?=.)','match'); 
            if ~strcmp(filterName, 'original') && ~strcmp(filterName, 'class')
                f = filterImage(im, filterType, sigmas);      
                fCell = [fCell, {filterName; f}];
            end
        end
    else
        for o = 1:numFilters
            filterName = filters{o};
            filterType = regexp(filterName,'\D*(?=_\d)', 'match');
            filterType = filterType{1};
            sigmas = regexp(filterName,'(?<=_)\d(?=.)','match'); 
            if ~strcmp(filterName, 'original') && ~strcmp(filterName, 'class')
                f = filterImage(im, filterType, sigmas);      
                [~,x] = find(strcmp(fCell,'filterName'));
                if x==1
                    fCell{end+1, x} = f; 
                else
                    fCell{end, x} = f; 
                end                   
            end
        end
    end
    fCellLength = size(fCell,2);
    for p = 1:fCellLength
        arff{1, p} = fCell{1, p};
    end
    n = 2;
    for k = 1:rows
        for m = 1:cols
            for p = 1:fCellLength
                arff{n,p} = fCell{2, p}(k,m);                    
                n = n + 1;
            end
        end
    end
end

T = cell2table(arff);
writetable(T,[DropboxFolder,filesep,Prefix,filesep, Prefix,'.arff'])