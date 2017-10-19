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
arfftemp = []; 
for i = 1:numFrames
    imPath= [PreProcPath,filesep,Prefix,filesep,'stacks', filesep, iIndex(i,3),'.tif'];
    truthPath = [FISHPath,filesep,Prefix,'_', filesep, 'binary_masks', filesep, 'binary_stack_',Prefix,'_',iIndex(i,3),'.tiff'];
    im = zeros(rows,cols,zSlices);
    groundTruth = zeros(rows,cols,zSlices);
    for j = 2:zSlices+1 %we don't want the top and bottom slices
        im(:,:,j-1) = imread(imPath, j);
        groundTruth(:,:,j-1) = imread(truthPath, j);
    end
    if i==1
        fCell = [{'original'; im},{'class'; groundTruth}];
%         parfor o = 1:numFilters
        for o = 1:numFilters
            filterName = filters{o};
            filterType = regexp(filterName,'\D*(?=_\d)', 'match');
            filterType = filterType{1};
            sigmas = regexp(filterName,'(?<=_)\d(?=.)','match'); 
            if ~strcmp(filterName, 'original') && ~strcmp(filterName, 'class')
                f = filterImage(im, filterType, sigmas);
                arfftemp(lastnonzeroentry:stackSize, filterIndex) = f(:)'; %note for armando- this needs to be filled in so i can eliminate fcell entirely.
                %i still need to figure out how to define
                %lastnonzeroentry and filterIndex
                fCell = [fCell, {filterName; f}];
            end           
        end
        fNames = fCell(1,:);
    else
        for o = 1:numFilters
            filterName = filters{o};
            filterType = regexp(filterName,'\D*(?=_\d)', 'match');
            filterType = filterType{1};
            sigmas = regexp(filterName,'(?<=_)\d(?=.)','match'); 
            if ~strcmp(filterName, 'original') && ~strcmp(filterName, 'class')
                f = filterImage(im, filterType, sigmas);    
                ind = strcmp(fNames, filterName);
                fCell{end+1, ind} = f; 
            end
        end
    end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%  
%     fCellSmall =  fCell(~cellfun('isempty',fCell))
    fCellLength1 = size(fCell,2);
    fCellLength2 = size(fCell,1);
    for p = 1:fCellLength1
        arff{1, p} = fCell{1, p};
    end
%     n = 2;
    
   n = 1;
   arfftemp = [];
   stackSize = length(fCell{q,p}(:));
   for p = 1:fCellLength1
       for q = 2:fCellLength2
            stack = fCell{q,p};
            lastNonZeroEntry = 1;
%                      if ~isdempty(stack)
%                 for k = 1:rows
%                     for m = 1:cols
%                         arff{end+1:end+1+rows+cols+zSlices,p} = stack(:); 
%                            arff{:,p} = [arff{:,p}, stack(:)'];
            if ~isempty(stack)
%                 arfftemp(:,p) = vertcat(arfftemp,stack(:)'); 
%                 arfftemp = vertcat(arfftemp,stack(:)'); 
                n = n + 1;
                arfftemp(lastnonzeroentry:stackSize,p) = stack(:)';
            end
%                           arff{end+1,p} = stack{q}(k,m);
%                         n = n + 1;
%                     end
%                 end
%             end
       end
   end
csvwrite('C:\Users\ArmandoReimer\Desktop\arffout.arff', arfftemp);
% T = cell2table(arffcell, 'VariableNames', arff);
% writetable(T,[DropboxFolder,filesep,Prefix,filesep, Prefix,'.arff'])