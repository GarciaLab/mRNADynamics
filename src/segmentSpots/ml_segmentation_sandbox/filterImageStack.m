% Function to apply selected filters to specified frame of a time series
function featureTable = filterImageStack(Prefix,frame,sigmaVec,featureCell,varargin)
trainingFlag = 0;
featureTable = [];
for i = 1:numel(varargin)
    if strcmpi(varargin{i},'Training')
        trainingFlag = 1;
        featureTable = varargin{i+1};
    end
end
[SourcePath,ProcPath,DropboxFolder,MS2CodePath, PreProcPath,...
    Folder, Prefix, ExperimentType,Channel1,Channel2,OutputFolder,...
    Channel3] = readMovieDatabase(Prefix);

spotChannels = find(contains([Channel1,Channel2,Channel3],'CP'));
if length(spotChannels) > 1
    error('nope. talk to nick');
end

% determine size of image stack
load([DropboxFolder, filesep, Prefix, filesep, 'FrameInfo.mat']);
zDim = FrameInfo(1).NumberSlices + 2;
yDim = FrameInfo(1).LinesPerFrame;
xDim = FrameInfo(1).PixelsPerLine;
pixelSize = FrameInfo(1).PixelSize;

gp = gpuDevice;
if gp.AvailableMemory < 1E9
    gpuDevice(1);
end
% load stack
rawStack = zeros(yDim,xDim,zDim, 'single','gpuArray');
mcpChannel = find(contains([Channel1,Channel2,Channel3],'MCP'));
for z = 1:zDim
    fileName = [Prefix '_' sprintf('%03d',frame) '_z' sprintf('%02d',z)  '_ch' sprintf('%02d',mcpChannel) '.tif'];
    rawStack(:,:,z) = imread([PreProcPath '/' Prefix '/' fileName]);
end
rawStack = gpuArray(rawStack);
% apply filters
for i = 1:numel(featureCell)
    feature = featureCell{i};
    for j = 1:numel(sigmaVec)
        if strcmpi(feature,'Difference_of_Gaussian')
            ft_im = filterImage(rawStack, 'Difference_of_Gaussian', {sigmaVec(j), sigmaVec(j)*4});
        else
            ft_im = filterImage(rawStack, feature, {sigmaVec(j)});
        end
        
        ft_im = gather(ft_im);
        featureName = [feature '_s' num2str(sigmaVec(j))];
        if trainingFlag
            % initialize field if it's not a;ready present
            if ~ismember(featureName,featureTable.Properties.VariableNames)
                featureTable.(featureName) = NaN(size(featureTable,1),1);
            end
            % add filtered values
            frame_vec = featureTable.frame';
            z_vec = featureTable.z';
            index_vec = featureTable.linIndex';
            for z = unique(z_vec(frame_vec==frame))
                slice = ft_im(:,:,z);
                featureTable(frame_vec==frame&z_vec==z,:).(featureName) = reshape(slice(index_vec(z_vec==z&frame_vec==frame)),[],1);
            end
        else
            if i == 1 && j == 1
                featureTable = ft_im(:);
            else
                featureTable = [featureTable ft_im(:)];
            end
        end
    end
end
end