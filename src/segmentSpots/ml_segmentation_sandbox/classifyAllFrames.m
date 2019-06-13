function classifyAllFrames(Prefix, trainingFeatures, sigmaVec, beta)

[~,ProcPath,DropboxFolder,~, PreProcPath,...
    ~, Prefix, ~,~,~,~,...
    ~, spotChannels] = readMovieDatabase(Prefix);

if length(spotChannels) > 1
    error('nope. talk to nick');
end
% determine size of image stack
load([DropboxFolder, filesep, Prefix, filesep, 'FrameInfo.mat']);
zDim = FrameInfo(1).NumberSlices + 2;
yDim = FrameInfo(1).LinesPerFrame;
xDim = FrameInfo(1).PixelsPerLine;
pixelSize = FrameInfo(1).PixelSize / 1000;
nFrames = length(FrameInfo);
zMax = zDim - 2;

nameSuffix = ['_ch', iIndex(spotChannels, 2)];

filterImageFrames(Prefix,sigmaVec,trainingFeatures, beta);
% % label full stack
% pihat = mnrval(beta,tData);
% probStack = reshape(pihat(:,2),[yDim, xDim, zMax, nFrames]);
% 
% for currentFrame = 1:nFrames
%     for i = 1:zMax
%         p_name = ['prob', Prefix, '_', iIndex(currentFrame, 3), '_z', iIndex(i, 2), nameSuffix, '.tif'];
%         imwrite(uint16(probStack(:, :, i, currentFrame)), [ProcPath, filesep,Prefix, '_', filesep, p_name]);
%     end
% end

end
