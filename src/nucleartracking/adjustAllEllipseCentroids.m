function Ellipses = adjustAllEllipseCentroids(Prefix, varargin)


%function description- we want to use existing ellipses, but mildly adjust the position and
%diameter to better model the nuclei

nWorkers = 1;
displayFigures = false;

for i = 1:(numel(varargin)-1)
    if i ~= numel(varargin)
        if ~ischar(varargin{i+1})
            eval([varargin{i} '=varargin{i+1};']);
        end
    end
end

[~,~,DropboxFolder,~,PreProcPath]=...
    DetermineLocalFolders(Prefix);

load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'], 'FrameInfo');

[xDim, yDim, ~, ~, ~,...
    nFrames, ~, ~] = getFrameInfoParams(FrameInfo);

load([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'], 'Ellipses');

hisFile = [PreProcPath, filesep, Prefix, filesep, Prefix, '_hisMat.mat'];
hisMat = loadHisMat(hisFile);

hisMat = double(hisMat);

d = struct;

EllipsesNew = Ellipses;

for f = 1:nFrames
    
    d(f).hisImage = squeeze(hisMat(:, :, f));
        
    rad0 = median(Ellipses{f}(:, 3), 1);  %3 and 4 are the semimajor/minor axes of the ellipses
    
    d(f).ellipseMask0 = double(makeNuclearMask(Ellipses{f}, [yDim, xDim])); %make a mask with the initial ellipse configuration
    
    smoothSigma = rad0/2;
    d(f).hisSmooth = imgaussfilt(d(f).hisImage, smoothSigma, 'Padding',0); %denoise a little
    
    d(f).smoothMasked = d(f).hisSmooth.*d(f).ellipseMask0;
    
    d(f).hisMasked = d(f).hisImage.*d(f).ellipseMask0;
    
    d(f).label0 = bwlabel(d(f).ellipseMask0);
    
    nRegions = max(max(d(f).label0));
        
    for r = 1:nRegions
        
        maskTemp = (d(f).label0 == r) .* d(f).hisSmooth;
        hold on;
        [~,ind] = max(maskTemp,[],'all','linear');
        [y1, x1] = ind2sub([yDim, xDim],ind);
        
       if displayFigures
            figure(1);
            imagesc(maskTemp, [min(min(maskTemp(maskTemp>0))), max(maskTemp(:))]);
            hold on;
            plot(y1, x1, 'xk', 'LineWidth', 2);
            hold off;
       end
        
        EllipsesNew{f}(r, 1) = x1;
        EllipsesNew{f}(r, 2) = y1;
        
    end
    
    d(f).ellipseNewMask = double(makeNuclearMask(EllipsesNew{f}, [yDim, xDim])); %make a mask with the initial ellipse configuration
    
    d(f).smoothNewMasked = d(f).hisSmooth.*d(f).ellipseNewMask;
       
    d(f).hisNewMasked = d(f).hisImage.*d(f).ellipseNewMask;

    
    if displayFigures
        ims = fieldnames(d);
        imageCell = {};
        for i = 1:length(ims)
            imageCell{i} = d(f).(ims{i});
        end
        imageTile(imageCell, ims);
    end
    
    
end

Ellipses = EllipsesNew;

save([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'], 'Ellipses', '-v6');

end