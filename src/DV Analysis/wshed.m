function maskOut = wshed(maskIn, varargin)

arguments
    maskIn (:,:) logical
end
arguments (Repeating)
    varargin
end

%wrapper function for watersheeding. main input parameter would be
%distThresh;

displayFigures = false;
distThresh = 3;

for i = 1:(numel(varargin)-1)
    if i ~= numel(varargin)
        if ~ischar(varargin{i+1})
            eval([varargin{i} '=varargin{i+1};']);
        end
    end
end

maskInOriginal = maskIn;
%invert the image polarity if necessary. not sure if this is the best
%metric. size or euler number could also work
% stats =  regionprops(maskIn, 'EulerNumber');
% statsInvert =  regionprops(~maskIn, 'EulerNumber');
% maskIn = gpuArray(maskIn);
if chooseKLabel(maskIn)
    maskIn=~maskIn;
end

% if min([statsInvert.EulerNumber]) < min([stats.EulerNumber])
%     maskIn=~maskIn;
% end

%get distance transform
D = -bwdist(maskIn);

%quality control
D = imimposemin(D, imextendedmin( D, distThresh ));

%perform the actual watershed
maskOut = ~~watershed(D);
maskOut(maskIn) = 0;
% 
% if ~chooseKLabel(maskOut)
%     maskOut=~maskOut;
% end

if displayFigures
    imshowpair(maskInOriginal, maskOut, 'montage');    
end

end