function I_mask=GetEmbryoMaskWithThreshold(I,blurSigma, threshold)
    % Presumably, onbe of the corners of the image is outside of the embryo.
%     corners = [max(flatten(I(1:10,1:10))), ...
%                max(flatten(I(end-9:end,1:10))), ...
%                max(flatten(I(1:10,end-9:end))), ...
%                max(flatten(I(end-9:end,end-9:end)))];
%     minCorner = min(corners)+10;
%     % disregard the corners that are inside the embryo
%     corners(corners > 10* minCorner) = minCorner; 
%     threshold = 1.3 * (max(corners)+40);
if ~exist('threshold', 'var') || isempty(threshold)
    if max(I(:))>2^12
        % probably 16 bit
        threshold = 2000;
    else
        threshold = 50; %200
    end
end
bwfill=imfill(I>threshold,'holes');
I_inside = bwselect(bwfill,round(size(bwfill,2)/2),round(size(bwfill,1)/2));
I_inside=uint16((2^16-1)*I_inside);
I_blurred = imfilter(I_inside,...
    fspecial('gaussian',2*blurSigma,blurSigma),'symmetric','conv');
level = graythresh(I_blurred);
I_mask = im2bw(I_blurred,level);
end
