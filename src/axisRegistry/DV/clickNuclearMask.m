function clickNuclearMask ()

surfmax = imread('E:\Armando\LivemRNA\Data\Dropbox\DorsalSyntheticsDropbox\2019-07-18-1Dg_Zeiss\DV\surf_max.tif');
surfmaxmch= imread('E:\Armando\LivemRNA\Data\Dropbox\DorsalSyntheticsDropbox\2019-07-18-1Dg_Zeiss\DV\surfmaxmCh.tif');
surfmaxvenus= imread('E:\Armando\LivemRNA\Data\Dropbox\DorsalSyntheticsDropbox\2019-07-18-1Dg_Zeiss\DV\surfmaxVenus.tif');

embmask = imread('E:\Armando\LivemRNA\Data\Dropbox\DorsalSyntheticsDropbox\2019-07-18-1Dg_Zeiss\APDetection\midMask.tif');
embmask2 = getEmbryoMaskLive(surfmaxvenus, .5);
SE = strel('disk',20);
a = imerode(embmask2, SE);
figure()
im = imcomplement(surfmaxmch).*uint16(a);
minv = min(im(im ~= 0));


rad = 4;

Overlay = figure;
ax = axes(Overlay);
% z = zoom(Overlay);
% z.ButtonDownFilter = @mycallback;
% z.ButtonDownFilter = true;

imshow(im, [minv, max(im(:))], 'Parent', ax)

%%make corrections
cc=1;

centers = [];
radii = [];
while (cc~='x')

    ct=waitforbuttonpress;
    cc=get(Overlay,'currentcharacter');
    cm=get(ax,'CurrentPoint');

    if (ct==0)&(strcmp(get(Overlay,'SelectionType'),'normal'))
        cc=1;
        if cm(1,2)>0 & cm(1,1)>0 & cm(1,2) <= size(im, 1) & cm(1,1) <= size(im, 2)
            centers = vertcat(centers, [cm(1, 1), cm(1, 2)]);
            radii = vertcat(radii, rad);
        end
        imshow(im, [minv, max(im(:))], 'Parent', ax)
        viscircles(ax, centers, radii,'Color','r', 'LineWidth', .2);

    elseif (ct==0)&(strcmp(get(Overlay,'SelectionType'),'alt'))
        cc=1;
        point = [cm(1,1), cm(1,2)];
        if (cm(1,2)>0)&(cm(1,1)>0)&(cm(1,2)<=size(im, 1))& (cm(1,1)<=size(im, 2))
            dif = centers - point;
            [~, nearestNuc] = min(abs(dif));
            centers(nearestNuc(2), :) = [];
            radii(nearestNuc(2), :) = [];
            imshow(im, [minv, max(im(:))], 'Parent', ax)
            viscircles(ax, centers, radii,'Color','r', 'LineWidth', .2);
        end
%     elseif  (ct~=0)&(cc==' ')
% %         if strcmp(z.Enable,'on')
% %             z.Direction = 'in';
% %       if strcmp(z.Enable,'off')
% %           z.Enable = 'on';
% zoom on
% %       end
    elseif (ct~=0)&(cc=='9')    %Debug mode
        keyboard
    end
end



end

% function [res] = mycallback(obj,event_obj)
% % obj          handle to the object clicked on
% % event_obj    struct for event data (empty in this release)
% % res [output] a logical flag determines whether the zoom
% %              operation should take place(for 'res' set
% %              to 'false' or the 'ButtonDownFcn' property
% %              of the object should take precedence (when
% %              'res' is 'true')
% res = false;
% end