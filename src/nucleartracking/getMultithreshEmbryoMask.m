function [embryoMask]=getMultithreshEmbryoMask(FrameInfo, names, diameters)
%%
lastframe = length(names);
nucleusDiameter = diameters(end);
I = imread(names{lastframe});
f_sigma = round(nucleusDiameter / FrameInfo(1).PixelSize);
I_blurred = imfilter(I,...
     fspecial('gaussian',2*f_sigma,f_sigma),'symmetric','conv');
embryoMask = true(size(imread(names{lastframe})));
goodmask = 0;
threshcount = 0;
while ~goodmask
    [xm, ym] = find(embryoMask == 1);
    b_points = boundary(xm,ym); 
    figure(1)
    imagesc(I)
    hold on 
    ax = gca();
    scatter(ax, ym(b_points),xm(b_points),60,'r','filled');
    title('Sample Frame with embryo mask indicated in red')
    prompt = ['Boundary of the current embryo mask indicated in red.\n',...
    'Does this mask properly exclude the background in the image\n',...
    'while keeping all nuclei (y/n)? '];
    ID = input(prompt,'s');
   
    if ID == 'y'       
        goodmask = 1;
        close all
        break
    elseif ID == 'n'
        threshcount = threshcount + 1;
    else
        error('Must indicate "y" or "n"')
    end
    if threshcount == 1
        embryoMask = imbinarize(I_blurred);
    elseif threshcount > 1
        levels = multithresh(I_blurred,threshcount);
        embryoMask = I_blurred > levels(1);
    end
    close all
end
%     
end