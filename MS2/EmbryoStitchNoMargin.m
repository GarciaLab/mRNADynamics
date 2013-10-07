function [ImageOutput,xo1,yo1]=EmbryoStitchNoMargin(leftfilename,rightfilename,FFfilename,xo1,yo1)

    %HG: Modified EmbryoStitch in order to get rid of the issues at the
    %margins.
    
    
%Parameters:
Margin=10;      %Pixels to get rid of at the left and margins
    

   % Read in flat field, filter, convert from uint16 to double
   FF=imread(FFfilename);
   FF=imfilter(double(FF), fspecial('disk', 30), 'replicate', 'same');
   FF=FF/mean(FF(:));
   FF=(FF-1)*1+1;
   
   
   %See if we have a histone channel. It assumes that it's the second one
   %if so.
   ImInfo=imfinfo(leftfilename);
   
   if length(ImInfo)==2
       % Read in left and right images
       left=imread(leftfilename,2);
       right=imread(rightfilename,2);
   else
       % Read in left and right images
       left=imread(leftfilename);
       right=imread(rightfilename);
   end
   
   % x and y limits
   x_min=280; x_max=450; %y_min=-40; y_max=40;
   y_min=-50; y_max=50;
   
   
%Crop the images
FF=FF(Margin+1:end-Margin,Margin+1:end-Margin);
left=left(Margin+1:end-Margin,Margin+1:end-Margin);
right=right(Margin+1:end-Margin,Margin+1:end-Margin);
   

   % loop for what to do if there is no xo1 ??? or maybe solves for xo1 yo1
if ~exist('xo1') || isempty(xo1)
        
     [Dummy, xo1 yo1]=autoStitch(left,right,y_min:1:y_max,x_min:1:x_max,1,[1 2]); 
     
     if (xo1-x_min)<2 || (x_max-xo1)<2 || sum([yo1-y_min y_max-yo1]<2)>=1
      %[Dummy, xo1 yo1]=autoStitch(left,right,-50:1:50,10:1:511,1,[1 2]);   
     end

    else
     [Dummy, xo1 yo1]=autoStitch(left,right,y_min:1:y_max,x_min:1:x_max,1,[1 2]);
     
     if (xo1-x_min)<2 || (x_max-xo1)<2 || sum([yo1-y_min y_max-yo1]<2)>=1
      [Dummy, xo1 yo1]=autoStitch(left,right,-50:1:50,10:1:511,1,[1 2]); 
     end
end
% end xo1 loop
    
%Then, do this    
    left=double(left)./FF;
    right=double(right)./FF; 
    im2=imstitch(left,right,xo1,yo1,[1 2]);
    [h w]=size(left);
    imwhole=im2(:,1:2*w+1-xo1);
    
    %cropping
    cr=10; % amount of right picture to crop away


    
    imwhole(y_max:h-y_max,w-xo1:w-xo1+10) = left(y_max+yo1:h-y_max+yo1,w-xo1:w-xo1+cr);
    %imwhole=im2; 
    
    
    %Output
%     filename='output.tif';
%     imwrite(uint16(imwhole),filename,'compression','none');
    ImageOutput=uint16(imwhole);
   
    
    %Account for the Margin in the shift that is given
    xo1=xo1+Margin*2;
    
end

function [imm2 xo1 yo1] = autoStitch(varargin)

    if size(varargin,2)==6
        left = varargin{1};
        right = varargin{2};
        YRange = varargin{3};
        XRange = varargin{4};
        fig = varargin{5};
        orientation = varargin{6};
    else
        warning('Number of inputs must be 6');
    end

        %% calculate normalized images
        %??????
    leftN = normalizer(left);
    rightN = normalizer(right);

        %% calculate offsets by looking for the row where the embryo is same width
    if size(varargin,2)==6
        [xo1, yo1] = matchManual(leftN, rightN, YRange, XRange, fig);
    end
    
        %% stitch
    imm1 = imstitch(left,right, xo1, yo1,orientation);
    [h, w] = size(left);
     imm2=imm1(:,1:2*w+1-xo1);
    %  imm2=imm1;
    
end

%Works
function [xoffset, yoffset] = matchManual(leftN, rightN, yRange, xRange, fig) 
    % goes through where the widths are close, and match the two normalized
    % images to see how good they match

    ys = repmat(yRange, 1, length(xRange));
    xs = reshape(repmat(xRange, length(yRange), 1), 1, []);
    
    % claculate scores in region and find the best
    scores = arrayfun(@(x, y) picDiff(leftN, rightN, x, y), xs, ys);
    best = find(scores==min(scores), 1, 'first');
    
    if fig>1,
        figure('windowstyle', 'docked');
        hold on;
%         plot(ys, scores, 'b');
%         plot(ys(best), scores(best), 'r.');
        plot(scores, 'b');
        plot(best, scores(best), 'r.');
    end
    
    xoffset = xs(best);
    yoffset = ys(best);
end

%Works
function score = picDiff(leftN, rightN, x, y)
    % returns how good the two normalized images match.
    [h, w] = size(leftN);
    if y>=0,
        leftWindow = leftN(y+1:h,w-x+1:w);
        rightWindow = rightN(1:h-y, 1:x);
        area = (h-y)*x;
    else
        leftWindow = leftN(1:h+y,w-x+1:w);
        rightWindow = rightN(-y+1:h,1:x);
        area = (h+y)*x;
    end
    
%         if x>=0,
%         leftWindow = leftN(h-y+1:h, x+1:w);
%         rightWindow = rightN(1:y, 1:w-x);
%         area = (w-x)*y;
%     else
%         leftWindow = leftN(h-y+1:h, 1:w+x);
%         rightWindow = rightN(1:y, -x+1:w);
%         area = (w+x)*y;
%     end
    
    diffSquared = (leftWindow - rightWindow) .^ 2;
    score = sum(diffSquared(:))/area;
end

%Works
function im2 = normalizer(im, width) 
    if ~exist('width', 'var'), width = 80; end
    im = double(im);
    im = im - imfilter(im, fspecial('average', width)); %mean normalize
    im = im ./ (sqrt(imfilter(im .* im, fspecial('average', width))));
    im = im / std(im(:));
    im2 = double(im); % > mean(im(:)) + 2*std(im(:)));
end


function imm2 = imstitch(left,right, xo1, yo1,order)
    %% stitches images given offsets. you specify the drawing order
    % [xo1 yo1] = offsets
    % order = draw order {1=left, 2 = right}. drawn left to right
    
    [h, w] = size(left);
    imm2 = zeros(h, 2*w);

    % ?????
        for i=order
            switch i
                case 1
                    imm2(:, 1:w) = left;
                case 2
                    imm2(:, w+1-xo1:2*w-xo1) = circshift(right(1:h, :), [yo1, 0]);
            end
        end
    %end
  %if abs(xo1)>abs(xo2)
      yshift=yo1;
  %else
  %    xshift=xo2;
  %end
  imm2=circshift(imm2(:,:),[-yshift, 0]);
end