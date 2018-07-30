function [ x_shift, y_shift ] = measure_shift( frame1, frame2, varargin )
%MEASURE_SHIFT This function tries to estimate the shift between two
%frames.
%   Algorithm tries to detect 'bulk shifts' from one frame to the other and
%   to return its x and y magnitude.
%   It computes the intercorrelation between the two frames and then detect
%   the maximum of the output. The intercorrelation is just calculated for
%   a limited amount of frames, based on the assumption that the shift does
%   not exceed a radius of RADIUS pixels.
%
%   Inputs :
%   frame1 : frame before shift occurs
%   frame2 : frame after shift occurred
%   RADIUS : optional parameter defining the maximal possible shift.
%   Default is 8.
%
%   Outputs :
%   x_shift : shift along vertical axis of image(?!). Positive values mean a shift to
%   the bottom occurred.
%   y_shift : shift along horipadarray(frame1,2,'both'))zontal axis image(?!). Positive values mean a shift to
%   the right occurred.

if nargin > 2
    RADIUS = varargin{1};
else
    RADIUS = 15;
end;

intC = interCorrelation(double(frame1), double(rot90(frame2,2)), RADIUS);

% figure;
% imagesc(intC);

[dummy,i] = max(intC(:));
[x_shift, y_shift] = ind2sub(size(intC),i);
% x_shift = -x_shift + RADIUS + 1;
% y_shift = -y_shift + RADIUS + 1;

x_shift = -x_shift + size(intC,1)/2 ;
y_shift = -y_shift + size(intC,2)/2 ;

function out = interCorrelation(frame1, frame2, radius)
% Calculate the intercorrelation in a square of side 2*'radius'+1 aroung
% origins.
PADDING = 10;
out = ifftshift(abs(ifft2(fft2(padarray(frame1,[PADDING PADDING],'both')).*fft2(padarray(frame2,[PADDING PADDING],'both')),'symmetric')));%imfilter(frame1,frame2);%zeros(2*radius+1);

out = out(1+PADDING:end-PADDING,1+PADDING:end-PADDING);
% 
% center_correction = radius+1;
% 
% [nx, ny] = size(frame1);
% 
% for x = 1:2*center_correction - 1
%     
%     shiftX = x - center_correction;
%     lower_x1 = max(1, 1 + shiftX);
%     lower_x2 = max(1, 1 - shiftX);
%     upper_x1 = min(nx, nx + shiftX);
%     upper_x2 = min(nx, nx - shiftX);
%     
%     for y = 1:2*center_correction - 1
%         
%         shiftY = y - center_correction;
%         
%         if shiftX^2+shiftY^2 > RADIUS^2
%             continue
%         end
%         
%         lower_y1 = max(1, 1 + shiftY);
%         lower_y2 = max(1, 1 - shiftY);
%         upper_y1 = min(ny, ny + shiftY);
%         upper_y2 = min(ny, ny - shiftY);
% 
%     
%         out(x,y) = sum(sum(frame1(lower_x1:upper_x1,lower_y1:upper_y1).*frame2(lower_x2:upper_x2,lower_y2:upper_y2)));
%         
%     end
% end
% 
end


end

