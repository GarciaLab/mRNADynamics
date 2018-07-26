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
