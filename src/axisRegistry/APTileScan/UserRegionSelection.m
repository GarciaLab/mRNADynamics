% UserRegionSelection.m
% author: Gabriella Martini
% date created: 7/30/20
% date last modified: 8/12/20


function [MinRow, MaxRow, MinColumn, MaxColumn] = UserRegionSelection(normed_scores,scores,areas,...
    tile_array, tArRange, tAcRange, tA_ind, tB_ind)
close all
MinColumn = 1;
MaxColumn = size(normed_scores, 2);
MinRow = 1;
MaxRow = size(normed_scores, 1);
tileA = tile_array.tiles{tA_ind};
tileB = tile_array.tiles{tB_ind};
zdim_scores = size(scores, 3);
TileImages = figure(1);
subplot(3,2,1);
imagesc(tileA);
title(['Tile A: ', num2str(tA_ind)])
subplot(3,2,2);
imagesc(tileB);
title(['Tile B: ', num2str(tB_ind)])
tilearrayAx=subplot(3,2,[3,4,5,6]);
ScoresImageFig = figure(2);
scAx = axes(ScoresImageFig);



%Now, do the correction
cc=1;
% imagesc(tileA)

% figure(2) 
% imagesc(tileB)


while (cc~='x')
    scores_copy = normed_scores;
    scores_mask = zeros(size(normed_scores));
    scores_mask(MinRow:MaxRow, MinColumn:MaxColumn) = 1;
    scores_copy(scores_mask == 0) = max(max(scores_copy));
    [minr, minc] = find(scores_copy == min(min(scores_copy)));
    tile_array_copy = tile_array;
    newr = tArRange(minr);
    newc = tAcRange(minc);
    tile_array_copy.rows{tA_ind} = newr;
    tile_array_copy.cols{tA_ind} = newc;
    rmins = [tile_array_copy.rows{:}];
    top_limit = min(rmins);
    for rr =1:length(tile_array_copy.rows)
        tile_array_copy.rows{rr} = tile_array_copy.rows{rr} + (1-top_limit);
    end

    cmins = [tile_array_copy.cols{:}]; 
    left_limit = min(cmins);
    for cc =1:length(tile_array_copy.cols)
        tile_array_copy.cols{cc} = tile_array_copy.cols{cc} + (1-left_limit);
    end
    imagesc(imstitchTile(tile_array_copy), 'Parent', tilearrayAx)
    imagesc(normed_scores, 'Parent', scAx)
    %imshow(APImage,DisplayRange)
    
    axis image
    %axis off
    title('Boundaries show in red and minimum score');
        
    hold on
    
    
    scatter(minc, minr, 100, 'r.') 
    colorbar
   
    xline(MinColumn,'r')

    xline(MaxColumn,'r')
    yline(MinRow,'r')

    yline(MaxRow,'r')


    
    hold off
    
    counter = 3;
    for i=1:zdim_scores
        figure(counter)
        imagesc(scores(:,:,i)./areas(:,:,i))
        counter = counter + 1;
    end
    
    figure(ScoresImageFig)
    ct=waitforbuttonpress;
    cc=get(ScoresImageFig,'currentcharacter');
    cm=get(scAx,'CurrentPoint');
    
    
    if (ct~=0)&(cc=='c')        %Clear all AP information
        MinColumn = 1;
        MaxColumn = size(normed_scores, 2);
        MinRow = 1;
        MaxRow = size(normed_scores, 1);
    elseif (ct~=0)&(cc=='l')	%Select anterior end
        [MinColumn,coordLy]=ginputc(1,'Color',[1,1,1]);
    elseif (ct~=0)&(cc=='r')    %Select posterior end
        [MaxColumn,coordRy]=ginputc(1,'Color',[1,1,1]);
    elseif (ct~=0)&(cc=='t')    %Select posterior end
        [coordTx,MinRow]=ginputc(1,'Color',[1,1,1]);
    elseif (ct~=0)&(cc=='b')    %Select posterior end
        [coordBx,MaxRow]=ginputc(1,'Color',[1,1,1]);
    end
    
    MaxRow = uint16(round(MaxRow));
    MinRow = uint16(round(MinRow));
    MaxColumn = uint16(round(MaxColumn));
    MinColumn = uint16(round(MinColumn));
end

close all

end
