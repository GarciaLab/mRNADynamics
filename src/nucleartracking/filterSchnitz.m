function schnitzcells = filterSchnitz(schnitzcells, imSize)

% Disapprove spurious nuclei that are too near edges
% and don't last long.

%% time thresholding
lengths = {[], [], []};

for nc = 12:14
    
    snc = schnitzcells([schnitzcells.cycle] == nc);
    for k = 1:length(snc)
        lengths{nc-11} = [lengths{nc-11}, length(snc(k).frames)];
    end
    
    %if a nucleus exists for less than half the median lifetime, disapprove it.
    frameCutoff(nc-11) = .5*median(lengths{nc-11}); 
    
    
end


%% size thresholding
ymax = imSize(1);
xmax = imSize(2);

thresh = .5; %     50% of the nuclear radius


for s = 1:length(schnitzcells)
    
    schnitz = schnitzcells(s);
    
    for frameIndex = 1:length(schnitz.frames)
        
        radius  = double(schnitz.len(frameIndex));
        
        if ( (schnitz.cenx(frameIndex) + radius*thresh ) > xmax ||...
                (schnitz.ceny(frameIndex) + radius*thresh ) > ymax ||...
                (schnitz.cenx(frameIndex) - radius*thresh  ) < 0 ||...
                (schnitz.ceny(frameIndex) - radius*thresh  ) < 0 )
            
            schnitzcells(s).FrameApproved(frameIndex) = false;
        else
            schnitzcells(s).FrameApproved(frameIndex) = true;
        end
        
    end
    
    for nc = 12:14
        if schnitzcells(s).cycle == nc &&...
                length(schnitzcells(s).frames) < frameCutoff(nc-11)
            schnitzcells(s).Approved = false;
        else
            schnitzcells(s).Approved = true;
        end
    end
    
    %disapprove an entire nucleus if less than half the frames are good.
    if sum(schnitzcells(s).FrameApproved) < .5*length(schnitzcells(s).frames)
        schnitzcells(s).Approved = false;
    else
        schnitzcells(s).Approved = true;
    end
    
end

end