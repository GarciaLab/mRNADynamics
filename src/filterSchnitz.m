function schnitzcells = filterSchnitz(schnitzcells, imSize)

%armando's criteria for good nuclei to include in dorsal analysis. 

%% time thresholding
lengths = {[], [], []};

for nc = 12:14
    
    snc = schnitzcells([schnitzcells.cycle] == nc);
    for i = 1:length(snc)
        lengths{nc-11} = [lengths{nc-11}, length(snc(i).frames)];
    end
    
    frameCutoff(nc-11) = .5*median(lengths{nc-11}); %if a nucleus exists for less than half the median lifetime, disapprove it.
    
    
end


%% size thresholding
xmax = imSize(1);
ymax = imSize(2);

thresh = .75; % 75% of the nuclear radius


for s = 1:length(schnitzcells)
    
    schnitz = schnitzcells(s);
    
    for f = 1:length(schnitz.frames)
        
        r  = schnitz.len(f);
        
        if ( (schnitz.cenx(f) + r*thresh ) > xmax |...
                (schnitz.ceny(f) + r*thresh ) > ymax |...
                (schnitz.cenx(f) - r*thresh  ) < 0 |...
                (schnitz.ceny(f) - r*thresh  ) < 0 )
            
            schnitzcells(s).FrameApproved(f) = false;
            
        end
        
    end
    
    for nc = 12:14
        if schnitzcells(s).cycle == nc & length(schnitzcells(s).frames) < frameCutoff(nc-11)
            schnitzcells(s).Approved = false;
        end
    end
    
    %disapprove an entire nucleus if less than half the frames are good. 
    if sum(schnitzcells(s).FrameApproved) < .5*length(schnitzcells(s).frames)
        schnitzcells(s).Approved = false;
    end
    
end

end