function schnitzcells = filterSchnitz(schnitzcells, imSize)


lengths = {[], [], []};

for nc = 12:14
    
    snc = schnitzcells([schnitzcells.cycle] == nc)
    for i = 1:length(snc)
        lengths{nc-11} = [lengths{nc-11}, length(snc(i).frames)];
    end
    
    frameCutoff(nc-11) = .5*median(lengths{nc-11}); %if a nucleus exists for less than half the median lifetime, disapprove it.
    
    
end

xmax = imSize(1);
ymax = imSize(2);

thresh = .75;


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
    
    if sum(schnitzcells(s).FrameApproved) < .5*length(schnitzcells(s).frames)
        schnitzcells(s).Approved = false;
    end
    
end

end