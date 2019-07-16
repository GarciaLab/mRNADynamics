function [schnitzcells, Ellipses] = breakUpSchnitzesAtMitoses(schnitzcells, Ellipses, ncs, nFrames)

cycleFrames = nan(1,nFrames);

for i = 1:length(ncs)
    if i==14
        cycleFrames(ncs(i):end) = 14;
    elseif ncs(i) == 0
        cycleFrames(1:ncs(i+1)) = i;
    else
        cycleFrames(ncs(i):ncs(i+1)) = i;
    end
end


tempSchnitzcells = schnitzcells;
nNuclei = length(schnitzcells)
j = 1;
for s = 1:nNuclei
    
    sc  = schnitzcells(s);  
    
    hasncs = unique(cycleFrames(sc.frames));
    
    if length(hasncs) > 1
        
%         keyboard
        for i = 1:length(hasncs)
            newInd = nNuclei + j;
            tempSchnitzcells(newInd) = sc;
            newFrames = cycleFrames(sc.frames) == hasncs(i)

            tempSchnitzcells(newInd).frames = sc.frames(newFrames)
            tempSchnitzcells(newInd).cenx = sc.cenx(newFrames);
            tempSchnitzcells(newInd).ceny = sc.ceny(newFrames);
            tempSchnitzcells(newInd).len = sc.len(newFrames);
            tempSchnitzcells(newInd).cellno = sc.cellno(newFrames);
            tempSchnitzcells(newInd).Fluo = sc.Fluo(newFrames, :);
            tempSchnitzcells(newInd).cellno = sc.cellno(newFrames);
            tempSchnitzcells(newInd).APpos = sc.APpos(newFrames);
            tempSchnitzcells(newInd).DVpos = sc.DVpos(newFrames);
             tempSchnitzcells(newInd).FrameApproved = sc.FrameApproved(newFrames);
             try
                tempSchnitzcells(newInd).FluoTimeTrace = sc.FluoTimeTrace(newFrames);
             end
            j = j+1;
        end
        
        tempSchnitzcells(s) = [];
        j = j-1;
        
        
    end
    
end

schnitzcells = tempSchnitzcells;

Ellipses = addSchnitzIndexToEllipses(Ellipses, schnitzcells)

end