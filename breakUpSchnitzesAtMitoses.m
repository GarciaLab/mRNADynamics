function schnitzcells = breakUpSchnitzesAtMitoses(schnitzcells, ncs, nFrames)

% ncs = [zeros(1,8), nc9, nc10, nc11, nc12, nc13, nc14];
nFrames = 500;
cycleFrames = nan(1,500);

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
    hasncs = [];
    
    for nc = 1:length(ncs)
        hasnc = nc*logical(sum(sc.frames' >= ncs(nc) & sc.frames' < ncs(nc) + 1));
        if hasnc
            hasncs = [hasncs, hasnc]
        end
    end
    
    hasncs = unique(cycleFrames(sc.frames))
    
    if length(hasncs) > 1
        
        keyboard
        for i = 1:length(hasncs)
%             cycleFrames = find(sc.frames <  hasncs(i) & sc.frames > ncs(hasncs(i)-1));
            newInd = nNuclei + j;
            tempSchnitzcells(newInd) = sc;
            newFrames = (cycleFrames(sc.frames) == i);
            tempSchnitzcells(newInd).frames = sc.frames(newFrames);
            tempSchnitzcells(newInd).cenx = sc.cenx(newFrames);
            tempSchnitzcells(newInd).ceny = sc.ceny(newFrames);
            tempSchnitzcells(newInd).len = sc.len(newFrames);
            tempSchnitzcells(newInd).cellno = sc.cellno(newFrames);
            tempSchnitzcells(newInd).Fluo = sc.Fluo(newFrames, :);
            tempSchnitzcells(newInd).cellno = sc.cellno(newFrames);
            tempSchnitzcells(newInd).APPos = sc.APPos(newFrames);
            tempSchnitzcells(newInd).DVPos = sc.DVPos(newFrames);
             tempSchnitzcells(newInd).FrameApproved = sc.FrameApproved(newFrames);
             try
             tempSchnitzcells(newInd).FluoTimeTrace = sc.FluoTimeTrace(newFrames);
             end
            j = j+1;
        end
        
        
    end
    
end

end